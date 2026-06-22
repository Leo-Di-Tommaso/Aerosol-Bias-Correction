# =============================================================================
# MethaneSAT XCO2 L3 — Aerosol Bias Correction Pipeline
# =============================================================================
#
# Structure
# ---------
#   SECTION 0  Libraries and constants
#   SECTION 1  User settings
#   SECTION 2  Helper functions
#              A. L3 file discovery (from plotting_xco2_l3_leo)
#              B. L3 data reader  (from plotting_xco2_l3_leo)
#              C. VIIRS reader    (from plotting_xco2_l3_leo)
#              D. Spatial helpers (aggregation, upsampling, grid utils)
#              E. Plot helpers    (theme, axes, color scales)
#              F. Stats helpers   (VIF pruning, z-score, interpolation)
#   SECTION 3  Main loop — data extraction (L3 + VIIRS co-location)
#   SECTION 4  Statistical modeling (VIF + AICc dredge, pooled)
#   SECTION 5  XCO2 bias correction (full, non-scene-mean-preserving)
#   SECTION 6  Output plots and CSV export
#   SECTION 7  Final summary
#
# Data flow
# ---------
#   For every (TARGET_ID, SCENE_DATE) in `runs`:
#     1. Locate the highest-p L3 regrid NetCDF (45 m grid, MSAT L3 prefix).
#     2. Read L3: XCO2, dP, albedo, geometry, 2-D lat/lon.
#     3. Read matching VIIRS AERDB granule(s) from VIIRS_DATA_DIR.
#     4. Rasterize both datasets onto a shared 0.05° grid.
#     5. Merge, compute derived predictors (AOD_1286, cos angles).
#   After all scenes are extracted the pooled dataset goes through:
#     - VIF pruning  → AICc model selection (MuMIn::dredge)
#     - Full bias correction: xco2_corrected = xco2 − ŷ + reference_xco2
#     - Within-swath interpolation (akima) to fill VIIRS cloud gaps.
#   Outputs: CSV tables, before/after maps, bias maps, histograms.
# =============================================================================


# =============================================================================
# SECTION 0 — Libraries
# =============================================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(MASS)
  library(broom)
  library(MuMIn)
  library(car)
  library(akima)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(writexl)
  library(ggplot2)
  library(scales)
  library(tibble)
  library(tools)
  library(stringr)
  library(grid)
  library(gridExtra)
})


# =============================================================================
# SECTION 1 — User Settings
# =============================================================================

# ---- Scene list ----
# Set selected_targets = TRUE to process a hand-picked list.
# Set to FALSE to discover all scenes under ROOT_L3_DIR.
selected_targets <- TRUE
starting_id      <- NA   # only used when selected_targets = FALSE

if (selected_targets) {
  runs <- tibble::tribble(
    ~TARGET_ID, ~SCENE_DATE,
    56,  "20241223",
    349, "20250516",
    349, "20241230"
  )
} else {
  all_ids    <- list.files("/mnt/gcs/msat-prod-data-methanesat-level3/regrid")
  target_ids <- as.integer(sub("^t", "", all_ids))
  target_ids <- target_ids[which(target_ids == starting_id):length(target_ids)]
  runs_list  <- list()
  for (tid in target_ids) {
    scene_dates <- list.files(
      file.path("/mnt/gcs/msat-prod-data-methanesat-level3/regrid",
                paste0("t", tid))
    )
    for (sd in scene_dates)
      runs_list[[length(runs_list) + 1]] <-
        data.frame(TARGET_ID = tid, SCENE_DATE = sd, stringsAsFactors = FALSE)
  }
  runs <- as_tibble(do.call(rbind, runs_list))
}

# ---- Directory paths ----
ROOT_L3_DIR   <- "/mnt/gcs/msat-prod-data-methanesat-level3/regrid"
VIIRS_DATA_DIR <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/viirs_data2"
OUT_ROOT      <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/scripts/L3_results"
OUTPUT_DIR    <- file.path(OUT_ROOT,
                           paste0("bias_correction_", format(Sys.Date(), "%Y%m%d")))

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ---- Color palette ----
SSEC_PALETTE_FILE <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/scripts/ssec.256.txt"
ssec.256 <- if (file.exists(SSEC_PALETTE_FILE)) readLines(SSEC_PALETTE_FILE) else NULL

# ---- Binning / display ----
AGG_R        <- 10       # row block size for L3 aggregation
AGG_C        <- 10       # col block size
AGG_MIN_FRAC <- 0.1      # minimum fill fraction to keep a block
UPSAMPLE_ENABLE <- TRUE
UPSAMPLE_FACTOR <- 4
FONT_SCALE   <- 1.7

# ---- Modeling constants ----
VIF_THRESHOLD       <- 5
TOP_N_DREDGE        <- 10
TEMPORAL_WINDOW_HOURS <- 1.5
INTERP_MAX_GAP_DEG  <- 0.15   # max distance to fill via Delaunay interpolation

# ---- Quality thresholds ----
MIN_OVERLAP_PIXELS <- 100
AE_MIN             <- -0.5
AE_MAX             <-  3.0


# =============================================================================
# SECTION 2A — L3 File Discovery Helpers
# (verbatim from plotting_xco2_l3_leo)
# =============================================================================

find_l3_regrid_nc <- function(root_dir, target_id, scene_date,
                               grid_dir = "45m", pattern = "\\.nc$",
                               prefer_msat_prefix = TRUE,
                               keep_only_highest_p = TRUE) {
  base_dir   <- file.path(root_dir, paste0("t", target_id), scene_date)
  if (!dir.exists(base_dir)) return(character(0))
  candidates <- list.files(base_dir, recursive = TRUE,
                            full.names = TRUE, pattern = pattern)
  if (length(candidates) == 0) return(character(0))
  candidates <- candidates[grepl(paste0("/p\\d+/", grid_dir, "/"), candidates)]
  if (length(candidates) == 0) return(character(0))
  if (isTRUE(prefer_msat_prefix)) {
    pref <- candidates[grepl("^MSAT_L3_", basename(candidates))]
    if (length(pref) > 0) candidates <- pref
  }
  if (isTRUE(keep_only_highest_p)) {
    pp <- vapply(candidates, p_number, integer(1))
    if (any(is.finite(pp))) {
      pmax_val <- max(pp[is.finite(pp)])
      candidates <- candidates[pp == pmax_val]
    }
  }
  candidates
}

p_number <- function(path) {
  m <- regexpr("p\\d+", path)
  if (m[1] == -1) return(NA_integer_)
  as.integer(sub("^p", "", regmatches(path, m)))
}

pick_latest_mtime <- function(paths) {
  if (length(paths) <= 1) return(paths)
  mt <- file.info(paths)$mtime
  paths[which.max(mt)]
}


# =============================================================================
# SECTION 2B — L3 Data Reader
# (verbatim from plotting_xco2_l3_leo)
# =============================================================================

read_msat_l3_regrid <- function(nc_path) {
  stopifnot(file.exists(nc_path))
  nc <- ncdf4::nc_open(nc_path)
  on.exit(try(ncdf4::nc_close(nc), silent = TRUE), add = TRUE)

  tryGet_any <- function(names, default = NA_real_) {
    for (nm in names) {
      out <- try(ncdf4::ncvar_get(nc, nm), silent = TRUE)
      if (!inherits(out, "try-error")) return(out)
    }
    default
  }

  to_latlon <- function(x_lonlat) {
    if (is.null(dim(x_lonlat)) || length(dim(x_lonlat)) != 2) return(x_lonlat)
    aperm(x_lonlat, c(2, 1))
  }

  make_empty_like <- function(lat_len, lon_len)
    matrix(NA_real_, nrow = lat_len, ncol = lon_len)

  # 1-D coordinate vectors → 2-D matrices
  lat_vec <- tryGet_any(c("lat", "latitude"))
  lon_vec <- tryGet_any(c("lon", "longitude"))
  if (!is.numeric(lat_vec) || !is.numeric(lon_vec) ||
      length(lat_vec) < 2  || length(lon_vec) < 2)
    stop("Could not read 1D lat/lon from: ", nc_path)

  lat_mat <- matrix(rep(lat_vec, times = length(lon_vec)),
                    nrow = length(lat_vec), ncol = length(lon_vec),
                    byrow = FALSE)
  lon_mat <- matrix(rep(lon_vec, each = length(lat_vec)),
                    nrow = length(lat_vec), ncol = length(lon_vec),
                    byrow = FALSE)
  empty2d <- make_empty_like(length(lat_vec), length(lon_vec))

  # Geophysical fields
  dp      <- to_latlon(tryGet_any(
    c("o2dp_fit_diagnostics/bias_corrected_delta_pressure",
      "o2dp_fit_diagnostics/delta_pressure", "delta_pressure"),
    default = empty2d))

  alb1606 <- to_latlon(tryGet_any(
    c("co2proxy_fit_diagnostics/retrieved_albedo_1606nm",
      "apriori_data/albedo_ch4band", "albedo_ch4band"),
    default = empty2d))

  alb1286 <- to_latlon(tryGet_any(
    c("o2dp_fit_diagnostics/retrieved_albedo_1286nm"),
    default = empty2d))

  psfc <- to_latlon(tryGet_any(
    c("co2proxy_fit_diagnostics/retrieved_surface_pressure"),
    default = empty2d))

  orog <- to_latlon(tryGet_any(
    c("geolocation/terrain_height"),
    default = empty2d))

  sza  <- to_latlon(tryGet_any(c("geolocation/sza"), default = empty2d))
  vza  <- to_latlon(tryGet_any(c("geolocation/vza"), default = empty2d))
  saa  <- to_latlon(tryGet_any(c("geolocation/saa"), default = empty2d))
  vaa  <- to_latlon(tryGet_any(c("geolocation/vaa"), default = empty2d))

  # Quality flag (1 = good)
  qf <- to_latlon(tryGet_any(
    c("product_o2proxy/main_quality_flag"),
    default = empty2d))

  # XCO2 via O2-proxy
  co2_vcd_bc <- to_latlon(tryGet_any(
    c("co2proxy_fit_diagnostics/bias_corrected_co2_vcd",
      "bias_corrected_co2_vcd"), default = empty2d))
  o2_vcd_bc  <- to_latlon(tryGet_any(
    c("o2dp_fit_diagnostics/bias_corrected_o2_vcd",
      "bias_corrected_o2_vcd"), default = empty2d))
  ch4_vcd_bc <- to_latlon(tryGet_any(
    c("co2proxy_fit_diagnostics/bias_corrected_ch4_vcd",
      "bias_corrected_ch4_vcd"), default = empty2d))
  xch4_native <- to_latlon(tryGet_any(c("xch4"), default = empty2d))

  xco2_o2 <- empty2d
  ok_o2   <- is.finite(o2_vcd_bc) & (o2_vcd_bc != 0) &
             is.finite(co2_vcd_bc) & is.finite(ch4_vcd_bc)
  xco2_o2[ok_o2]  <- (co2_vcd_bc[ok_o2] / o2_vcd_bc[ok_o2]) * 0.2095
  xco2_o2_ppm     <- xco2_o2 * 1e6

  # Range masks
  xco2_o2_ppm[xco2_o2_ppm < 300 | xco2_o2_ppm > 600] <- NA
  if (!all(is.na(dp)))     dp[abs(dp) > 1e10]             <- NA
  if (!all(is.na(alb1286))) alb1286[alb1286 < 0 | alb1286 > 1] <- NA
  if (!all(is.na(alb1606))) alb1606[alb1606 < 0 | alb1606 > 1] <- NA
  if (!all(is.na(psfc)))   psfc[psfc < 50000 | psfc > 110000]   <- NA
  if (!all(is.na(orog)))   orog[abs(orog) > 10000]              <- NA
  if (!all(is.na(sza)))    sza[sza < 0 | sza > 90]              <- NA
  if (!all(is.na(vza)))    vza[vza < 0 | vza > 90]              <- NA
  if (!all(is.na(saa)))    saa[saa < 0 | saa > 360]             <- NA
  if (!all(is.na(vaa)))    vaa[vaa < 0 | vaa > 360]             <- NA

  # Apply quality flag
  if (!all(is.na(qf))) {
    bad <- !is.finite(qf) | (qf != 1)
    xco2_o2_ppm[bad] <- NA
    dp[bad]      <- NA
    alb1286[bad] <- NA
    alb1606[bad] <- NA
    psfc[bad]    <- NA
    orog[bad]    <- NA
    sza[bad]     <- NA
    vza[bad]     <- NA
    saa[bad]     <- NA
    vaa[bad]     <- NA
    message(sprintf("    Quality filter: %d / %d pixels pass QF=1",
                    sum(!bad), length(bad)))
  }

  list(
    lat_mat     = lat_mat,
    lon_mat     = lon_mat,
    lat_vec     = lat_vec,
    lon_vec     = lon_vec,
    xco2_o2_ppm = xco2_o2_ppm,
    dp          = dp,
    alb1286     = alb1286,
    alb1606     = alb1606,
    psfc        = psfc,
    orog        = orog,
    sza         = sza,
    vza         = vza,
    saa         = saa,
    vaa         = vaa,
    xch4_native = xch4_native   # used for valid-pixel mask in plotting
  )
}


# =============================================================================
# SECTION 2C — VIIRS AERDB Reader
# (verbatim from plotting_xco2_l3_leo, extended with AE reading from 3_4)
# =============================================================================

read_viirs_aerdb_full <- function(nc_path) {
  # Returns list with: aod, ae, lat, lon (all vectors, NAs masked)
  stopifnot(file.exists(nc_path))
  nc <- ncdf4::nc_open(nc_path)
  on.exit(try(ncdf4::nc_close(nc), silent = TRUE), add = TRUE)

  safe_get <- function(name) {
    out <- try(ncdf4::ncvar_get(nc, name), silent = TRUE)
    if (inherits(out, "try-error")) NULL else out
  }

  aod_raw <- safe_get("Aerosol_Optical_Thickness_550_Land_Best_Estimate")
  ae_raw  <- safe_get("Angstrom_Exponent_Land_Best_Estimate")
  lat_raw <- safe_get("Latitude")
  lon_raw <- safe_get("Longitude")

  if (is.null(aod_raw))
    stop("Could not read AOD variable from: ", nc_path)

  # ncdf4 returns (XDim, YDim); transpose to (YDim, XDim)
  d_aod <- dim(aod_raw)
  aod_mat <- if (length(d_aod) == 2) t(as.matrix(aod_raw)) else as.matrix(aod_raw)
  nr <- nrow(aod_mat); nc_aod <- ncol(aod_mat)

  fix_dims <- function(x, target_nrow, target_ncol) {
    if (is.null(x)) return(matrix(NA_real_, nrow = target_nrow, ncol = target_ncol))
    d <- dim(x)
    if (length(d) == 2) {
      if (d[1] == target_nrow && d[2] == target_ncol) return(as.matrix(x))
      if (d[1] == target_ncol && d[2] == target_nrow) return(t(as.matrix(x)))
    }
    matrix(NA_real_, nrow = target_nrow, ncol = target_ncol)
  }

  ae_mat  <- fix_dims(ae_raw,  nr, nc_aod)
  lat_mat <- fix_dims(lat_raw, nr, nc_aod)
  lon_mat <- fix_dims(lon_raw, nr, nc_aod)

  # Fill value handling for AOD
  fill_aod <- ncdf4::ncatt_get(nc, "Aerosol_Optical_Thickness_550_Land_Best_Estimate",
                                "_FillValue")
  if (fill_aod$hasatt)
    aod_mat[abs(aod_mat - fill_aod$value) < 1e-5] <- NA_real_
  aod_mat[aod_mat <= -999 | aod_mat < -0.1 | aod_mat > 5] <- NA_real_

  # Fill value handling for AE
  if (!is.null(ae_raw)) {
    fill_ae <- ncdf4::ncatt_get(nc, "Angstrom_Exponent_Land_Best_Estimate", "_FillValue")
    if (fill_ae$hasatt)
      ae_mat[abs(ae_mat - fill_ae$value) < 1e-5] <- NA_real_
    ae_mat[ae_mat <= -999 | ae_mat < AE_MIN | ae_mat > AE_MAX] <- NA_real_
  }

  list(
    aod  = as.vector(aod_mat),
    ae   = as.vector(ae_mat),
    lat  = as.vector(lat_mat),
    lon  = as.vector(lon_mat),
    fname = basename(nc_path)
  )
}

# Infer UTC acquisition time from AERDB filename (e.g. A2024357.2100)
extract_viirs_time <- function(viirs_file) {
  fn         <- basename(viirs_file)
  time_match <- regmatches(fn, regexpr("A(\\d{7})\\.(\\d{4})", fn))
  if (length(time_match) == 0) return(NULL)
  parts  <- strsplit(time_match, "\\.")[[1]]
  yyyydoy <- substr(parts[1], 2, 8)
  hhmm    <- parts[2]
  year    <- as.integer(substr(yyyydoy, 1, 4))
  doy     <- as.integer(substr(yyyydoy, 5, 7))
  hour    <- as.integer(substr(hhmm, 1, 2))
  minute  <- as.integer(substr(hhmm, 3, 4))
  date_part <- as.Date(doy - 1, origin = paste0(year, "-01-01"))
  as.POSIXct(paste(date_part, sprintf("%02d:%02d:00", hour, minute)), tz = "UTC")
}

# Find VIIRS granules matching a scene date and within the temporal window
find_viirs_files <- function(scene_date, scene_time_median) {
  all_files <- list.files(VIIRS_DATA_DIR, pattern = "\\.nc$",
                           full.names = TRUE, recursive = TRUE)
  if (length(all_files) == 0) {
    warning("No VIIRS files found in: ", VIIRS_DATA_DIR)
    return(character(0))
  }

  # Filter by date
  day_files <- Filter(function(f) {
    fn    <- basename(f)
    m     <- regmatches(fn, regexpr("A(\\d{7})", fn))
    if (length(m) == 0) return(FALSE)
    yyyydoy <- substr(m, 2, 8)
    yr  <- as.integer(substr(yyyydoy, 1, 4))
    doy <- as.integer(substr(yyyydoy, 5, 7))
    gd  <- as.Date(doy - 1, origin = paste0(yr, "-01-01"))
    !is.na(gd) && gd == scene_date
  }, all_files)

  if (length(day_files) == 0) return(character(0))

  # Filter by temporal window
  Filter(function(f) {
    vt <- extract_viirs_time(f)
    if (is.null(vt)) return(FALSE)
    abs(as.numeric(difftime(vt, scene_time_median, units = "hours"))) <=
      TEMPORAL_WINDOW_HOURS
  }, day_files)
}


# =============================================================================
# SECTION 2D — Spatial / Matrix Helpers
# (verbatim from plotting_xco2_l3_leo)
# =============================================================================

collapse_to_2d <- function(x, how = c("mean", "first")) {
  how <- match.arg(how)
  if (is.matrix(x)) return(x)
  if (!is.array(x) || is.null(dim(x))) return(x)
  d <- dim(x)
  if (length(d) == 2) return(as.matrix(x))
  if (how == "first") {
    idx <- c(list(TRUE), rep(list(1), length(d) - 1))
    return(as.matrix(do.call(`[`, c(list(x), idx, list(drop = TRUE)))))
  } else {
    return(as.matrix(apply(x, c(1, 2), mean, na.rm = TRUE)))
  }
}

aggregate_matrix <- function(m, row_block_size = AGG_R, col_block_size = AGG_C,
                              min_frac = AGG_MIN_FRAC) {
  if (!is.matrix(m) || all(is.na(m))) return(m)
  nr <- nrow(m); nc_m <- ncol(m)
  block_n <- row_block_size * col_block_size

  if (nr %% row_block_size == 0 && nc_m %% col_block_size == 0) {
    nr_out  <- nr  / row_block_size
    nc_out  <- nc_m / col_block_size
    a        <- array(m, dim = c(row_block_size, nr_out, col_block_size, nc_out))
    out_mean <- apply(a, c(2, 4), mean, na.rm = TRUE)
    a_valid  <- array(is.finite(m),
                      dim = c(row_block_size, nr_out, col_block_size, nc_out))
    out_nval <- apply(a_valid, c(2, 4), sum)
    out_mean[out_nval < min_frac * block_n] <- NA_real_
    return(out_mean)
  }

  r        <- terra::rast(m)
  agg      <- terra::aggregate(r, fact = c(row_block_size, col_block_size),
                                fun = "mean", na.rm = TRUE)
  result   <- terra::as.matrix(agg, wide = TRUE)
  r_valid  <- terra::rast(is.finite(m) * 1L)
  agg_nval <- terra::aggregate(r_valid,
                                fact = c(row_block_size, col_block_size),
                                fun = "sum", na.rm = TRUE)
  nval_mat <- terra::as.matrix(agg_nval, wide = TRUE)
  result[nval_mat < min_frac * block_n] <- NA_real_
  result
}

mat_to_df <- function(m, name = "val") {
  if (!is.matrix(m)) return(tibble())
  nr <- nrow(m); nc_m <- ncol(m)
  ij <- expand.grid(i = seq_len(nr), j = seq_len(nc_m))
  tibble(i = ij$i, j = ij$j, !!name := as.vector(m))
}

upsample_matrix <- function(mat, factor = 4, method = "bilinear") {
  if (!is.matrix(mat)) return(NULL)
  nr <- nrow(mat); nc_m <- ncol(mat)
  r   <- terra::rast(mat)
  terra::ext(r) <- terra::ext(1, nr, 1, nc_m)
  r_t <- terra::rast(nrows = nr * factor, ncols = nc_m * factor)
  terra::ext(r_t) <- terra::ext(r)
  r_up <- terra::resample(r, r_t, method = method)
  terra::as.matrix(r_up, wide = TRUE)
}

prep_plot_matrix <- function(m, kind = c("cont", "near"), factor = UPSAMPLE_FACTOR) {
  kind <- match.arg(kind)
  if (!is.matrix(m)) return(list(mat = m, factor = 1L))
  if (!isTRUE(UPSAMPLE_ENABLE) || !is.finite(factor) || factor <= 1)
    return(list(mat = m, factor = 1L))
  method <- if (kind == "near") "near" else "bilinear"
  m_up   <- upsample_matrix(m, factor = as.integer(factor), method = method)
  if (!is.matrix(m_up)) return(list(mat = m, factor = 1L))
  list(mat = m_up, factor = as.integer(factor))
}

# Compute scattering angle (degrees) from SZA, VZA, SAA, VAA
calc_primary_angle <- function(sza, vza, saa, vaa) {
  dtor <- pi / 180
  cos_phase <- cos(sza * dtor) * cos(vza * dtor) +
               sin(sza * dtor) * sin(vza * dtor) * cos((saa - vaa) * dtor)
  acos(pmax(-1, pmin(1, cos_phase))) * (180 / pi)
}


# =============================================================================
# SECTION 2E — Plot Helpers
# (verbatim from plotting_xco2_l3_leo)
# =============================================================================

pct_limits <- function(x, lo = 0.05, hi = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NULL)
  as.numeric(stats::quantile(x, probs = c(lo, hi), na.rm = TRUE))
}

add_fill_scale_ssec <- function(lims, colors = ssec.256) {
  if (is.null(colors)) colors <- viridisLite::viridis(256)
  if (is.null(lims) || any(!is.finite(lims)) || lims[1] >= lims[2])
    scale_fill_gradientn(colors = colors, na.value = "transparent")
  else
    scale_fill_gradientn(colors = colors, na.value = "transparent",
                         limits = lims, oob = scales::squish)
}

add_fill_scale <- function(lims) {
  if (is.null(lims) || any(!is.finite(lims)) || lims[1] >= lims[2])
    scale_fill_viridis_c(option = "C", na.value = "transparent")
  else
    scale_fill_viridis_c(option = "C", na.value = "transparent",
                         limits = lims, oob = scales::squish)
}

fmt_deg <- function(x, digits = 2) paste0(round(as.numeric(x), digits = digits), "\u00B0")

add_latlon_axes_plot <- function(lon_b, lat_b, plot_nr, plot_nc,
                                  n_ticks = 5, digits = 2,
                                  agg_fun = median, na_rm = TRUE) {
  if (!is.matrix(lon_b) || !is.matrix(lat_b) ||
      !is.finite(plot_nr) || !is.finite(plot_nc)) {
    return(list(scale_x_continuous(expand = c(0, 0)),
                scale_y_continuous(expand = c(0, 0))))
  }
  nr0 <- nrow(lon_b); nc0 <- ncol(lon_b)
  f_row <- if (nr0 > 0) plot_nr / nr0 else 1
  f_col <- if (nc0 > 0) plot_nc / nc0 else 1
  f     <- as.integer(max(1, round(min(f_row, f_col))))
  x0    <- unique(pmax(1, pmin(nc0, pretty(seq_len(nc0), n = n_ticks))))
  y0    <- unique(pmax(1, pmin(nr0, pretty(seq_len(nr0), n = n_ticks))))
  x_breaks <- (x0 - 1) * f + 1
  y_breaks <- (y0 - 1) * f + 1
  lon_by_col <- apply(lon_b, 2, function(v) {
    if (na_rm) v <- v[is.finite(v)]
    if (length(v) == 0) NA_real_ else agg_fun(v, na.rm = TRUE)
  })
  lat_by_row <- apply(lat_b, 1, function(v) {
    if (na_rm) v <- v[is.finite(v)]
    if (length(v) == 0) NA_real_ else agg_fun(v, na.rm = TRUE)
  })
  list(
    scale_x_continuous(expand = c(0, 0), breaks = x_breaks,
                       labels = fmt_deg(lon_by_col[x0], digits)),
    scale_y_continuous(expand = c(0, 0), breaks = y_breaks,
                       labels = fmt_deg(lat_by_row[y0], digits))
  )
}

theme_msat_map <- function(base_size = 12) {
  theme_minimal(base_size = base_size * FONT_SCALE) +
    theme(
      plot.title.position = "plot",
      plot.title    = element_text(face = "bold",  size = rel(1.1)),
      plot.subtitle = element_text(size = rel(0.95), margin = margin(b = 6)),
      plot.caption  = element_text(size = rel(0.8), color = "grey35", hjust = 0),
      axis.title    = element_text(size = rel(0.9)),
      axis.text     = element_text(size = rel(0.85), color = "grey20"),
      panel.grid.major = element_line(linewidth = 0.25, color = "grey92"),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(size = rel(0.9)),
      legend.text      = element_text(size = rel(0.85)),
      legend.key.height = unit(28 * FONT_SCALE, "pt"),
      legend.key.width  = unit(10 * FONT_SCALE, "pt")
    )
}

plot_caption <- function(TARGET_ID, SCENE_DATE, nc_file, extra = NULL) {
  cap <- paste0("Target t", TARGET_ID, " | Scene ", SCENE_DATE,
                " | File: ", basename(nc_file))
  if (!is.null(extra) && nzchar(extra)) cap <- paste0(cap, " | ", extra)
  cap
}

plot_matrix_map <- function(mat, nm, lon_b, lat_b,
                             fill_lab = "", use_ssec = FALSE,
                             caption_text = NULL,
                             upsample_for_plot = FALSE,
                             upsample_kind = c("cont", "near")) {
  if (!is.matrix(mat) || !is.matrix(lon_b) || !is.matrix(lat_b)) return(NULL)
  upsample_kind <- match.arg(upsample_kind)
  pp       <- if (isTRUE(upsample_for_plot))
    prep_plot_matrix(mat, kind = upsample_kind) else list(mat = mat, factor = 1L)
  mat_plot <- pp$mat
  d        <- mat_to_df(mat_plot, "val")
  fill_sc  <- if (use_ssec) add_fill_scale_ssec(pct_limits(d$val)) else
                             add_fill_scale(pct_limits(d$val))
  ggplot(d, aes(x = j, y = i, fill = val)) +
    geom_raster(interpolate = FALSE) +
    coord_equal(expand = FALSE) +
    fill_sc +
    add_latlon_axes_plot(lon_b, lat_b,
                         plot_nr = nrow(mat_plot), plot_nc = ncol(mat_plot),
                         n_ticks = 5, digits = 2) +
    labs(title   = nm,
         fill    = fill_lab,
         x       = "Approx. Longitude",
         y       = "Approx. Latitude",
         caption = if (!is.null(caption_text)) caption_text else "") +
    theme_msat_map()
}

# Side-by-side before/after map for corrected XCO2 (lon/lat geo-coordinates)
create_map_pair_geo <- function(df, var1, var2, title1, title2, main_title,
                                 tile_size = 0.05) {
  all_vals    <- c(df[[var1]], df[[var2]])
  scene_mean  <- mean(all_vals, na.rm = TRUE)
  scene_sd    <- sd(all_vals,   na.rm = TRUE)
  color_limits <- c(
    max(min(all_vals, na.rm = TRUE), scene_mean - 3 * scene_sd),
    min(max(all_vals, na.rm = TRUE), scene_mean + 3 * scene_sd)
  )

  make_p <- function(var, title) {
    ggplot(df, aes(x = x, y = y, fill = .data[[var]])) +
      geom_tile(width = tile_size, height = tile_size) +
      scale_fill_viridis_c(option = "magma", name = "ppm",
                           limits = color_limits, oob = scales::squish) +
      coord_fixed() +
      labs(x = "Longitude", y = "Latitude", title = title) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position  = "bottom",
        panel.grid       = element_blank(),
        plot.title       = element_text(hjust = 0.5, face = "bold"),
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(0.4, "cm")
      )
  }

  grob <- arrangeGrob(
    make_p(var1, title1), make_p(var2, title2), ncol = 2,
    top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
  )
  grob
}


# =============================================================================
# SECTION 2F — Statistics Helpers
# (verbatim from 3_4_script.R)
# =============================================================================

z_score_apply <- function(x, mean_val, sd_val) {
  if (sd_val < 1e-10) return(rep(NA_real_, length(x)))
  (x - mean_val) / sd_val
}

max_vif_num <- function(df, xvars) {
  if (length(xvars) < 2) return(1.0)
  form    <- as.formula(paste("xco2 ~", paste(xvars, collapse = "+")))
  dat_fit <- if (nrow(df) > 5000) df[sample(nrow(df), 5000), ] else df
  tryCatch({
    mod  <- lm(form, data = dat_fit)
    vifs <- car::vif(mod)
    if (is.null(dim(vifs))) max(vifs) else max(vifs[, 1])
  }, error = function(e) NA_real_)
}

vif_prune <- function(df, xvars, threshold = VIF_THRESHOLD) {
  keep <- xvars
  repeat {
    if (length(keep) < 2) break
    curr_max <- max_vif_num(df, keep)
    if (is.na(curr_max) || curr_max <= threshold) break
    form    <- as.formula(paste("xco2 ~", paste(keep, collapse = "+")))
    dat_fit <- if (nrow(df) > 5000) df[sample(nrow(df), 5000), ] else df
    mod     <- lm(form, data = dat_fit)
    v       <- car::vif(mod)
    val_to_check <- if (is.null(dim(v))) v else v[, 1]
    keep    <- setdiff(keep, names(which.max(val_to_check)))
  }
  list(keep = keep)
}

# Within-swath spatial interpolation (fills VIIRS cloud gaps)
interpolate_within_swath <- function(df, var, max_gap_deg = INTERP_MAX_GAP_DEG) {
  df_valid <- df %>% filter(!is.na(.data[[var]]))

  if (nrow(df_valid) < 3) {
    message(sprintf("    Skipping interpolation for %s: too few valid points (%d)",
                    var, nrow(df_valid)))
    df[[paste0(var, "_interpolated")]] <- FALSE
    return(df)
  }

  na_rows <- which(is.na(df[[var]]))
  if (length(na_rows) == 0) {
    df[[paste0(var, "_interpolated")]] <- FALSE
    return(df)
  }

  interp_result <- tryCatch(
    akima::interp(
      x      = df_valid$x,
      y      = df_valid$y,
      z      = df_valid[[var]],
      xo     = sort(unique(df$x)),
      yo     = sort(unique(df$y)),
      linear = TRUE,
      extrap = FALSE
    ),
    error = function(e) {
      message(sprintf("    Interpolation failed for %s: %s", var, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(interp_result)) {
    df[[paste0(var, "_interpolated")]] <- FALSE
    return(df)
  }

  interp_df <- expand.grid(x = interp_result$x, y = interp_result$y)
  interp_df$z_interp <- as.vector(interp_result$z)
  interp_df <- interp_df %>% filter(!is.na(z_interp))

  valid_coords <- as.matrix(df_valid[, c("x", "y")])
  df[[paste0(var, "_interpolated")]] <- FALSE
  filled <- 0

  for (i in na_rows) {
    px <- df$x[i]; py <- df$y[i]
    dists <- sqrt((valid_coords[, 1] - px)^2 + (valid_coords[, 2] - py)^2)
    if (min(dists) > max_gap_deg) next
    idx <- which(abs(interp_df$x - px) < 1e-9 & abs(interp_df$y - py) < 1e-9)
    if (length(idx) == 0 || is.na(interp_df$z_interp[idx[1]])) next
    df[[var]][i] <- interp_df$z_interp[idx[1]]
    df[[paste0(var, "_interpolated")]][i] <- TRUE
    filled <- filled + 1
  }

  message(sprintf("    Interpolated %d / %d missing %s pixels (max gap %.2f°)",
                  filled, length(na_rows), var, max_gap_deg))
  df
}


# =============================================================================
# SECTION 3 — Main Extraction Loop
# =============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("MethaneSAT XCO2 L3 — Aerosol Bias Correction Pipeline")
message(paste(rep("=", 70), collapse = ""))
message(sprintf("Output directory: %s\n", OUTPUT_DIR))

all_scenes_list      <- list()
extraction_metadata  <- list()
extraction_errors    <- list()
raster_attributes    <- list()

for (k in seq_len(nrow(runs))) {
  TARGET_ID  <- runs$TARGET_ID[k]
  SCENE_DATE <- runs$SCENE_DATE[k]
  scene_id   <- sprintf("t%d_%s", TARGET_ID, SCENE_DATE)

  message(sprintf("\n[%d/%d] Scene: %s", k, nrow(runs), scene_id))

  # ---- 3a. Locate L3 NetCDF ----
  nc_files <- find_l3_regrid_nc(ROOT_L3_DIR, TARGET_ID, SCENE_DATE, grid_dir = "45m")
  if (length(nc_files) == 0) {
    msg <- sprintf("No L3 regrid NC found for %s", scene_id)
    message("  SKIP: ", msg)
    extraction_errors[[scene_id]] <- msg
    next
  }
  nc_path <- pick_latest_mtime(nc_files)
  message("  L3 file: ", basename(nc_path))

  res <- tryCatch({

    # ---- 3b. Read L3 data ----
    message("  Reading L3 data...")
    l3 <- read_msat_l3_regrid(nc_path)

    # Scene centre time — parse from filename (YYYYMMDDTHHMMSSZ)
    time_match <- regmatches(basename(nc_path),
                              regexpr("\\d{8}T\\d{6}Z", basename(nc_path)))
    if (length(time_match) > 0) {
      scene_time_median <- as.POSIXct(time_match[1], format = "%Y%m%dT%H%M%SZ",
                                       tz = "UTC")
    } else {
      scene_time_median <- as.POSIXct(paste0(SCENE_DATE, " 12:00:00"), tz = "UTC")
      message("  WARNING: Could not parse time from filename, using noon UTC.")
    }
    scene_date_obj <- as.Date(SCENE_DATE, format = "%Y%m%d")

    message(sprintf("  Scene time: %s",
                    format(scene_time_median, "%Y-%m-%d %H:%M:%S UTC")))

    # ---- 3c. Find and load VIIRS granules ----
    message("  Finding VIIRS granules...")
    viirs_files <- find_viirs_files(scene_date_obj, scene_time_median)

    if (length(viirs_files) == 0)
      stop(sprintf("No VIIRS granules within ±%.1fh of overpass", TEMPORAL_WINDOW_HOURS))

    message(sprintf("  Found %d VIIRS granule(s)", length(viirs_files)))

    viirs_dfs <- lapply(viirs_files, function(vf) {
      vdata <- tryCatch(read_viirs_aerdb_full(vf), error = function(e) {
        message("  WARNING reading ", basename(vf), ": ", conditionMessage(e))
        NULL
      })
      if (is.null(vdata)) return(NULL)
      data.frame(lon   = vdata$lon,
                 lat   = vdata$lat,
                 AOD550 = vdata$aod,
                 AE    = vdata$ae) %>%
        filter(complete.cases(.))
    })
    viirs_combined <- bind_rows(Filter(Negate(is.null), viirs_dfs))

    if (nrow(viirs_combined) == 0)
      stop("No valid VIIRS pixels after loading all granules")

    message(sprintf("  VIIRS valid pixels: %d", nrow(viirs_combined)))

    # ---- 3d. Build flat MSAT point data frame ----
    message("  Building MSAT point data frame...")
    lat_v  <- as.vector(l3$lat_mat)
    lon_v  <- as.vector(l3$lon_mat)
    xco2_v <- as.vector(l3$xco2_o2_ppm)

    # Cosine-transformed geometry
    cos_sza_v <- cos_vza_v <- cos_pa_v <- rep(NA_real_, length(lat_v))
    if (!all(is.na(l3$sza))) cos_sza_v <- as.vector(cos(l3$sza * pi / 180))
    if (!all(is.na(l3$vza))) cos_vza_v <- as.vector(cos(l3$vza * pi / 180))
    if (!all(is.na(l3$sza)) && !all(is.na(l3$vza)) &&
        !all(is.na(l3$saa)) && !all(is.na(l3$vaa))) {
      pa_v    <- calc_primary_angle(as.vector(l3$sza), as.vector(l3$vza),
                                     as.vector(l3$saa), as.vector(l3$vaa))
      cos_pa_v <- cos(pa_v * pi / 180)
    }

    msat_df <- data.frame(
      lon               = lon_v,
      lat               = lat_v,
      xco2              = xco2_v,
      dp                = as.vector(l3$dp),
      alb1286           = as.vector(l3$alb1286),
      alb1606           = as.vector(l3$alb1606),
      psfc              = as.vector(l3$psfc),
      orog              = as.vector(l3$orog),
      cos_sza           = cos_sza_v,
      cos_vza           = cos_vza_v,
      cos_primary_angle = cos_pa_v
    ) %>% filter(complete.cases(lon, lat, xco2))

    if (nrow(msat_df) == 0) stop("MSAT has no valid pixels after filtering")
    message(sprintf("  MSAT valid pixels: %d", nrow(msat_df)))

    # ---- 3e. Rasterize to common 0.05° grid ----
    message("  Rasterizing to common 0.05° grid...")
    buf_lon <- (max(msat_df$lon) - min(msat_df$lon)) * 0.1
    buf_lat <- (max(msat_df$lat) - min(msat_df$lat)) * 0.1
    res_deg <- 0.05

    r_template <- terra::rast(
      nrows = ceiling(((max(msat_df$lat) + buf_lat) - (min(msat_df$lat) - buf_lat)) / res_deg),
      ncols = ceiling(((max(msat_df$lon) + buf_lon) - (min(msat_df$lon) - buf_lon)) / res_deg),
      xmin  = min(msat_df$lon) - buf_lon,
      xmax  = max(msat_df$lon) + buf_lon,
      ymin  = min(msat_df$lat) - buf_lat,
      ymax  = max(msat_df$lat) + buf_lat,
      crs   = "EPSG:4326"
    )

    viirs_vect <- terra::vect(viirs_combined, geom = c("lon", "lat"), crs = "EPSG:4326")
    v_raster   <- terra::rasterize(viirs_vect, r_template,
                                    field = c("AOD550", "AE"), fun = "mean", na.rm = TRUE)

    msat_vars  <- c("xco2", "dp", "alb1286", "alb1606", "psfc", "orog",
                    "cos_sza", "cos_vza", "cos_primary_angle")
    msat_pts   <- terra::vect(msat_df, geom = c("lon", "lat"), crs = "EPSG:4326")
    m_raster   <- terra::rasterize(msat_pts, r_template,
                                    field = msat_vars, fun = "mean", na.rm = TRUE)

    stack_all <- c(v_raster, m_raster)
    expected  <- c("AOD550", "AE", msat_vars)
    if (terra::nlyr(stack_all) == length(expected))
      names(stack_all) <- expected

    # ---- 3f. Convert to data frame, filter, derive AOD_1286 ----
    df_scene <- as.data.frame(stack_all, xy = TRUE, na.rm = FALSE) %>%
      filter(!is.na(AOD550) & !is.na(AE) & !is.na(xco2)) %>%
      filter(AE >= AE_MIN & AE <= AE_MAX)

    if (nrow(df_scene) < MIN_OVERLAP_PIXELS)
      stop(sprintf("Insufficient co-located pixels: %d < %d required",
                   nrow(df_scene), MIN_OVERLAP_PIXELS))

    # Derived predictor: AOD at 1286 nm via Angstrom exponent
    df_scene$AOD550[df_scene$AOD550 <= 0] <- NA
    df_scene$AOD_1286 <- df_scene$AOD550 * (1286 / 550)^(-df_scene$AE)
    df_scene$AOD_1286[!is.finite(df_scene$AOD_1286)] <- NA
    df_scene$scene_id <- scene_id

    # Attach native-resolution data and raster as attributes for plotting
    attr(df_scene, "native_msat_data") <- msat_df
    attr(df_scene, "native_viirs_data") <- viirs_combined
    attr(df_scene, "viirs_raster_aod")  <- stack_all[[1]]
    attr(df_scene, "n_viirs_granules")  <- length(viirs_files)
    attr(df_scene, "grid_dims")         <- c(terra::nrow(r_template), terra::ncol(r_template))
    attr(df_scene, "time_median")       <- scene_time_median
    attr(df_scene, "l3_info")           <- list(
      lat_mat = l3$lat_mat, lon_mat = l3$lon_mat,
      xco2 = l3$xco2_o2_ppm, dp = l3$dp,
      alb1606 = l3$alb1606, xch4_native = l3$xch4_native
    )

    message(sprintf("  SUCCESS: %d co-located pixels on %dx%d grid",
                    nrow(df_scene), attr(df_scene, "grid_dims")[1],
                    attr(df_scene, "grid_dims")[2]))
    df_scene

  }, error = function(e) {
    message(sprintf("  FAILED: %s", conditionMessage(e)))
    extraction_errors[[scene_id]] <<- conditionMessage(e)
    NULL
  })

  if (!is.null(res)) {
    all_scenes_list[[length(all_scenes_list) + 1]] <- res
    raster_attributes[[scene_id]] <- list(
      viirs_raster_aod = attr(res, "viirs_raster_aod"),
      n_viirs_granules = attr(res, "n_viirs_granules"),
      grid_dims        = attr(res, "grid_dims"),
      l3_info          = attr(res, "l3_info")
    )
    extraction_metadata[[scene_id]] <- list(
      n_pixels         = nrow(res),
      n_viirs_granules = attr(res, "n_viirs_granules"),
      grid_dims        = attr(res, "grid_dims"),
      time_median      = as.character(attr(res, "time_median"))
    )
  }
}

# ---- Save extraction errors ----
if (length(extraction_errors) > 0) {
  write_csv(
    data.frame(scene_id = names(extraction_errors),
               error    = unlist(extraction_errors),
               stringsAsFactors = FALSE),
    file.path(OUTPUT_DIR, "extraction_errors.csv")
  )
  message(sprintf("\nSaved extraction errors (%d failures).", length(extraction_errors)))
}

if (length(all_scenes_list) == 0)
  stop("No scenes extracted successfully. Check paths and VIIRS files.")

all_scenes_df <- bind_rows(all_scenes_list)

write_csv(all_scenes_df, file.path(OUTPUT_DIR, "unified_database.csv"))
message(sprintf("\nUnified database: %d rows from %d scene(s).",
                nrow(all_scenes_df), length(all_scenes_list)))

meta_df <- bind_rows(lapply(names(extraction_metadata), function(sid) {
  m <- extraction_metadata[[sid]]
  data.frame(scene_id = sid, n_pixels = m$n_pixels,
             n_viirs_granules = m$n_viirs_granules,
             grid_rows = m$grid_dims[1], grid_cols = m$grid_dims[2],
             time_median = m$time_median, stringsAsFactors = FALSE)
}))
write_csv(meta_df, file.path(OUTPUT_DIR, "extraction_metadata.csv"))


# =============================================================================
# SECTION 4 — Statistical Modeling (pooled VIF + AICc dredge)
# =============================================================================

message("\n", paste(rep("-", 60), collapse = ""))
message("SECTION 4: Statistical Modeling")
message(paste(rep("-", 60), collapse = ""))

essential_predictors <- c("AOD_1286", "AE")
predictors_all <- c("dp", "alb1286", "alb1606", "psfc", "orog",
                    "AOD_1286", "AE", "cos_primary_angle")

scaling_params <- list()   # populated below; used in Section 5
best_pool_mod  <- NULL
df_pool_z      <- NULL     # z-scored training frame; used in Section 5

if (nrow(all_scenes_df) < 500) {
  message("  Skipping: fewer than 500 co-located pixels in pooled dataset.")

} else {

  df_pool_clean <- all_scenes_df %>%
    dplyr::select(any_of(c("xco2", "scene_id", predictors_all))) %>%
    filter(complete.cases(across(any_of(c("xco2", "scene_id", essential_predictors)))))

  message(sprintf("  Rows after requiring essential columns: %d", nrow(df_pool_clean)))

  preds_pool <- intersect(names(df_pool_clean), predictors_all)
  preds_pool <- preds_pool[sapply(preds_pool, function(p) {
    sum(!is.na(df_pool_clean[[p]])) >= nrow(df_pool_clean) * 0.5
  })]
  message(sprintf("  Predictors with >=50%% valid data: %s", paste(preds_pool, collapse = ", ")))

  n_scenes <- length(unique(df_pool_clean$scene_id))
  message(sprintf("  Unique scenes: %d", n_scenes))

  if (n_scenes < 2 || length(preds_pool) < 2) {
    message("  Skipping: need >=2 scenes and >=2 valid predictors.")

  } else {

    df_pool_clean$scene_id <- as.factor(df_pool_clean$scene_id)

    df_pool_z <- df_pool_clean %>%
      dplyr::select(any_of(c("xco2", "scene_id", preds_pool))) %>%
      drop_na()

    message(sprintf("  Rows for modeling: %d", nrow(df_pool_z)))

    if (nrow(df_pool_z) < 500) {
      message("  Skipping: fewer than 500 complete rows.")

    } else {

      # Save pooled scaling parameters before z-scoring (reused in correction)
      for (v in preds_pool) {
        scaling_params[[v]] <- list(
          mean = mean(df_pool_z[[v]], na.rm = TRUE),
          sd   = sd(df_pool_z[[v]],   na.rm = TRUE)
        )
      }
      message(sprintf("  Saved scaling parameters for %d predictors.", length(scaling_params)))

      # Z-score all predictors using pooled statistics
      for (v in preds_pool) {
        df_pool_z[[v]] <- z_score_apply(df_pool_z[[v]],
                                         scaling_params[[v]]$mean,
                                         scaling_params[[v]]$sd)
      }

      # VIF pruning
      vp_pool    <- vif_prune(df_pool_z, preds_pool)
      kept_pool  <- vp_pool$keep
      message(sprintf("  Predictors after VIF pruning: %s", paste(kept_pool, collapse = ", ")))

      if (length(kept_pool) == 0) {
        message("  Skipping: all predictors removed by VIF.")

      } else {

        # Drop scenes with fewer than 10 pixels
        scene_counts  <- df_pool_z %>% group_by(scene_id) %>%
          summarise(n = n(), .groups = "drop") %>% arrange(n)
        message("\n  Scene pixel counts:"); print(scene_counts)

        df_pool_z <- df_pool_z %>%
          filter(scene_id %in% scene_counts$scene_id[scene_counts$n >= 10]) %>%
          mutate(scene_id = droplevels(scene_id))

        n_scenes_final <- length(unique(df_pool_z$scene_id))
        message(sprintf("  Final scenes with >=10 pixels: %d", n_scenes_final))

        if (n_scenes_final < 2) {
          message("  Skipping: need >=2 scenes with >=10 pixels.")

        } else {

          message("  Running AICc model selection...")
          form_global <- as.formula(paste("xco2 ~ scene_id +",
                                          paste(kept_pool, collapse = "+")))
          mod_global <- lm(form_global, data = df_pool_z, na.action = "na.fail")

          dredge_res <- MuMIn::dredge(mod_global, rank = "AICc",
                                       fixed = "scene_id", trace = FALSE,
                                       m.lim = c(1, length(kept_pool)))

          write_csv(as.data.frame(dredge_res),
                    file.path(OUTPUT_DIR, "pooled_dredge_results.csv"))

          best_pool_mod <- MuMIn::get.models(dredge_res, subset = 1)[[1]]
          r2_adj <- summary(best_pool_mod)$adj.r.squared

          final_preds <- names(coef(best_pool_mod))
          final_preds <- final_preds[!grepl("Intercept|scene_id", final_preds)]
          message(sprintf("  Best model adj. R² = %.3f", r2_adj))
          message(sprintf("  Final predictors: %s", paste(final_preds, collapse = ", ")))

          write_csv(broom::tidy(best_pool_mod),
                    file.path(OUTPUT_DIR, "best_model_coefficients.csv"))

          png(file.path(OUTPUT_DIR, "model_diagnostics.png"),
              width = 1200, height = 900, res = 120)
          par(mfrow = c(2, 2))
          plot(best_pool_mod, main = "Model Diagnostics")
          dev.off()
          message("  Model diagnostics saved.")
        }
      }
    }
  }
}


# =============================================================================
# SECTION 5 — XCO2 Bias Correction (Full Correction)
# =============================================================================

message("\n", paste(rep("-", 60), collapse = ""))
message("SECTION 5: XCO2 Bias Correction")
message(paste(rep("-", 60), collapse = ""))

all_scenes_corrected <- list()
correction_stats     <- list()

if (is.null(best_pool_mod)) {
  message("  Skipping: no valid model available.")
  all_scenes_corrected <- setNames(
    lapply(unique(all_scenes_df$scene_id), function(sid) {
      df <- all_scenes_df %>% filter(scene_id == sid)
      df$xco2_corrected <- df$xco2
      df
    }),
    unique(all_scenes_df$scene_id)
  )

} else {

  coef_vec        <- coef(best_pool_mod)
  predictor_coefs <- coef_vec[!grepl("^scene_id|^\\(Intercept\\)", names(coef_vec))]

  # Reference level: mean of z-scored XCO2 in training data
  reference_xco2  <- mean(df_pool_z$xco2, na.rm = TRUE)
  training_scenes <- unique(df_pool_z$scene_id)

  message(sprintf("  Predictor coefficients used: %s",
                  paste(names(predictor_coefs), collapse = ", ")))
  message(sprintf("  Reference XCO2 (training z-score mean): %.4f", reference_xco2))

  for (sid in unique(all_scenes_df$scene_id)) {

    df_scene <- all_scenes_df %>% filter(scene_id == sid)
    if (nrow(df_scene) < MIN_OVERLAP_PIXELS) {
      message(sprintf("  [%s] Skipping: too few pixels (%d).", sid, nrow(df_scene)))
      next
    }

    scene_coef_name <- paste0("scene_id", sid)
    if (!scene_coef_name %in% names(coef_vec)) {
      message(sprintf("  [%s] WARNING: scene not in training data — skipping.", sid))
      next
    }

    # Compute linear predictor for each pixel
    y_hat <- rep(coef_vec["(Intercept)"], nrow(df_scene))
    for (v in names(predictor_coefs)) {
      if (v %in% names(df_scene) && v %in% names(scaling_params)) {
        z_vals <- z_score_apply(df_scene[[v]],
                                 scaling_params[[v]]$mean,
                                 scaling_params[[v]]$sd)
        if (!any(is.na(z_vals)))
          y_hat <- y_hat + predictor_coefs[[v]] * z_vals
        else
          message(sprintf("  [%s] WARNING: NAs in z-scored %s — skipped.", sid, v))
      }
    }
    y_hat <- y_hat + coef_vec[scene_coef_name]

    if (any(is.infinite(y_hat)) || all(is.na(y_hat))) {
      warning(sprintf("  [%s] ERROR: degenerate predictions — skipping.", sid))
      next
    }

    # Full bias correction (not scene-mean-preserving)
    bias                   <- y_hat - reference_xco2
    df_scene$xco2_corrected <- df_scene$xco2 - bias

    correction_stats[[sid]] <- data.frame(
      scene_id         = sid,
      mean_bias        = mean(bias, na.rm = TRUE),
      sd_bias          = sd(bias,   na.rm = TRUE),
      min_bias         = min(bias,   na.rm = TRUE),
      max_bias         = max(bias,   na.rm = TRUE),
      xco2_mean_before = mean(df_scene$xco2,           na.rm = TRUE),
      xco2_mean_after  = mean(df_scene$xco2_corrected, na.rm = TRUE),
      xco2_sd_before   = sd(df_scene$xco2,             na.rm = TRUE),
      xco2_sd_after    = sd(df_scene$xco2_corrected,   na.rm = TRUE),
      variance_reduction = 1 - sd(df_scene$xco2_corrected, na.rm = TRUE) /
        sd(df_scene$xco2, na.rm = TRUE),
      mean_shift = mean(df_scene$xco2_corrected, na.rm = TRUE) -
        mean(df_scene$xco2, na.rm = TRUE),
      mean_shift_pct = 100 * (mean(df_scene$xco2_corrected, na.rm = TRUE) -
                                mean(df_scene$xco2,          na.rm = TRUE)) /
        mean(df_scene$xco2, na.rm = TRUE)
    )

    if (sid %in% names(raster_attributes)) {
      attr(df_scene, "viirs_raster_aod") <- raster_attributes[[sid]]$viirs_raster_aod
      attr(df_scene, "n_viirs_granules") <- raster_attributes[[sid]]$n_viirs_granules
      attr(df_scene, "grid_dims")        <- raster_attributes[[sid]]$grid_dims
      attr(df_scene, "l3_info")          <- raster_attributes[[sid]]$l3_info
    }
    all_scenes_corrected[[sid]] <- df_scene

    message(sprintf(
      "  [%s] mean bias = %.2f ppm | shift = %.2f ppm (%.2f%%) | var reduction = %.1f%%",
      sid,
      correction_stats[[sid]]$mean_bias,
      correction_stats[[sid]]$mean_shift,
      correction_stats[[sid]]$mean_shift_pct,
      correction_stats[[sid]]$variance_reduction * 100
    ))
  }

  if (length(correction_stats) > 0) {
    correction_summary <- bind_rows(correction_stats)
    write_csv(correction_summary, file.path(OUTPUT_DIR, "correction_statistics.csv"))
    message(sprintf("\nCorrection statistics saved (%d scenes).", nrow(correction_summary)))

    message("\n--- Bias Correction Summary ---")
    message(sprintf("  Avg bias removed:    %.2f ppm",
                    mean(correction_summary$mean_bias, na.rm = TRUE)))
    message(sprintf("  Avg mean shift:      %.2f ppm (%.2f%%)",
                    mean(correction_summary$mean_shift, na.rm = TRUE),
                    mean(correction_summary$mean_shift_pct, na.rm = TRUE)))
    message(sprintf("  Bias range:          [%.2f, %.2f] ppm",
                    min(correction_summary$mean_bias, na.rm = TRUE),
                    max(correction_summary$mean_bias, na.rm = TRUE)))
    message(sprintf("  Avg variance reduction: %.1f%%",
                    mean(correction_summary$variance_reduction, na.rm = TRUE) * 100))
  }

  # Rebuild all_scenes_df with corrected values
  all_scenes_df <- bind_rows(all_scenes_corrected)
}


# =============================================================================
# SECTION 6 — Output Plots
# =============================================================================
#
# Primary output per scene
# ------------------------
#   *_four_panel.png  — the main figure: four spatially registered panels on a
#                       shared lon/lat coordinate system, all at 0.05° resolution:
#                         Panel 1: VIIRS AOD at 550 nm
#                         Panel 2: MethaneSAT albedo at 1606 nm
#                         Panel 3: XCO2 (O2-proxy, uncorrected)
#                         Panel 4: XCO2 (O2-proxy, aerosol-corrected)
#
# Secondary diagnostics (retained from original pipeline)
# --------------------------------------------------------
#   *_bias_map.png    — spatial map of removed bias (xco2 − xco2_corrected)
#   *_histogram.png   — XCO2 distribution before/after (observed pixels only)
#   *_scatter.png     — pixel-level correction relationship
#
# Coordinate strategy
# -------------------
# VIIRS AOD and the two XCO2 fields are already on the 0.05° raster grid
# (columns "x", "y" from as.data.frame(stack_all, xy=TRUE)).
# Albedo (1606 nm) lives on the native L3 lat/lon matrix. We convert it to a
# tidy lon/lat data frame at its native resolution so all four panels share the
# same geographic axes via geom_tile(), making direct spatial comparison easy.
# =============================================================================

message("\n", paste(rep("-", 60), collapse = ""))
message("SECTION 6: Plots")
message(paste(rep("-", 60), collapse = ""))

FIGS_DIR <- file.path(OUTPUT_DIR, "figures")
if (!dir.exists(FIGS_DIR)) dir.create(FIGS_DIR, recursive = TRUE)

# ---- Shared ggplot2 theme for all four panels --------------------------------
theme_panel <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = base_size),
      plot.subtitle    = element_text(hjust = 0.5, size = base_size - 1,
                                      color = "grey30"),
      axis.title       = element_text(size = base_size - 1),
      axis.text        = element_text(size = base_size - 2),
      legend.position  = "bottom",
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 2),
      legend.key.width = unit(1.8, "cm"),
      legend.key.height = unit(0.35, "cm"),
      panel.grid       = element_blank()
    )
}

# ---- Helper: build one geo-tile panel ----------------------------------------
# df       : data frame with columns lon, lat, val
# title    : panel title string
# fill_lab : legend label
# colors   : colour ramp vector
# lims     : c(lo, hi) for the fill scale; NULL = auto from 2nd/98th percentile
make_geo_panel <- function(df, title, fill_lab, colors,
                            lims = NULL, tile_deg = 0.05) {
  df <- df %>% filter(is.finite(val))
  if (nrow(df) == 0) return(NULL)

  if (is.null(lims)) {
    lims <- as.numeric(quantile(df$val, probs = c(0.02, 0.98), na.rm = TRUE))
  }
  if (lims[1] >= lims[2]) lims <- range(df$val, na.rm = TRUE)

  n_valid <- nrow(df)
  mean_v  <- mean(df$val, na.rm = TRUE)

  ggplot(df, aes(x = lon, y = lat, fill = val)) +
    geom_tile(width = tile_deg, height = tile_deg) +
    scale_fill_gradientn(
      colors   = colors,
      limits   = lims,
      oob      = scales::squish,
      na.value = "grey85",
      name     = fill_lab
    ) +
    coord_fixed(expand = FALSE) +
    labs(
      title    = title,
      subtitle = sprintf("n = %s  |  mean = %s %s",
                         format(n_valid, big.mark = ","),
                         round(mean_v, 3), fill_lab),
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_panel()
}

# ---- Null-coalescing helper (base R lacks %||%) ------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Colour ramps ------------------------------------------------------------
# AOD:       green → yellow → red  (low to high aerosol load)
# Albedo:    grey → white → bright (dark to bright surface)
# XCO2:      uses ssec.256 if available, otherwise viridis magma
cols_aod    <- c("#1a9850", "#ffffbf", "#d73027")
cols_albedo <- c("#2c2c2c", "#888888", "#f5f5f5", "#ffffff")
cols_xco2   <- if (!is.null(ssec.256)) ssec.256 else
               scales::viridis_pal(option = "magma")(256)


# ==============================================================================
# Main per-scene plotting loop
# ==============================================================================

for (sid in unique(all_scenes_df$scene_id)) {

  message(sprintf("\n  [%s] Generating plots...", sid))

  # Pull the co-located data frame for this scene
  df_scene_raw <- all_scenes_df %>% filter(scene_id == sid)
  l3_info      <- raster_attributes[[sid]]$l3_info

  has_corrected <- "xco2_corrected" %in% names(df_scene_raw)
  has_coords    <- all(c("x", "y") %in% names(df_scene_raw))

  if (!has_coords) {
    message(sprintf("  [%s] Skipping: missing x/y coordinate columns.", sid))
    next
  }

  # ============================================================
  # 6a. Within-swath interpolation (fills VIIRS cloud gaps)
  # ============================================================
  df_scene <- df_scene_raw %>%
    filter(is.finite(xco2), xco2 > 300, xco2 < 600)

  if (has_corrected)
    df_scene <- df_scene %>%
      filter(is.finite(xco2_corrected),
             xco2_corrected > 300, xco2_corrected < 600)

  if (nrow(df_scene) < MIN_OVERLAP_PIXELS) {
    message(sprintf("  [%s] Skipping: too few valid pixels after filtering.", sid))
    next
  }

  message(sprintf("  [%s] Applying within-swath interpolation...", sid))
  df_scene <- interpolate_within_swath(df_scene, "xco2")
  if (has_corrected)
    df_scene <- interpolate_within_swath(df_scene, "xco2_corrected")

  # Rename raster grid coordinates to lon/lat for geo panels
  df_scene <- df_scene %>% rename(lon = x, lat = y)

  # ============================================================
  # 6b. Albedo from native L3 matrix → tidy lon/lat data frame
  # ============================================================
  df_alb <- NULL
  if (!is.null(l3_info)) {
    lat_mat  <- l3_info$lat_mat
    lon_mat  <- l3_info$lon_mat
    alb_mat  <- l3_info$alb1606
    xch4_mat <- l3_info$xch4_native

    # Apply the same valid-pixel mask used in the L3 reader
    mask_mat <- ifelse(is.finite(xch4_mat), 1L, 0L)
    alb_mat[mask_mat == 0L] <- NA_real_
    alb_mat[alb_mat < 0 | alb_mat > 1] <- NA_real_

    # Convert to tidy frame — use the native L3 resolution (not binned)
    # so the albedo panel has the finest spatial detail available.
    df_alb <- data.frame(
      lon = as.vector(lon_mat),
      lat = as.vector(lat_mat),
      val = as.vector(alb_mat)
    ) %>% filter(is.finite(val))

    # Tile size matches the native L3 grid spacing (degrees per cell)
    # Inferred from the coordinate vectors; fall back to 0.05° if unclear.
    lat_vec  <- l3_info$lat_vec %||% sort(unique(df_alb$lat))
    alb_tile <- if (length(lat_vec) > 1)
      abs(diff(lat_vec[1:2])) else 0.05
  }

  # ============================================================
  # 6c. Build the four panels
  # ============================================================

  # --- Panel 1: VIIRS AOD 550 nm ---
  df_aod <- df_scene %>%
    filter(is.finite(AOD550), AOD550 >= 0, AOD550 <= 3) %>%
    select(lon, lat, val = AOD550)

  p1 <- make_geo_panel(
    df       = df_aod,
    title    = "VIIRS AOD at 550 nm",
    fill_lab = "AOD",
    colors   = cols_aod,
    tile_deg = 0.05
  )

  # --- Panel 2: MethaneSAT Albedo at 1606 nm ---
  p2 <- if (!is.null(df_alb) && nrow(df_alb) > 0) {
    make_geo_panel(
      df       = df_alb,
      title    = "MSAT Albedo at 1606 nm",
      fill_lab = "Albedo",
      colors   = cols_albedo,
      tile_deg = alb_tile
    )
  } else {
    # Fallback: use the albedo column from the raster grid if available
    if ("alb1606" %in% names(df_scene)) {
      df_alb_fb <- df_scene %>%
        filter(is.finite(alb1606), alb1606 >= 0, alb1606 <= 1) %>%
        select(lon, lat, val = alb1606)
      make_geo_panel(df_alb_fb, "MSAT Albedo at 1606 nm",
                     "Albedo", cols_albedo, tile_deg = 0.05)
    } else NULL
  }

  # --- Panels 3 & 4: XCO2 — use a shared colour scale so the two maps are
  #     directly comparable. Compute limits from both fields together.
  xco2_vals <- c(
    df_scene$xco2[is.finite(df_scene$xco2)],
    if (has_corrected) df_scene$xco2_corrected[is.finite(df_scene$xco2_corrected)]
    else numeric(0)
  )
  xco2_mean <- mean(xco2_vals, na.rm = TRUE)
  xco2_sd   <- sd(xco2_vals,   na.rm = TRUE)
  xco2_lims <- c(
    max(min(xco2_vals, na.rm = TRUE), xco2_mean - 2.5 * xco2_sd),
    min(max(xco2_vals, na.rm = TRUE), xco2_mean + 2.5 * xco2_sd)
  )

  # --- Panel 3: XCO2 uncorrected ---
  df_xco2_raw <- df_scene %>%
    filter(is.finite(xco2)) %>%
    select(lon, lat, val = xco2)

  p3 <- make_geo_panel(
    df       = df_xco2_raw,
    title    = "XCO\u2082 (O\u2082-proxy)",
    fill_lab = "ppm",
    colors   = cols_xco2,
    lims     = xco2_lims,
    tile_deg = 0.05
  )

  # --- Panel 4: XCO2 corrected ---
  p4 <- if (has_corrected) {
    df_xco2_cor <- df_scene %>%
      filter(is.finite(xco2_corrected)) %>%
      select(lon, lat, val = xco2_corrected)
    make_geo_panel(
      df       = df_xco2_cor,
      title    = "XCO\u2082 (O\u2082-proxy, aerosol-corrected)",
      fill_lab = "ppm",
      colors   = cols_xco2,
      lims     = xco2_lims,
      tile_deg = 0.05
    )
  } else {
    # No correction available — repeat uncorrected with a note
    p3 + labs(title = "XCO\u2082 — no correction (model unavailable)")
  }

  # ============================================================
  # 6d. Assemble and save the four-panel figure
  # ============================================================
  panels       <- list(p1, p2, p3, p4)
  valid_panels <- Filter(Negate(is.null), panels)
  n_panels     <- length(valid_panels)

  if (n_panels == 0) {
    message(sprintf("  [%s] No valid panels to plot — skipping.", sid))
    next
  }

  main_title <- textGrob(
    sprintf("%s  —  VIIRS AOD | MSAT Albedo | XCO\u2082 | XCO\u2082 Corrected", sid),
    gp = gpar(fontsize = 13, fontface = "bold")
  )

  four_panel_path <- file.path(FIGS_DIR, sprintf("%s_four_panel.png", sid))
  grDevices::png(four_panel_path,
                  width  = 400 * n_panels,   # 400 px per panel
                  height = 520,
                  res    = 100)
  gridExtra::grid.arrange(
    grobs = valid_panels,
    nrow  = 1,
    ncol  = n_panels,
    top   = main_title
  )
  grDevices::dev.off()
  message(sprintf("  [%s] Saved four-panel figure: %s", sid, basename(four_panel_path)))

  # ============================================================
  # 6e. Secondary diagnostics
  # ============================================================

  if (!has_corrected) next

  # --- Spatial bias map ---
  df_scene$bias_ppm <- df_scene$xco2 - df_scene$xco2_corrected
  bias_mean  <- mean(df_scene$bias_ppm, na.rm = TRUE)
  bias_sd    <- sd(df_scene$bias_ppm,   na.rm = TRUE)
  bias_range <- range(df_scene$bias_ppm, na.rm = TRUE)
  df_filt    <- df_scene %>%
    filter(abs(bias_ppm - bias_mean) < 5 * bias_sd)

  if (nrow(df_filt) > 50) {
    bias_lim <- max(abs(df_filt$bias_ppm - bias_mean), na.rm = TRUE)
    p_bias   <- ggplot(df_filt, aes(x = lon, y = lat, fill = bias_ppm)) +
      geom_tile(width = 0.05, height = 0.05) +
      coord_fixed() +
      scale_fill_gradient2(
        low      = "#2166ac", mid = "white", high = "#d6604d",
        midpoint = bias_mean,
        limits   = c(bias_mean - bias_lim, bias_mean + bias_lim),
        name     = "Bias (ppm)"
      ) +
      labs(
        title    = sprintf("%s \u2014 Removed Bias", sid),
        subtitle = sprintf("Mean: %.2f ppm | SD: %.2f ppm | Range: [%.2f, %.2f] ppm",
                           bias_mean, bias_sd, bias_range[1], bias_range[2]),
        x = "Longitude", y = "Latitude"
      ) +
      theme_panel()
    ggsave(file.path(FIGS_DIR, sprintf("%s_bias_map.png", sid)),
           p_bias, width = 10, height = 8, dpi = 300, bg = "white")
    message(sprintf("  [%s] Saved bias map.", sid))
  }

  # --- Histogram: observed pixels only ---
  obs_raw <- df_scene$xco2[!df_scene$xco2_interpolated]
  obs_cor <- df_scene$xco2_corrected[!df_scene$xco2_corrected_interpolated]

  if (length(obs_raw) > 0 && length(obs_cor) > 0) {
    df_hist <- data.frame(
      value = c(obs_raw, obs_cor),
      type  = c(rep("Uncorrected", length(obs_raw)),
                rep("Corrected",   length(obs_cor)))
    )
    p_hist <- ggplot(df_hist, aes(x = value, fill = type)) +
      geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
      scale_fill_manual(
        values = c("Uncorrected" = "#E64B35", "Corrected" = "#4DBBD5"),
        name   = NULL
      ) +
      labs(
        title    = sprintf("%s \u2014 XCO\u2082 Distribution (observed pixels)", sid),
        subtitle = sprintf("Mean: %.2f \u2192 %.2f ppm  (\u0394 = %.2f ppm)",
                           mean(obs_raw, na.rm = TRUE),
                           mean(obs_cor, na.rm = TRUE),
                           mean(obs_cor, na.rm = TRUE) - mean(obs_raw, na.rm = TRUE)),
        x = "XCO\u2082 (ppm)", y = "Frequency"
      ) +
      theme_panel()
    ggsave(file.path(FIGS_DIR, sprintf("%s_histogram.png", sid)),
           p_hist, width = 8, height = 6, dpi = 300, bg = "white")
    message(sprintf("  [%s] Saved histogram.", sid))
  }

  # --- Scatter: observed pixels only ---
  df_scatter <- df_scene %>%
    filter(!xco2_interpolated & !xco2_corrected_interpolated,
           is.finite(xco2), is.finite(xco2_corrected))

  if (nrow(df_scatter) > 10 && nrow(df_scatter) <= 10000) {
    p_scatter <- ggplot(df_scatter, aes(x = xco2, y = xco2_corrected)) +
      geom_point(alpha = 0.3, size = 0.8, color = "#4DBBD5") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                  color = "red", linewidth = 0.8) +
      geom_smooth(method = "lm", se = TRUE,
                  color = "black", linewidth = 1) +
      labs(
        title    = sprintf("%s \u2014 Correction Relationship", sid),
        subtitle = "Red: 1:1 line | Black: regression fit | Observed pixels only",
        x = "Uncorrected XCO\u2082 (ppm)",
        y = "Corrected XCO\u2082 (ppm)"
      ) +
      theme_panel()
    ggsave(file.path(FIGS_DIR, sprintf("%s_scatter.png", sid)),
           p_scatter, width = 8, height = 7, dpi = 300, bg = "white")
    message(sprintf("  [%s] Saved scatter plot.", sid))
  }
}


# =============================================================================
# SECTION 7 — Final Summary
# =============================================================================

message("\n", paste(rep("=", 70), collapse = ""))
message("PIPELINE COMPLETE")
message(paste(rep("=", 70), collapse = ""))
message(sprintf("Output directory: %s", OUTPUT_DIR))
message(sprintf("Figures:          %s", FIGS_DIR))
message("\nGenerated files:")
message("  unified_database.csv            — co-located pixel data (all scenes)")
message("  extraction_metadata.csv         — per-scene metadata")
message("  pooled_dredge_results.csv       — AICc model selection table")
message("  best_model_coefficients.csv     — final model coefficients")
message("  correction_statistics.csv       — per-scene correction metrics")
message("  model_diagnostics.png           — lm residual diagnostics")
message("  figures/<scene>_four_panel.png  — VIIRS AOD | MSAT Albedo | XCO2 | XCO2 corrected")
message("  figures/<scene>_bias_map.png    — spatial map of removed bias")
message("  figures/<scene>_histogram.png   — XCO2 distribution shift (observed only)")
message("  figures/<scene>_scatter.png     — pixel-level correction relationship")
message("\nAnalysis complete.")