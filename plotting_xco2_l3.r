# Plotting for CO2 and VIIRS 

# ===============================

# ===============================
# MethaneSAT XCO2 / dP / Albedo — Minimal Plot Script
#
# Plots exactly three maps per scene:
#   1. Binned XCO2 (O2-proxy, ppm)
#   2. Delta-pressure (dP, hPa)
#   3. Albedo at 1606 nm
#
# All reading, binning, smoothing, and plotting methods are taken
# verbatim from CO2_flags_02172026.r — nothing has been changed.
# Only the output section is reduced to the three maps above.
# ===============================

suppressPackageStartupMessages({
  library(ncdf4)
  library(MASS)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(tibble)
  library(tools)
  library(tidyr)
  library(terra)
  library(grid)
  library(gridExtra)
})

# ===============================
# ---- USER SETTINGS ----
# ===============================

selected_targets <- TRUE
starting_id <- NA

if (selected_targets == TRUE) {
  runs <- tibble::tribble(
    ~TARGET_ID, ~SCENE_DATE,
    56, "20241223",
    349, "20250516",
    349, "20241230"
  )
} else {
  all_ids    <- list.files("/mnt/gcs/msat-prod-data-methanesat-level3/regrid")
  target_ids <- as.integer(sub("^t", "", all_ids))
  target_ids <- target_ids[which(target_ids == starting_id):length(target_ids)]
  runs_list  <- list()
  for (tid in target_ids) {
    scene_dates <- list.files(file.path("/mnt/gcs/msat-prod-data-methanesat-level3/regrid",
                                        paste0("t", tid)))
    for (sd in scene_dates)
      runs_list[[length(runs_list) + 1]] <- data.frame(TARGET_ID = tid, SCENE_DATE = sd,
                                                        stringsAsFactors = FALSE)
  }
  runs <- as_tibble(do.call(rbind, runs_list))
}

ROOT_DIR <- "/mnt/gcs/msat-prod-data-methanesat-level3/regrid"
OUT_ROOT <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/scripts/L3_results"

# Binning
AGG_R        <- 10
AGG_C        <- 10
AGG_MIN_FRAC <- 0.1

# ---- Pre-compute color palette ----
ssec.256 <- readLines("/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/scripts/ssec.256.txt")

# Upsampling
UPSAMPLE_ENABLE <- TRUE
UPSAMPLE_FACTOR <- 4

FONT_SCALE <- 1.7

# ===============================
# ---- Helpers (unchanged from original)
# ===============================

find_l3_regrid_nc <- function(root_dir, target_id, scene_date,
                               grid_dir = "45m", pattern = "\\.nc$",
                               prefer_msat_prefix = TRUE,
                               keep_only_highest_p = TRUE) {
  base_dir <- file.path(root_dir, paste0("t", target_id), scene_date)
  if (!dir.exists(base_dir)) return(character(0))
  candidates <- list.files(base_dir, recursive = TRUE, full.names = TRUE, pattern = pattern)
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
      pmax <- max(pp[is.finite(pp)])
      candidates <- candidates[pp == pmax]
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
  nr <- nrow(m); nc <- ncol(m)
  block_n <- row_block_size * col_block_size

  if (nr %% row_block_size == 0 && nc %% col_block_size == 0) {
    nr_out <- nr / row_block_size
    nc_out <- nc / col_block_size
    a        <- array(m, dim = c(row_block_size, nr_out, col_block_size, nc_out))
    out_mean <- apply(a, c(2, 4), mean, na.rm = TRUE)
    a_valid  <- array(is.finite(m), dim = c(row_block_size, nr_out, col_block_size, nc_out))
    out_nval <- apply(a_valid, c(2, 4), sum)
    out_mean[out_nval < min_frac * block_n] <- NA_real_
    return(out_mean)
  }

  r   <- terra::rast(m)
  agg <- terra::aggregate(r, fact = c(row_block_size, col_block_size),
                           fun = "mean", na.rm = TRUE)
  result    <- terra::as.matrix(agg, wide = TRUE)
  r_valid   <- terra::rast(is.finite(m) * 1L)
  agg_nval  <- terra::aggregate(r_valid, fact = c(row_block_size, col_block_size),
                                 fun = "sum", na.rm = TRUE)
  nval_mat  <- terra::as.matrix(agg_nval, wide = TRUE)
  result[nval_mat < min_frac * block_n] <- NA_real_
  result
}

mat_to_df <- function(m, name = "val") {
  if (!is.matrix(m)) return(tibble())
  nr <- nrow(m); nc <- ncol(m)
  ij <- expand.grid(i = seq_len(nr), j = seq_len(nc))
  tibble(i = ij$i, j = ij$j, !!name := as.vector(m))
}

pct_limits <- function(x, lo = 0.05, hi = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NULL)
  as.numeric(stats::quantile(x, probs = c(lo, hi), na.rm = TRUE))
}

add_fill_scale_ssec <- function(lims, colors = ssec.256) {
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

upsample_matrix <- function(mat, factor = 4, method = "bilinear") {
  if (!is.matrix(mat)) return(NULL)
  nr <- nrow(mat); nc <- ncol(mat)
  r <- terra::rast(mat)
  terra::ext(r) <- terra::ext(1, nr, 1, nc)
  r_t <- terra::rast(nrows = nr * factor, ncols = nc * factor)
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
  m_up <- upsample_matrix(m, factor = as.integer(factor), method = method)
  if (!is.matrix(m_up)) return(list(mat = m, factor = 1L))
  list(mat = m_up, factor = as.integer(factor))
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
    scale_x_continuous(expand = c(0,0), breaks = x_breaks, labels = fmt_deg(lon_by_col[x0], digits)),
    scale_y_continuous(expand = c(0,0), breaks = y_breaks, labels = fmt_deg(lat_by_row[y0], digits))
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

plot_caption <- function(TARGET_ID, SCENE_DATE, l2_pp, extra = NULL) {
  cap <- paste0("Target t", TARGET_ID, " | Scene ", SCENE_DATE,
                " | File: ", basename(l2_pp))
  if (!is.null(extra) && nzchar(extra)) cap <- paste0(cap, " | ", extra)
  cap
}

plot_matrix_map_originalgrid <- function(mat, nm, lon_b, lat_b,
                                          fill_lab = "", use_ssec = FALSE,
                                          caption_text = NULL,
                                          upsample_for_plot = FALSE,
                                          upsample_kind = c("cont", "near")) {
  if (!is.matrix(mat) || !is.matrix(lon_b) || !is.matrix(lat_b)) return(NULL)
  upsample_kind <- match.arg(upsample_kind)
  pp       <- if (isTRUE(upsample_for_plot)) prep_plot_matrix(mat, kind = upsample_kind) else list(mat = mat, factor = 1L)
  mat_plot <- pp$mat
  d        <- mat_to_df(mat_plot, "val")
  fill_sc  <- if (use_ssec) add_fill_scale_ssec(pct_limits(d$val)) else add_fill_scale(pct_limits(d$val))
  ggplot(d, aes(x = j, y = i, fill = val)) +
    geom_raster(interpolate = FALSE) +
    coord_equal(expand = FALSE) +
    fill_sc +
    add_latlon_axes_plot(lon_b, lat_b,
                         plot_nr = nrow(mat_plot), plot_nc = ncol(mat_plot),
                         n_ticks = 5, digits = 2) +
    labs(title   = paste0("RAW: ", nm),
         fill    = fill_lab,
         x       = "Approx. Longitude",
         y       = "Approx. Latitude",
         caption = if (!is.null(caption_text)) caption_text else "") +
    theme_msat_map()
}

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

  make_empty_like <- function(lat_len, lon_len) matrix(NA_real_, nrow = lat_len, ncol = lon_len)

  lat_vec <- tryGet_any(c("lat", "latitude"))
  lon_vec <- tryGet_any(c("lon", "longitude"))
  if (!is.numeric(lat_vec) || !is.numeric(lon_vec) ||
      length(lat_vec) < 2 || length(lon_vec) < 2)
    stop("Could not read 1D lat/lon from: ", nc_path)

  lat_mat <- matrix(rep(lat_vec, times = length(lon_vec)),
                    nrow = length(lat_vec), ncol = length(lon_vec), byrow = FALSE)
  lon_mat <- matrix(rep(lon_vec, each  = length(lat_vec)),
                    nrow = length(lat_vec), ncol = length(lon_vec), byrow = FALSE)
  empty2d <- make_empty_like(length(lat_vec), length(lon_vec))

  dp      <- to_latlon(tryGet_any(c("o2dp_fit_diagnostics/bias_corrected_delta_pressure",
                                     "o2dp_fit_diagnostics/delta_pressure", "delta_pressure"),
                                   default = empty2d))
  alb1606 <- to_latlon(tryGet_any(c("apriori_data/albedo_ch4band", "albedo_ch4band"),
                                   default = empty2d))
  co2_vcd_bc <- to_latlon(tryGet_any(c("co2proxy_fit_diagnostics/bias_corrected_co2_vcd",
                                        "bias_corrected_co2_vcd"), default = empty2d))
  ch4_vcd_bc <- to_latlon(tryGet_any(c("co2proxy_fit_diagnostics/bias_corrected_ch4_vcd",
                                        "bias_corrected_ch4_vcd"), default = empty2d))
  o2_vcd_bc  <- to_latlon(tryGet_any(c("o2dp_fit_diagnostics/bias_corrected_o2_vcd",
                                        "bias_corrected_o2_vcd"), default = empty2d))
  xch4_native <- to_latlon(tryGet_any(c("xch4"), default = empty2d))

  xco2_o2 <- empty2d
  ok_o2   <- is.finite(o2_vcd_bc) & (o2_vcd_bc != 0) &
              is.finite(co2_vcd_bc) & is.finite(ch4_vcd_bc)
  xco2_o2[ok_o2] <- (co2_vcd_bc[ok_o2] / o2_vcd_bc[ok_o2]) * 0.2095
  xco2_o2_ppm    <- xco2_o2 * 1e6

  list(
    lat = lat_mat, lon = lon_mat,
    dp = dp, alb1606 = alb1606,
    xco2_o2_ppm = xco2_o2_ppm,
    xch4_co2_ppb = xch4_native    # used only to derive the valid-pixel mask
  )
}

# ===============================
# ---- Main loop
# ===============================

for (k in seq_len(nrow(runs))) {
  TARGET_ID  <- runs$TARGET_ID[k]
  SCENE_DATE <- runs$SCENE_DATE[k]

  OUT_DIR <- file.path(OUT_ROOT, "figures_xco2_dp_alb")
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

  pdf_pat   <- paste0("^maps_xco2_dp_alb_t", TARGET_ID, "_", SCENE_DATE, "_.*\\.pdf$")
  scene_pdf <- list.files(OUT_DIR, pattern = pdf_pat, full.names = TRUE)
  if (length(scene_pdf) > 0) {
    message("Skipping (already processed): t", TARGET_ID, " ", SCENE_DATE)
    next
  }

  nc_files <- find_l3_regrid_nc(ROOT_DIR, TARGET_ID, SCENE_DATE, grid_dir = "45m")
  if (length(nc_files) == 0) {
    message("No L3 regrid .nc found for t", TARGET_ID, " ", SCENE_DATE, " (skipping)")
    next
  }
  nc_files <- pick_latest_mtime(nc_files)
  message("Using NetCDF: ", nc_files)

  for (l2_pp in nc_files) {
    message("\n=== Processing: ", basename(l2_pp), " ===")

    l3 <- read_msat_l3_regrid(l2_pp)

    lat         <- l3$lat
    lon         <- l3$lon
    dp          <- l3$dp
    alb1606     <- l3$alb1606
    xco2        <- l3$xco2_o2_ppm
    xch4_co2    <- l3$xch4_co2_ppb

    # Valid-pixel mask driven by xch4_co2 (same logic as original)
    mask_raw <- ifelse(is.finite(xch4_co2), 1, 0)

    # Bin all fields
    xco2_b    <- aggregate_matrix(xco2)
    dp_b      <- aggregate_matrix(dp)
    alb1606_b <- aggregate_matrix(alb1606)
    lat_b     <- aggregate_matrix(lat)
    lon_b     <- aggregate_matrix(lon)
    mask_b    <- aggregate_matrix(mask_raw)
    mask_b    <- ifelse(is.finite(mask_b) & mask_b >= AGG_MIN_FRAC, 1, 0)

    # Apply mask
    xco2_b[mask_b    == 0] <- NA_real_
    dp_b[mask_b      == 0] <- NA_real_
    alb1606_b[mask_b == 0] <- NA_real_

    cap <- plot_caption(TARGET_ID, SCENE_DATE, l2_pp)

    # ---- Three maps ----
    p_xco2 <- plot_matrix_map_originalgrid(
      xco2_b, "Binned XCO2 (O2-proxy)",
      lon_b = lon_b, lat_b = lat_b,
      fill_lab = "ppm", use_ssec = TRUE,
      upsample_for_plot = TRUE,
      caption_text = cap
    )

    p_dp <- plot_matrix_map_originalgrid(
      dp_b, "Delta-pressure (dP)",
      lon_b = lon_b, lat_b = lat_b,
      fill_lab = "hPa", use_ssec = FALSE,
      upsample_for_plot = TRUE,
      caption_text = cap
    )

    p_alb <- plot_matrix_map_originalgrid(
      alb1606_b, "Albedo at 1606 nm",
      lon_b = lon_b, lat_b = lat_b,
      fill_lab = "", use_ssec = FALSE,
      upsample_for_plot = TRUE,
      caption_text = cap
    )

    # ---- Write one-page PDF (1 row × 3 cols) ----
    pdf_path <- file.path(OUT_DIR,
      paste0("maps_xco2_dp_alb_t", TARGET_ID, "_", SCENE_DATE, "_",
             tools::file_path_sans_ext(basename(l2_pp)), ".pdf"))

    grDevices::pdf(pdf_path, width = 21, height = 8, onefile = TRUE)
    gridExtra::grid.arrange(p_xco2, p_dp, p_alb, nrow = 1L, ncol = 3L)
    grDevices::dev.off()

    message("Wrote: ", pdf_path)
  }
}

message("\nAll scenes complete.")


# ===============================
# ---- VIIRS Deep Blue AOD Section ----
#
# Reads VIIRS AERDB_L2_VIIRS_SNPP NetCDF4 files and plots:
#   AOD at 550 nm (Deep Blue, land best estimate)
#
# File format:
#   - Variable: Aerosol_Optical_Thickness_550_Land_Best_Estimate
#   - Grid: YDim × XDim (geographic, GCTP_GEO)
#   - Coordinate variables: Latitude, Longitude (2D arrays, YDim × XDim)
#   - Fill values: any value <= -999 (e.g. -999.9, -999.8, … -999)
#   - Valid range: [0, ~5] (physically meaningful AOD)
#
# Output: one PDF per file in figures_viirs_aod/
# ===============================

# ---- VIIRS: user settings ----
# Directory containing VIIRS AERDB_L2_VIIRS_SNPP*.nc files.
# Set to a single path or a vector of paths.
VIIRS_NC_DIR  <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/viirs_data2"
VIIRS_NC_GLOB <- "AERDB_L2_VIIRS_SNPP.*\\.nc$"   # regex passed to list.files()

# If you want to process only specific files, list them here; otherwise
# set to NULL to discover all matching files under VIIRS_NC_DIR.
VIIRS_NC_FILES_OVERRIDE <- NULL   # e.g. c("/path/to/file1.nc", "/path/to/file2.nc")

VIIRS_OUT_DIR <- file.path(OUT_ROOT, "figures_viirs_aod")

# ---- VIIRS: reader ----
read_viirs_aerdb <- function(nc_path) {
  # Returns a list with:
  #   $aod   – matrix [nrow=YDim, ncol=XDim], NA where fill/out-of-range
  #   $lat   – matrix [nrow=YDim, ncol=XDim], degrees N
  #   $lon   – matrix [nrow=YDim, ncol=XDim], degrees E
  #   $fname – basename of file (for labelling)

  stopifnot(file.exists(nc_path))
  nc <- ncdf4::nc_open(nc_path)
  on.exit(try(ncdf4::nc_close(nc), silent = TRUE), add = TRUE)

  safe_get <- function(name) {
    out <- try(ncdf4::ncvar_get(nc, name), silent = TRUE)
    if (inherits(out, "try-error")) NULL else out
  }

  # AOD variable
  aod_raw <- safe_get("Aerosol_Optical_Thickness_550_Land_Best_Estimate")
  if (is.null(aod_raw))
    stop("Could not read AOD variable from: ", nc_path)

  # Coordinate variables — AERDB files store them as 2-D arrays (YDim, XDim)
  # ncdf4 returns them transposed (XDim, YDim), so we aperm below.
  lat_raw <- safe_get("Latitude")
  lon_raw <- safe_get("Longitude")

  # Helper: transpose if needed so result is [YDim, XDim] = [nrow, ncol]
  fix_dims <- function(x, target_nrow, target_ncol) {
    if (is.null(x)) return(matrix(NA_real_, nrow = target_nrow, ncol = target_ncol))
    d <- dim(x)
    if (length(d) == 2) {
      if (d[1] == target_nrow && d[2] == target_ncol) return(as.matrix(x))
      if (d[1] == target_ncol && d[2] == target_nrow) return(t(as.matrix(x)))
    }
    # 1-D vectors (degenerate case): try to build 2-D from bounding box
    return(matrix(NA_real_, nrow = target_nrow, ncol = target_ncol))
  }

  # AOD is stored [XDim, YDim] by ncdf4 (Fortran order), transpose to [YDim, XDim]
  d_aod <- dim(aod_raw)
  if (length(d_aod) == 2) {
    # ncdf4 reads as (XDim=1199, YDim=500); transpose -> (500, 1199)
    aod_mat <- t(as.matrix(aod_raw))
  } else {
    aod_mat <- as.matrix(aod_raw)
  }
  nr <- nrow(aod_mat); nc_aod <- ncol(aod_mat)

  lat_mat <- fix_dims(lat_raw, nr, nc_aod)
  lon_mat <- fix_dims(lon_raw, nr, nc_aod)

  # If lat/lon are missing, reconstruct from bounding box metadata
  if (all(is.na(lat_mat)) || all(is.na(lon_mat))) {
    message("  [VIIRS] Coordinate variables missing; reconstructing from global attributes.")
    get_att <- function(nm) {
      v <- try(ncdf4::ncatt_get(nc, 0, nm)$value, silent = TRUE)
      if (inherits(v, "try-error") || is.null(v)) NA_real_ else as.numeric(v)
    }
    lon_min <- get_att("geospatial_lon_min")
    lon_max <- get_att("geospatial_lon_max")
    lat_min <- get_att("geospatial_lat_min")
    lat_max <- get_att("geospatial_lat_max")
    if (all(is.finite(c(lon_min, lon_max, lat_min, lat_max)))) {
      lon_vec <- seq(lon_min, lon_max, length.out = nc_aod)
      lat_vec <- seq(lat_max, lat_min, length.out = nr)   # N→S
      lon_mat <- matrix(rep(lon_vec, each = nr), nrow = nr, ncol = nc_aod)
      lat_mat <- matrix(rep(lat_vec, times = nc_aod), nrow = nr, ncol = nc_aod)
    }
  }

  # Mask fill values (all fill codes are <= -999)
  aod_mat[aod_mat <= -999] <- NA_real_
  # Mask physically implausible values
  aod_mat[aod_mat < -0.2 | aod_mat > 10] <- NA_real_

  list(aod = aod_mat, lat = lat_mat, lon = lon_mat,
       fname = basename(nc_path))
}

# ---- VIIRS: plotting helper (reuses existing infrastructure) ----
plot_viirs_aod <- function(viirs, caption_text = NULL) {
  # viirs is the list returned by read_viirs_aerdb()
  aod <- viirs$aod
  lat <- viirs$lat
  lon <- viirs$lon

  if (!is.matrix(aod) || all(is.na(aod))) {
    message("  [VIIRS] No valid AOD data to plot in: ", viirs$fname)
    return(NULL)
  }

  # Upsample for smoother display (reuses existing prep_plot_matrix)
  pp       <- prep_plot_matrix(aod, kind = "cont")
  aod_plot <- pp$mat
  d        <- mat_to_df(aod_plot, "val")

  # Use viridis "A" (magma) — works well for AOD (yellow=high, dark=low)
  lims     <- pct_limits(d$val, lo = 0.02, hi = 0.98)
  fill_sc  <- if (!is.null(lims) && all(is.finite(lims)) && lims[1] < lims[2])
    scale_fill_viridis_c(option = "A", direction = -1,
                         na.value = "transparent",
                         limits = lims, oob = scales::squish,
                         name = "AOD 550 nm")
  else
    scale_fill_viridis_c(option = "A", direction = -1,
                         na.value = "transparent",
                         name = "AOD 550 nm")

  ggplot(d, aes(x = j, y = i, fill = val)) +
    geom_raster(interpolate = FALSE) +
    coord_equal(expand = FALSE) +
    fill_sc +
    add_latlon_axes_plot(lon, lat,
                         plot_nr = nrow(aod_plot), plot_nc = ncol(aod_plot),
                         n_ticks = 5, digits = 1) +
    labs(
      title    = "VIIRS Deep Blue AOD at 550 nm (Land, Best Estimate)",
      subtitle = viirs$fname,
      x        = "Approx. Longitude",
      y        = "Approx. Latitude",
      caption  = if (!is.null(caption_text)) caption_text else
        paste0("Source: ", viirs$fname,
               " | Valid pixels: ",
               format(sum(is.finite(viirs$aod)), big.mark = ","), " / ",
               format(length(viirs$aod), big.mark = ","))
    ) +
    theme_msat_map()
}

# ---- VIIRS: main loop ----
if (!dir.exists(VIIRS_OUT_DIR)) dir.create(VIIRS_OUT_DIR, recursive = TRUE)

viirs_files <- if (!is.null(VIIRS_NC_FILES_OVERRIDE)) {
  VIIRS_NC_FILES_OVERRIDE
} else {
  list.files(VIIRS_NC_DIR, pattern = VIIRS_NC_GLOB,
             full.names = TRUE, recursive = FALSE)
}

if (length(viirs_files) == 0) {
  message("\n[VIIRS] No files found under: ", VIIRS_NC_DIR,
          "\n  → Set VIIRS_NC_DIR or VIIRS_NC_FILES_OVERRIDE and re-run.")
} else {
  message("\n[VIIRS] Processing ", length(viirs_files), " file(s)...")

  for (vf in viirs_files) {
    out_pdf <- file.path(VIIRS_OUT_DIR,
      paste0("viirs_aod_", tools::file_path_sans_ext(basename(vf)), ".pdf"))

    if (file.exists(out_pdf)) {
      message("  Skipping (already exists): ", basename(out_pdf))
      next
    }

    message("  Reading: ", basename(vf))
    viirs <- tryCatch(read_viirs_aerdb(vf), error = function(e) {
      message("  ERROR reading ", basename(vf), ": ", conditionMessage(e))
      NULL
    })
    if (is.null(viirs)) next

    n_valid <- sum(is.finite(viirs$aod))
    message("  Valid AOD pixels: ", format(n_valid, big.mark = ","),
            " / ", format(length(viirs$aod), big.mark = ","),
            "  (", round(100 * n_valid / length(viirs$aod), 1), "%)")

    p_aod <- plot_viirs_aod(viirs)
    if (is.null(p_aod)) next

    grDevices::pdf(out_pdf, width = 14, height = 7, onefile = TRUE)
    print(p_aod)
    grDevices::dev.off()

    message("  Wrote: ", out_pdf)
  }
}

message("\n[VIIRS] Done.")
