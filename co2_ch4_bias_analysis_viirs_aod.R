library(ncdf4)
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mgcv)   # for spline smoothing via gam

output_dir <- "./results_bias_analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

VIIRS_DATA_DIR  <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/viirs_data2"
PYTHON_SCRIPT_PATH <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/download_viirs_deep_blue2.py"

if (!dir.exists(VIIRS_DATA_DIR)) dir.create(VIIRS_DATA_DIR, recursive = TRUE)

# ---- ALL SCENES TO PROCESS ----
scene_files <- c(
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t349/20250516/c0BC915D0/p3838/MSAT_L2_c0BC915D0_p3838_v04001006_20250516T212530Z_212602Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t348/20250430/c0A9D15C0/p3665/MSAT_L2_c0A9D15C0_p3665_v04001006_20250430T132730Z_132802Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t361/20250320/c09041690/p1403/MSAT_L2_c09041690_p1403_v03004000_20250320T063412Z_063444Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t380/20250319/c08F617C0/p1524/MSAT_L2_c08F617C0_p1524_v03005001_20250319T080703Z_080735Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t381/20250426/c0A4217D0/p1618/MSAT_L2_c0A4217D0_p1618_v03005001_20250426T144316Z_144348Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t25/20241119/c04350190/p898/MSAT_L2_c04350190_p898_v03002000_20241119T122526Z_122558Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t26/20241124/c049501A0/p748/MSAT_L2_c049501A0_p748_v03002000_20241124T143538Z_143610Z.nc"
)

# Aerosol_Type_Land_Ocean integer codes from VIIRS Deep Blue product
AEROSOL_TYPE_LABELS <- c(
  "0" = "dust",
  "1" = "smoke",
  "2" = "high alt. smoke",
  "3" = "pyrocumulonimbus",
  "4" = "non-smoke fine",
  "5" = "mixed",
  "6" = "background",
  "7" = "fine dominated"
)

# Short ID from path, e.g. "t349_c0BC915D0_p3838_20250516" for titles/filenames
scene_id_fn <- function(path) {
  base  <- tools::file_path_sans_ext(basename(path))
  parts <- strsplit(base, "_")[[1]]
  coll  <- parts[3]; pass <- parts[4]
  date  <- sub("T.*", "", parts[6])
  track <- regmatches(path, regexpr("t[0-9]+", path))
  track <- if (length(track) == 0) "" else track
  paste(track, coll, pass, date, sep = "_")
}

# =============================================================================
# VIIRS ACQUISITION HELPERS
# (adapted from MSAT_VIIRS_model_pipeline.R)
# =============================================================================

# Extracts bounding box, date string, and median UTC time from an MSAT NetCDF.
get_bbox_and_time <- function(msat_file) {
  nc <- try(nc_open(msat_file), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)
  lon       <- try(ncvar_get(nc, "geolocation/longitude"), silent = TRUE)
  lat       <- try(ncvar_get(nc, "geolocation/latitude"),  silent = TRUE)
  time_vals <- try(ncvar_get(nc, "geolocation/time"),       silent = TRUE)
  nc_close(nc)
  if (inherits(lon, "try-error") || inherits(lat, "try-error") || inherits(time_vals, "try-error")) return(NULL)
  origin   <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
  mid_time <- origin + median(time_vals, na.rm = TRUE)
  list(
    bbox = c(min(lon, na.rm=TRUE), min(lat, na.rm=TRUE), max(lon, na.rm=TRUE), max(lat, na.rm=TRUE)),
    date = as.character(as.Date(mid_time)),
    time = mid_time
  )
}

# Searches locally first, downloads via Python script if not found.
find_or_download_viirs <- function(msat_file) {
  scene_info <- get_bbox_and_time(msat_file)
  if (is.null(scene_info)) return(NA_character_)

  bbox  <- scene_info$bbox
  lon1  <- ifelse(bbox[1] > 180, bbox[1] - 360, bbox[1])
  lon2  <- ifelse(bbox[3] > 180, bbox[3] - 360, bbox[3])
  bbox_clean <- c(lon1, bbox[2], lon2, bbox[4])

  date_obj       <- as.Date(scene_info$date)
  year_doy       <- format(date_obj, "%Y%j")
  search_pattern <- sprintf("AERDB_L2_VIIRS_(SNPP|NOAA20|NOAA21).*A%s.*\\.(nc|nc4)$", year_doy)

  local_files <- list.files(VIIRS_DATA_DIR, pattern = search_pattern, full.names = TRUE, recursive = TRUE)
  if (length(local_files) > 0) {
    cat(sprintf(" -> Found %d existing VIIRS granule(s) for %s\n", length(local_files), scene_info$date))
    return(local_files[1])
  }

  cat(sprintf(" -> No local VIIRS file for %s. Downloading...\n", scene_info$date))
  if (!file.exists(PYTHON_SCRIPT_PATH)) {
    warning(sprintf("Downloader script not found at: %s", PYTHON_SCRIPT_PATH))
    return(NA_character_)
  }

  cmd_args <- c(
    PYTHON_SCRIPT_PATH,
    "--bbox",      as.character(bbox_clean[1]), as.character(bbox_clean[2]),
                   as.character(bbox_clean[3]), as.character(bbox_clean[4]),
    "--date",      scene_info$date,
    "--outdir",    VIIRS_DATA_DIR,
    "--platforms", "NOAA20,NOAA21"
  )

  res <- tryCatch(
    system2("python3", args = cmd_args, stdout = TRUE, stderr = TRUE),
    error = function(e) { warning("Download failed: ", conditionMessage(e)); NULL }
  )
  if (!is.null(res)) cat(paste(res, collapse = "\n"), "\n")

  downloaded <- list.files(VIIRS_DATA_DIR, pattern = search_pattern, full.names = TRUE, recursive = TRUE)
  if (length(downloaded) > 0) return(downloaded[1])

  cat(" -> VIIRS file could not be found or downloaded.\n")
  NA_character_
}

# Parses the UTC acquisition time embedded in a VIIRS granule filename (AYYYYDDD.HHMM).
get_viirs_granule_time <- function(fname) {
  m <- regmatches(fname, regexpr("A(\\d{7})\\.(\\d{4})", fname))
  if (length(m) == 0) return(as.POSIXct(NA, tz = "UTC"))
  yyyydoy <- substr(m, 2, 8)
  hhmm    <- substr(m, 10, 13)
  year <- as.integer(substr(yyyydoy, 1, 4))
  doy  <- as.integer(substr(yyyydoy, 5, 7))
  gdate <- as.Date(doy - 1, origin = paste0(year, "-01-01"))
  hh <- as.integer(substr(hhmm, 1, 2))
  mm <- as.integer(substr(hhmm, 3, 4))
  as.POSIXct(sprintf("%s %02d:%02d:00", format(gdate, "%Y-%m-%d"), hh, mm), tz = "UTC")
}

# Extracts mean AOD and top-3 aerosol types at the given MSAT pixel coordinates.
# Limits VIIRS granules to those within ±max_mins of the MSAT median overpass time.
#
# Fix: VIIRS Deep Blue is a curved swath product. The old approach of mapping the
# 2D array to a regular grid via min/max extents places values at wrong geographic
# positions. Instead, we read actual lon/lat per pixel, crop to the MSAT scene
# bounding box, then rasterize from those true coordinates before extracting.
extract_viirs_at_msat <- function(msat_file, msat_lon, msat_lat, max_mins = 500) {
  scene_info <- get_bbox_and_time(msat_file)
  if (is.null(scene_info)) return(list(mean_aod = NA_real_, top3_types = "no scene info"))

  scene_date <- as.Date(scene_info$date)
  msat_time  <- scene_info$time

  all_viirs <- list.files(VIIRS_DATA_DIR, pattern = "\\.(nc|nc4)$", full.names = TRUE)

  # Filter by calendar date
  viirs_day <- Filter(function(f) {
    fn <- basename(f)
    m  <- regmatches(fn, regexpr("A(\\d{7})", fn))
    if (length(m) == 0) return(FALSE)
    yyyydoy <- substr(m, 2, 8)
    yr  <- as.integer(substr(yyyydoy, 1, 4))
    doy <- as.integer(substr(yyyydoy, 5, 7))
    gdate <- as.Date(doy - 1, origin = paste0(yr, "-01-01"))
    !is.na(gdate) && gdate == scene_date
  }, all_viirs)

  # Filter to ±90 min; fall back to full day if none qualify
  viirs_close <- Filter(function(f) {
    gt <- get_viirs_granule_time(basename(f))
    if (is.na(gt)) return(FALSE)
    abs(as.numeric(difftime(gt, msat_time, units = "mins"))) <= max_mins
  }, viirs_day)

  if (length(viirs_close) == 0) {
    cat(sprintf("  -> No VIIRS granules within ±%d min; using all %d same-day granules.\n",
                max_mins, length(viirs_day)))
    viirs_close <- viirs_day
  }
  if (length(viirs_close) == 0) return(list(mean_aod = NA_real_, top3_types = "no VIIRS granules"))

  cat(sprintf("  -> Using %d VIIRS granule(s).\n", length(viirs_close)))

  # MSAT scene bounding box with 0.5° buffer for safe overlap
  bbox_buf <- c(
    min(msat_lon, na.rm = TRUE) - 0.5,
    min(msat_lat, na.rm = TRUE) - 0.5,
    max(msat_lon, na.rm = TRUE) + 0.5,
    max(msat_lat, na.rm = TRUE) + 0.5
  )

  # Read raw VIIRS pixels from each granule and crop to scene bbox.
  # Also read Angstrom exponent for Angstrom power-law scaling to 1606 nm.
  viirs_rows <- lapply(viirs_close, function(f) {
    nc     <- nc_open(f)
    lon_v  <- as.vector(ncvar_get(nc, "Longitude"))
    lat_v  <- as.vector(ncvar_get(nc, "Latitude"))
    aod_v  <- as.vector(ncvar_get(nc, "Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate"))
    ae_v   <- as.vector(ncvar_get(nc, "Angstrom_Exponent_Land_Ocean_Best_Estimate"))
    type_v <- as.vector(ncvar_get(nc, "Aerosol_Type_Land_Ocean"))
    aod_fill  <- ncatt_get(nc, "Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate", "_FillValue")
    ae_fill   <- ncatt_get(nc, "Angstrom_Exponent_Land_Ocean_Best_Estimate",             "_FillValue")
    type_fill <- ncatt_get(nc, "Aerosol_Type_Land_Ocean",                                "_FillValue")
    nc_close(nc)
    if (aod_fill$hasatt)  aod_v[aod_v   == aod_fill$value]  <- NA
    if (ae_fill$hasatt)   ae_v[ae_v     == ae_fill$value]   <- NA
    if (type_fill$hasatt) type_v[type_v == type_fill$value] <- NA
    in_box <- is.finite(lon_v) & is.finite(lat_v) &
              lon_v >= bbox_buf[1] & lon_v <= bbox_buf[3] &
              lat_v >= bbox_buf[2] & lat_v <= bbox_buf[4]
    if (sum(in_box) == 0) return(NULL)
    data.frame(lon  = lon_v[in_box],  lat  = lat_v[in_box],
               aod  = aod_v[in_box],  ae   = ae_v[in_box],
               type = type_v[in_box])
  })

  viirs_df <- do.call(rbind, Filter(Negate(is.null), viirs_rows))
  if (is.null(viirs_df) || nrow(viirs_df) == 0)
    return(list(type_stats = NULL, n_viirs_px = 0L))

  cat(sprintf("  -> %d VIIRS pixels in scene extent.\n", nrow(viirs_df)))

  # Per-type statistics for the top-3 aerosol types.
  # AOD at 1606 nm via Angstrom power law: AOD(λ) = AOD(550) * (λ/550)^(-α)
  SCALE_1606 <- (1606 / 550)           # wavelength ratio; exponent applied per pixel

  df_valid <- viirs_df[!is.na(viirs_df$type) & viirs_df$type >= 0 &
                         is.finite(viirs_df$aod) & is.finite(viirs_df$ae), ]

  if (nrow(df_valid) == 0) return(list(type_stats = NULL, n_viirs_px = nrow(viirs_df)))

  # Scale AOD per pixel before averaging (preserves the nonlinear scaling)
  df_valid$aod1606 <- df_valid$aod * SCALE_1606^(-df_valid$ae)

  type_key     <- as.character(round(df_valid$type))
  type_counts  <- sort(table(type_key), decreasing = TRUE)
  n_top        <- min(3L, length(type_counts))
  top_keys     <- names(type_counts)[seq_len(n_top)]
  top_labels   <- AEROSOL_TYPE_LABELS[top_keys]
  top_labels[is.na(top_labels)] <- paste0("type ", top_keys[is.na(top_labels)])
  top_pct      <- round(as.integer(type_counts[top_keys]) / nrow(df_valid) * 100)

  mean_aod550  <- sapply(top_keys, function(k) mean(df_valid$aod[type_key   == k], na.rm = TRUE))
  mean_ae      <- sapply(top_keys, function(k) mean(df_valid$ae[type_key    == k], na.rm = TRUE))
  mean_aod1606 <- sapply(top_keys, function(k) mean(df_valid$aod1606[type_key == k], na.rm = TRUE))

  type_stats <- data.frame(
    rank        = seq_len(n_top),
    label       = unname(top_labels),
    pct         = top_pct,
    mean_aod550 = mean_aod550,
    mean_ae     = mean_ae,
    mean_aod1606 = mean_aod1606,
    stringsAsFactors = FALSE
  )

  list(type_stats = type_stats, n_viirs_px = nrow(viirs_df))
}

# =============================================================================
# PER-SCENE PROCESSING
# =============================================================================
process_scene <- function(scene_file) {
  sid <- scene_id_fn(scene_file)
  message("Processing ", sid, " ...")

  nc <- nc_open(scene_file)
  on.exit(nc_close(nc), add = TRUE)

  get_flat <- function(name) {
    v    <- ncvar_get(nc, name)
    fill <- ncatt_get(nc, name, "_FillValue")
    if (fill$hasatt) v[v == fill$value] <- NA
    as.vector(v)
  }

  lon     <- get_flat("geolocation/longitude")
  lat     <- get_flat("geolocation/latitude")
  ch4_vcd <- get_flat("co2proxy_fit_diagnostics/retrieved_ch4_vcd")
  co2_vcd <- get_flat("co2proxy_fit_diagnostics/retrieved_co2_vcd")
  albedo  <- get_flat("co2proxy_fit_diagnostics/retrieved_albedo_1606nm")
  qflag   <- get_flat("product_co2proxy/main_quality_flag")  # 0 = good

  keep <- qflag == 0 & is.finite(albedo) & is.finite(ch4_vcd) & is.finite(co2_vcd)
  if (sum(keep) < 50) {
    warning("Scene ", sid, " has only ", sum(keep), " good pixels; skipping.")
    return(invisible(NULL))
  }

  # ---- VIIRS ACQUISITION & EXTRACTION ----
  cat(sprintf("  Acquiring VIIRS for %s...\n", sid))
  viirs_file <- tryCatch(find_or_download_viirs(scene_file), error = function(e) NA_character_)

  viirs_info <- list(type_stats = NULL, n_viirs_px = 0L)
  if (!is.na(viirs_file)) {
    viirs_info <- tryCatch(
      extract_viirs_at_msat(scene_file, lon[keep], lat[keep]),
      error = function(e) {
        message("  VIIRS extraction failed: ", conditionMessage(e))
        list(type_stats = NULL, n_viirs_px = 0L)
      }
    )
  }

  # Build annotation: per-type mean AOD at 550 nm and scaled to 1606 nm
  ts <- viirs_info$type_stats
  if (!is.null(ts) && nrow(ts) > 0) {
    header <- "Top aerosols  AOD550  AOD1606 % of pixels"
    rows   <- sprintf("%-18s %5.3f   %6.4f  (%d%%)",
                      ts$label, ts$mean_aod550, ts$mean_aod1606, ts$pct)
    annot_text <- paste(c(header, rows), collapse = "\n")
  } else {
    annot_text <- "VIIRS data unavailable"
  }

  # ---- NORMALIZED PERCENT BIAS ----
  pct_bias <- function(x) (x - mean(x)) / mean(x) * 100
  df <- data.frame(
    albedo = albedo[keep],
    CH4    = pct_bias(ch4_vcd[keep]),
    CO2    = pct_bias(co2_vcd[keep])
  ) %>%
    pivot_longer(c(CH4, CO2), names_to = "gas", values_to = "bias_pct")

  # ---- PLOT ----
  p <- ggplot(df, aes(albedo, bias_pct, color = gas, fill = gas)) +
    geom_point(alpha = 0.15, size = 0.6) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, linewidth = 1) +
    annotate(
      "label",
      x = Inf, y = -Inf,
      label    = annot_text,
      hjust    = 1, vjust = 0,
      size   = 2.8,
      fill   = "white",
      alpha  = 0.85,
      family = "mono"
    ) +
    labs(
      title    = paste0("Normalized % Bias vs Retrieved Albedo (1606 nm) — ", sid),
      subtitle = sprintf("N = %d good pixels | retrieved VCDs", sum(keep)),
      x        = "Retrieved albedo (CO2 band, 1606 nm)",
      y        = "Normalized bias (%)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  out_path <- file.path(output_dir, paste0("bias_vs_albedo_", sid, ".png"))
  ggsave(out_path, p, width = 9, height = 6, dpi = 200)
  message("  saved ", out_path)
  invisible(p)
}

# ---- LOOP OVER ALL SCENES ----
for (f in scene_files) {
  tryCatch(
    process_scene(f),
    error = function(e) message("FAILED ", basename(f), ": ", conditionMessage(e))
  )
}

message("Done. Plots are in ", normalizePath(output_dir))
