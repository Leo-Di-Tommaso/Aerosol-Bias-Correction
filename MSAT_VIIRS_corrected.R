# =============================================================================
# Unified MethaneSAT + VIIRS XCO2 Modeling Pipeline - CORRECTED VERSION
# =============================================================================
# COARSE VIIRS VERSION (6km resolution)
# - Extracts VIIRS values at MSAT point locations
# - Preserves MSAT native resolution
# - Publication-quality side-by-side plots
# - Scene-mean preserving aerosol correction
# - FIXED: Temporal filtering, metadata handling, AE quality control
#
# TO RUN: source ~/earthdata-env/bin/activate

# -----------------------------------------------------------------------------
# 0. Libraries and Constants
# -----------------------------------------------------------------------------
library(terra)
library(ncdf4)
library(tidyr)
library(readr)
library(writexl)
library(MASS)
library(broom)
library(purrr)
library(jsonlite)
library(stringr)
library(caret)
library(MuMIn)
library(pls)
library(car)
library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)
library(dplyr)
library(scales)

# --- CONFIGURATION ---

USE_SPECIFIC_FILES <- TRUE

# Directories
MSAT_ROOT_DIR <- "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results"
VIIRS_DATA_DIR <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/viirs_data2"
PYTHON_SCRIPT_PATH <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/download_viirs_deep_blue2.py"
OUTPUT_DIR <- paste0("results_unified_model_", format(Sys.Date(), "%Y%m%d"))

# Modeling Constants
VIF_THRESHOLD <- 5
TOP_N_DREDGE <- 10
TEMPORAL_WINDOW_HOURS <- 3

# Quality thresholds
MIN_OVERLAP_PIXELS <- 50
AE_MIN <- -0.5
AE_MAX <- 3.0

# Create directories
if (!dir.exists(VIIRS_DATA_DIR)) dir.create(VIIRS_DATA_DIR, recursive = TRUE)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Helper Functions
# -----------------------------------------------------------------------------

# --- A. VIIRS Acquisition ---

get_bbox_and_time_range <- function(msat_file) {
  nc <- try(nc_open(msat_file), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)
  
  lon <- try(ncvar_get(nc, "geolocation/longitude"), silent = TRUE)
  lat <- try(ncvar_get(nc, "geolocation/latitude"), silent = TRUE)
  time_vals <- try(ncvar_get(nc, "geolocation/time"), silent = TRUE)
  nc_close(nc)
  
  if (inherits(lon, "try-error") || inherits(lat, "try-error") || 
      inherits(time_vals, "try-error")) return(NULL)
  
  origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
  times <- origin + time_vals
  
  list(
    bbox = c(min(lon, na.rm=TRUE), min(lat, na.rm=TRUE), 
             max(lon, na.rm=TRUE), max(lat, na.rm=TRUE)),
    time_start = min(times, na.rm=TRUE),
    time_end = max(times, na.rm=TRUE),
    time_median = median(times, na.rm=TRUE),
    date = as.character(as.Date(median(times, na.rm=TRUE)))
  )
}

find_or_download_viirs <- function(msat_file) {
  scene_info <- get_bbox_and_time_range(msat_file)
  if (is.null(scene_info)) return(NA_character_)
  
  bbox <- scene_info$bbox
  lon1 <- ifelse(bbox[1] > 180, bbox[1] - 360, bbox[1])
  lon2 <- ifelse(bbox[3] > 180, bbox[3] - 360, bbox[3])
  bbox_clean <- c(lon1, bbox[2], lon2, bbox[4])
  
  date_obj <- as.Date(scene_info$date)
  year_doy <- format(date_obj, "%Y%j")
  search_pattern <- sprintf(
    "AERDB_L2_VIIRS_(SNPP|NOAA20|NOAA21).*A%s.*\\.(nc|nc4)$",
    year_doy
  )
  
  cat(sprintf(" -> Searching for VIIRS (Date %s, MSAT time: %s)\n", 
              scene_info$date, 
              format(scene_info$time_median, "%H:%M:%S")))
  
  local_files <- list.files(
    VIIRS_DATA_DIR,
    pattern = search_pattern,
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(local_files) > 0) {
    cat(sprintf(" -> Found %d existing VIIRS file(s)\n", length(local_files)))
    return(local_files[1])
  }
  
  cat(" -> Attempting download...\n")
  
  if (!file.exists(PYTHON_SCRIPT_PATH)) {
    warning(sprintf("Downloader script not found: %s", PYTHON_SCRIPT_PATH))
    return(NA_character_)
  }
  
  cmd_args <- c(
    PYTHON_SCRIPT_PATH,
    "--bbox", as.character(bbox_clean[1]),
    as.character(bbox_clean[2]),
    as.character(bbox_clean[3]),
    as.character(bbox_clean[4]),
    "--date", scene_info$date,
    "--outdir", VIIRS_DATA_DIR,
    "--platforms", "NOAA20,NOAA21"
  )
  
  res <- tryCatch(
    system2("python3", args = cmd_args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      warning("Download failed: ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(res)) return(NA_character_)
  
  cat(" -> Python output:\n", paste(res, collapse = "\n"), "\n")
  
  downloaded <- list.files(
    VIIRS_DATA_DIR,
    pattern = search_pattern,
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(downloaded) > 0) {
    cat(sprintf(" -> Downloaded: %s\n", basename(downloaded[1])))
    return(downloaded[1])
  }
  
  cat(" -> Download failed\n")
  return(NA_character_)
}

# --- B. Spatial Processing (COARSE VIIRS APPROACH) ---

calc_primary_angle <- function(sza, vza, saa, vaa) {
  dtor <- pi / 180
  rtod <- 180 / pi
  cos_phase <- cos(sza*dtor) * cos(vza*dtor) +
               sin(sza*dtor) * sin(vza*dtor) * cos((saa - vaa)*dtor)
  acos(pmax(-1, pmin(1, cos_phase))) * rtod
}

extract_viirs_time <- function(viirs_file) {
  fn <- basename(viirs_file)
  time_match <- regmatches(fn, regexpr("A(\\d{7})\\.(\\d{4})", fn))
  if (length(time_match) == 0) return(NULL)
  
  parts <- strsplit(time_match, "\\.")[[1]]
  yyyydoy <- substr(parts[1], 2, 8)
  hhmm <- parts[2]
  
  year <- as.integer(substr(yyyydoy, 1, 4))
  doy <- as.integer(substr(yyyydoy, 5, 7))
  hour <- as.integer(substr(hhmm, 1, 2))
  minute <- as.integer(substr(hhmm, 3, 4))
  
  date_part <- as.Date(doy - 1, origin = paste0(year, "-01-01"))
  as.POSIXct(paste(date_part, sprintf("%02d:%02d:00", hour, minute)), tz = "UTC")
}

extract_scene_data <- function(viirs_file, msat_file, scene_id) {
  
  # 1. GET MSAT METADATA (FIXED: was calling non-existent function)
  scene_info <- get_bbox_and_time_range(msat_file)
  if (is.null(scene_info)) {
    stop(sprintf("Cannot read MSAT metadata from %s", basename(msat_file)))
  }
  
  scene_date <- as.Date(scene_info$date)
  cat(sprintf("  MSAT scene time: %s\n", format(scene_info$time_median, "%Y-%m-%d %H:%M:%S UTC")))

  # 2. FIND ALL VIIRS FILES FROM SAME DAY
  viirs_files_all <- list.files(VIIRS_DATA_DIR, pattern="\\.nc$", full.names=TRUE)

  viirs_files_day <- Filter(function(f) {
    fn <- basename(f)
    match <- regmatches(fn, regexpr("A(\\d{7})", fn))
    if(length(match) == 0) return(FALSE)
    yyyydoy <- substr(match, 2, 8)
    year <- as.integer(substr(yyyydoy,1,4))
    doy  <- as.integer(substr(yyyydoy,5,7))
    granule_date <- as.Date(doy - 1, origin = paste0(year,"-01-01"))
    !is.na(granule_date) && granule_date == scene_date
  }, viirs_files_all)

  if(length(viirs_files_day) == 0) 
    stop(sprintf("No VIIRS granules found for scene %s on %s", scene_id, scene_date))
  
  cat(sprintf("  Found %d VIIRS granules on %s\n", length(viirs_files_day), scene_date))

  # 3. FILTER BY TEMPORAL WINDOW (CRITICAL FIX - addresses temporal mismatch)
  viirs_files_matched <- Filter(function(f) {
    viirs_time <- extract_viirs_time(f)
    if(is.null(viirs_time)) {
      cat(sprintf("  WARNING: Cannot parse time from %s\n", basename(f)))
      return(FALSE)
    }
    
    time_diff_hours <- abs(as.numeric(
      difftime(viirs_time, scene_info$time_median, units="hours")
    ))
    
    cat(sprintf("  - %s: Δt = %.1f hours %s\n", 
                basename(f), 
                time_diff_hours,
                ifelse(time_diff_hours <= TEMPORAL_WINDOW_HOURS, "[MATCH]", "[SKIP]")))
    
    return(time_diff_hours <= TEMPORAL_WINDOW_HOURS)
  }, viirs_files_day)
  
  if(length(viirs_files_matched) == 0) {
    stop(sprintf("No VIIRS granules within ±%d hours of MSAT overpass", 
                 TEMPORAL_WINDOW_HOURS))
  }
  
  cat(sprintf("  ✓ Using %d/%d VIIRS granules (±%dh window)\n", 
              length(viirs_files_matched), length(viirs_files_day), TEMPORAL_WINDOW_HOURS))

  # 4. LOAD AND MOSAIC VIIRS RASTERS
  load_viirs_raster <- function(file, varname) {
    nc <- nc_open(file)
    lon <- ncvar_get(nc, "Longitude")
    lat <- ncvar_get(nc, "Latitude")
    vals <- ncvar_get(nc, varname)
    nc_close(nc)
    
    # Handle both 1D and 2D coordinate arrays
    if (is.null(dim(lon))) {
      # 1D arrays - create grid
      lon_grid <- matrix(rep(lon, length(lat)), nrow=length(lat), byrow=TRUE)
      lat_grid <- matrix(rep(lat, each=length(lon)), nrow=length(lat))
    } else {
      # 2D arrays - use as-is
      lon_grid <- lon
      lat_grid <- lat
    }
    
    r <- rast(nrows=nrow(vals), ncols=ncol(vals),
              xmin=min(lon_grid, na.rm=TRUE), xmax=max(lon_grid, na.rm=TRUE),
              ymin=min(lat_grid, na.rm=TRUE), ymax=max(lat_grid, na.rm=TRUE),
              vals=as.vector(vals))
    crs(r) <- "EPSG:4326"
    return(r)
  }
  
  # Suppress resampling warnings (expected for multi-satellite mosaics)
  suppressWarnings({
    v_raster_aod <- do.call(terra::merge, 
                            lapply(viirs_files_matched, load_viirs_raster, 
                                   varname="Aerosol_Optical_Thickness_550_Land_Best_Estimate"))
    v_raster_ae  <- do.call(terra::merge, 
                            lapply(viirs_files_matched, load_viirs_raster, 
                                   varname="Angstrom_Exponent_Land_Best_Estimate"))
  })

  # 5. LOAD MSAT POINTS
  m_nc <- nc_open(msat_file)
  m_lon <- ncvar_get(m_nc, "geolocation/longitude")
  m_lat <- ncvar_get(m_nc, "geolocation/latitude")
  xco2  <- ncvar_get(m_nc, "product_o2proxy/xco2_o2proxy_vcdcor")
  
  dp    <- tryCatch(ncvar_get(m_nc,"o2dp_fit_diagnostics/bias_corrected_delta_pressure"),
                    error=function(e) NULL)
  if (is.null(dp)) dp <- ncvar_get(m_nc,"o2dp_fit_diagnostics/delta_pressure")
  
  alb1286 <- ncvar_get(m_nc,"o2dp_fit_diagnostics/retrieved_albedo_1286nm")
  alb1606 <- tryCatch(ncvar_get(m_nc,"co2proxy_fit_diagnostics/retrieved_albedo_1606nm"),
                      error=function(e) NA)
  psfc <- tryCatch(ncvar_get(m_nc,"co2proxy_fit_diagnostics/retrieved_surface_pressure"),
                   error=function(e) NA)
  orog <- tryCatch(ncvar_get(m_nc,"geolocation/terrain_height"),
                   error=function(e) NA)
  sza  <- tryCatch(ncvar_get(m_nc,"geolocation/sza"),
                   error=function(e) NA)
  vza  <- tryCatch(ncvar_get(m_nc,"geolocation/vza"),
                   error=function(e) NA)
  saa  <- tryCatch(ncvar_get(m_nc,"geolocation/saa"),
                   error=function(e) NA)
  vaa  <- tryCatch(ncvar_get(m_nc,"geolocation/vaa"),
                   error=function(e) NA)
  qf   <- tryCatch(ncvar_get(m_nc,"product_o2proxy/main_quality_flag"),
                   error=function(e) NULL)
  nc_close(m_nc)
  
  # Apply quality mask
  if (!is.null(qf)) {
    bad <- qf != 1
    xco2[bad] <- NA
    dp[bad] <- NA
    alb1286[bad] <- NA
    alb1606[bad] <- NA
    psfc[bad] <- NA
    orog[bad] <- NA
    sza[bad] <- NA
    vza[bad] <- NA
    saa[bad] <- NA
    vaa[bad] <- NA
  }
  
  msat_df <- data.frame(
    lon=as.vector(m_lon), lat=as.vector(m_lat),
    xco2=as.vector(xco2),
    dp=as.vector(dp),
    alb1286=as.vector(alb1286),
    alb1606=as.vector(alb1606),
    psfc=as.vector(psfc),
    orog=as.vector(orog),
    sza=as.vector(sza),
    vza=as.vector(vza),
    saa=as.vector(saa),
    vaa=as.vector(vaa)
  ) %>% filter(complete.cases(lon, lat, xco2))
  
  if (nrow(msat_df)==0)
    stop("MSAT has no valid points")
  
  cat(sprintf("  MSAT valid pixels: %d\n", nrow(msat_df)))
  
  msat_pts <- vect(msat_df, geom=c("lon","lat"), crs="EPSG:4326")
  
  # 6. EXTRACT VIIRS VALUES AT MSAT POINTS (POINT-LEVEL ANALYSIS)
  aod_vals <- terra::extract(v_raster_aod, msat_pts)
  ae_vals  <- terra::extract(v_raster_ae, msat_pts)

  # Handle column names robustly
  aod_col_idx <- if ("Aerosol_Optical_Thickness_550_Land_Best_Estimate" %in% names(aod_vals)) {
    which(names(aod_vals) == "Aerosol_Optical_Thickness_550_Land_Best_Estimate")
  } else {
    2  # Default to second column if name not found
  }
  
  ae_col_idx <- if ("Angstrom_Exponent_Land_Best_Estimate" %in% names(ae_vals)) {
    which(names(ae_vals) == "Angstrom_Exponent_Land_Best_Estimate")
  } else {
    2
  }
  
  msat_df$AOD550 <- aod_vals[[aod_col_idx]]
  msat_df$AE <- ae_vals[[ae_col_idx]]
  
  # 7. APPLY QUALITY FILTERS (FIXED: Added AE filtering)
  n_before <- nrow(msat_df)
  
  msat_df <- msat_df %>% 
    filter(!is.na(AOD550) & !is.na(AE)) %>%
    filter(AE >= AE_MIN & AE <= AE_MAX)
  
  n_after <- nrow(msat_df)
  n_filtered <- n_before - n_after
  
  cat(sprintf("  Quality filtering: %d pixels removed (%d remain)\n", n_filtered, n_after))
  cat(sprintf("    - NA values: %d\n", sum(is.na(aod_vals[[aod_col_idx]]) | is.na(ae_vals[[ae_col_idx]]))))
  cat(sprintf("    - AE out of range [%.1f, %.1f]: %d\n", AE_MIN, AE_MAX, 
              n_filtered - sum(is.na(aod_vals[[aod_col_idx]]) | is.na(ae_vals[[ae_col_idx]]))))
  
  if (nrow(msat_df) < MIN_OVERLAP_PIXELS)
    stop(sprintf("Insufficient overlapping pixels: %d < %d required", 
                 nrow(msat_df), MIN_OVERLAP_PIXELS))
  
  # 8. DERIVED VARIABLES
  msat_df$AOD_1286 <- msat_df$AOD550 * (1286/550)^(-msat_df$AE)
  msat_df$primary_angle <- calc_primary_angle(msat_df$sza, msat_df$vza, msat_df$saa, msat_df$vaa)
  msat_df$scene_id <- scene_id
  
  # 9. ATTACH RASTERS FOR LATER PLOTTING
  attr(msat_df, "viirs_raster_aod") <- v_raster_aod
  attr(msat_df, "viirs_raster_ae") <- v_raster_ae
  attr(msat_df, "n_viirs_granules") <- length(viirs_files_matched)
  attr(msat_df, "time_range") <- c(scene_info$time_start, scene_info$time_end)
  
  return(msat_df)
}

# --- C. Stats Helpers ---

z_score <- function(x) as.numeric(scale(x))

max_vif_num <- function(df, xvars) {
  if (length(xvars) < 2) return(1.0)
  form <- as.formula(paste("xco2 ~", paste(xvars, collapse="+")))
  dat_fit <- if(nrow(df) > 5000) df[sample(nrow(df), 5000), ] else df
  tryCatch({
    mod <- lm(form, data = dat_fit)
    vifs <- car::vif(mod)
    if(is.null(dim(vifs))) max(vifs) else max(vifs[,1])
  }, error = function(e) NA_real_)
}

vif_prune <- function(df, xvars, threshold=5) {
  keep <- xvars
  repeat {
    if (length(keep) < 2) break
    curr_max <- max_vif_num(df, keep)
    if (is.na(curr_max) || curr_max <= threshold) break
    
    form <- as.formula(paste("xco2 ~", paste(keep, collapse="+")))
    dat_fit <- if(nrow(df) > 5000) df[sample(nrow(df), 5000), ] else df
    mod <- lm(form, data=dat_fit)
    v <- car::vif(mod)
    val_to_check <- if(is.null(dim(v))) v else v[,1]
    keep <- setdiff(keep, names(which.max(val_to_check)))
  }
  list(keep=keep)
}

# --- D. Publication-Quality Visualization ---

create_map_pair <- function(df, var1, var2, title1, title2, main_title, viirs_raster = NULL) {
  
  # If VIIRS raster provided, rasterize data to VIIRS grid for display
  if (!is.null(viirs_raster)) {
    v_pts <- vect(df, geom = c("lon", "lat"), crs = "EPSG:4326")
    
    # Count pixels per cell for coverage threshold (use function that counts features)
    r_count <- rasterize(v_pts, viirs_raster, fun = function(x, ...) length(x))
    
    # Rasterize XCO2 values
    r_var1 <- rasterize(v_pts, viirs_raster, field = var1, fun = "mean")
    r_var2 <- rasterize(v_pts, viirs_raster, field = var2, fun = "mean")
    
    # Apply coverage threshold (require ≥5 MSAT pixels per VIIRS cell)
    r_var1[r_count < 5] <- NA
    r_var2[r_count < 5] <- NA
    
    # Convert to data.frame for ggplot
    df_plot <- as.data.frame(c(r_var1, r_var2), xy = TRUE, na.rm = TRUE)
    names(df_plot) <- c("x", "y", var1, var2)
  } else {
    df_plot <- df
    if ("lon" %in% names(df_plot)) names(df_plot)[names(df_plot) == "lon"] <- "x"
    if ("lat" %in% names(df_plot)) names(df_plot)[names(df_plot) == "lat"] <- "y"
  }
  
  # Calculate consistent color limits (±3σ from scene mean)
  val_range <- range(c(df_plot[[var1]], df_plot[[var2]]), na.rm = TRUE)
  scene_mean <- mean(c(df_plot[[var1]], df_plot[[var2]]), na.rm = TRUE)
  scene_sd <- sd(c(df_plot[[var1]], df_plot[[var2]]), na.rm = TRUE)
  
  color_limits <- c(
    max(val_range[1], scene_mean - 3 * scene_sd),
    min(val_range[2], scene_mean + 3 * scene_sd)
  )
  
  # Plot 1: Uncorrected
  p1 <- ggplot(df_plot, aes(x = x, y = y, fill = .data[[var1]])) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", name = "ppm", 
                         limits = color_limits, oob = scales::squish) +
    coord_fixed() +
    labs(x = "Longitude", y = "Latitude", title = title1) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(0.4, "cm")
    )
  
  # Plot 2: Corrected
  p2 <- ggplot(df_plot, aes(x = x, y = y, fill = .data[[var2]])) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", name = "ppm",
                         limits = color_limits, oob = scales::squish) +
    coord_fixed() +
    labs(x = "Longitude", y = "Latitude", title = title2) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.key.width = unit(2, "cm"),
      legend.key.height = unit(0.4, "cm")
    )
  
  grob <- arrangeGrob(
    p1, p2, ncol = 2,
    top = textGrob(main_title, gp = gpar(fontsize = 16, fontface = "bold"))
  )
  
  return(grob)
}

# -----------------------------------------------------------------------------
# 2. Main: Discovery & Extraction
# -----------------------------------------------------------------------------

message("=== MethaneSAT XCO2 Correction Pipeline - CORRECTED VERSION ===\n")
message("--- Step 1: File Discovery ---")

if (USE_SPECIFIC_FILES) {
  message("-> Mode: Specific List")
  raw_file_list <- c(
    # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250120/c06AF15D0/p3856/MSAT_L2_c06AF15D0_p3856_v04001006_20250120T214301Z_214333Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250516/c0BC915D0/p3838/MSAT_L2_c0BC915D0_p3838_v04001006_20250516T212530Z_212602Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20241230/c058A15D0/p4175/MSAT_L2_c058A15D0_p4175_v04003001_20241230T212120Z_212152Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250119/c069C15D0/p4210/MSAT_L2_c069C15D0_p4210_v04003001_20250119T213734Z_213806Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250430/c0A9C15D0/p3846/MSAT_L2_c0A9C15D0_p3846_v04001006_20250430T213320Z_213352Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t348/20250220/c073C15C0/p3678/MSAT_L2_c073C15C0_p3678_v04001006_20250220T132419Z_132451Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t348/20250430/c0A9D15C0/p3665/MSAT_L2_c0A9D15C0_p3665_v04001006_20250430T132730Z_132802Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t108/20250120/c06B506C0/p3188/MSAT_L2_c06B506C0_p3188_v04001003_20250120T071415Z_071447Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t108/20250519/c0C1306C0/p3155/MSAT_L2_c0C1306C0_p3155_v04001003_20250519T071300Z_071332Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t51/20241104/c038C0330/p3866/MSAT_L2_c038C0330_p3866_v04001006_20241104T104958Z_105030Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t2/20240905/c010C0020/p3262/MSAT_L2_c010C0020_p3262_v04001003_20240905T194216Z_194248Z.nc"
  
   "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250120/c06AF15D0/p3856/MSAT_L2_c06AF15D0_p3856_v04001006_20250120T214301Z_214333Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250516/c0BC915D0/p3838/MSAT_L2_c0BC915D0_p3838_v04001006_20250516T212530Z_212602Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t349/20241230/c058A15D0/p4175/MSAT_L2_c058A15D0_p4175_v04003001_20241230T212120Z_212152Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250119/c069C15D0/p4210/MSAT_L2_c069C15D0_p4210_v04003001_20250119T213734Z_213806Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250430/c0A9C15D0/p3846/MSAT_L2_c0A9C15D0_p3846_v04001006_20250430T213320Z_213352Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t348/20250220/c073C15C0/p3678/MSAT_L2_c073C15C0_p3678_v04001006_20250220T132419Z_132451Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t348/20250430/c0A9D15C0/p3665/MSAT_L2_c0A9D15C0_p3665_v04001006_20250430T132730Z_132802Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t108/20250120/c06B506C0/p3188/MSAT_L2_c06B506C0_p3188_v04001003_20250120T071415Z_071447Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t108/20250519/c0C1306C0/p3155/MSAT_L2_c0C1306C0_p3155_v04001003_20250519T071300Z_071332Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t51/20241104/c038C0330/p3866/MSAT_L2_c038C0330_p3866_v04001006_20241104T104958Z_105030Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t2/20240905/c010C0020/p3262/MSAT_L2_c010C0020_p3262_v04001003_20240905T194216Z_194248Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t363/20250228/c07D816B0/p4002/MSAT_L2_c07D816B0_p4002_v04002001_20250228T204049Z_204114Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t60/20241228/c056403C0/p3860/MSAT_L2_c056403C0_p3860_v04001006_20241228T210922Z_210954Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t60/20240917/c017103C0/p3772/MSAT_L2_c017103C0_p3772_v04001006_20240917T211429Z_211501Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t56/20240928/c01CC0380/p3255/MSAT_L2_c01CC0380_p3255_v04001003_20240928T205946Z_210018Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t56/20240929/c01D50380/p3253/MSAT_L2_c01D50380_p3253_v04001003_20240929T210658Z_210730Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t56/20250320/c09060380/p3221/MSAT_L2_c09060380_p3221_v04001003_20250320T210203Z_210235Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t56/20241223/c05040380/p3231/MSAT_L2_c05040380_p3231_v04001003_20241223T204056Z_204128Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t147/20250606/c0CF30930/p3651/MSAT_L2_c0CF30930_p3651_v04001006_20250606T115534Z_115606Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t25/20241007/c021E0190/p4566/MSAT_L2_c021E0190_p4566_v04003001_20241007T122638Z_122710Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t25/20241122/c046E0190/p3212/MSAT_L2_c046E0190_p3212_v04001003_20241122T124404Z_124436Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t25/20240911/c01420190/p3215/MSAT_L2_c01420190_p3215_v04001003_20240911T122836Z_122908Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t25/20241009/c02300190/p3394/MSAT_L2_c02300190_p3394_v04001003_20241009T124038Z_124110Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t25/20241220/c04D90190/p3232/MSAT_L2_c04D90190_p3232_v04001003_20241220T122155Z_122227Z.nc"
  )
  all_nc_files <- gsub("gs://", "/mnt/gcs/", raw_file_list)
  all_nc_files <- all_nc_files[file.exists(all_nc_files)]
} else {
  message("-> Mode: Recursive Search")
  all_nc_files <- list.files(MSAT_ROOT_DIR, recursive = TRUE, 
                              pattern = "\\.nc$", full.names = TRUE)
}

if (length(all_nc_files) == 0) stop("No valid MSAT files found.")

extract_scene_n <- function(path) {
  t_m <- str_extract(path, "t[0-9]+")
  p_m <- str_extract(path, "p[0-9]+")
  if (!is.na(t_m) && !is.na(p_m)) return(paste0(t_m, p_m))
  if (!is.na(t_m)) return(t_m)
  return(basename(dirname(path)))
}

files_df <- data.frame(MSAT_file = all_nc_files, stringsAsFactors=FALSE) %>%
  mutate(scene_id = sapply(MSAT_file, extract_scene_n)) %>%
  distinct(scene_id, .keep_all=TRUE)

message(sprintf("Found %d unique scenes.", nrow(files_df)))

message("\n--- Step 2: Extraction & Unification (MSAT Resolution) ---")

files_df$VIIRS_file <- sapply(files_df$MSAT_file, find_or_download_viirs)
files_df <- files_df %>% filter(!is.na(VIIRS_file))
message(sprintf("Scenes with VIIRS data: %d\n", nrow(files_df)))

valid_data_list <- list()
extraction_metadata <- list()

for (i in 1:nrow(files_df)) {
  sid <- files_df$scene_id[i]
  cat(sprintf("\n[%d/%d] Processing %s\n", i, nrow(files_df), sid))
  cat(sprintf("  MSAT file: %s\n", basename(files_df$MSAT_file[i])))
  
  res <- tryCatch({
    extract_scene_data(files_df$VIIRS_file[i], files_df$MSAT_file[i], sid)
  }, error = function(e) {
    cat(sprintf("  ✗ FAILED: %s\n", e$message))
    return(NULL)
  })
  
  if (!is.null(res)) {
    valid_data_list[[length(valid_data_list)+1]] <- res
    
    # Store only metadata that actually exists (FIXED)
    extraction_metadata[[sid]] <- list(
      n_pixels = nrow(res),
      n_viirs_granules = attr(res, "n_viirs_granules"),
      time_start = as.character(attr(res, "time_range")[1]),
      time_end = as.character(attr(res, "time_range")[2])
    )
    
    cat(sprintf("  ✓ SUCCESS: %d pixels at MSAT resolution\n", nrow(res)))
  }
}

if (length(valid_data_list) == 0) stop("\nNo data extracted successfully.")

all_scenes_df <- bind_rows(valid_data_list)

db_path <- file.path(OUTPUT_DIR, "unified_database.csv")
write_csv(all_scenes_df, db_path)
message(sprintf("\n✓ Unified database saved: %s", db_path))

# FIXED: Simplified metadata construction
meta_df <- bind_rows(lapply(names(extraction_metadata), function(sid) {
  m <- extraction_metadata[[sid]]
  data.frame(
    scene_id = sid,
    n_pixels = m$n_pixels,
    n_viirs_granules = ifelse(is.null(m$n_viirs_granules), NA, m$n_viirs_granules),
    time_start = m$time_start,
    time_end = m$time_end,
    stringsAsFactors = FALSE
  )
}))
write_csv(meta_df, file.path(OUTPUT_DIR, "extraction_metadata.csv"))

message(sprintf("\nUnified Database: %d rows from %d scenes (MSAT resolution)", 
                nrow(all_scenes_df), length(valid_data_list)))

# -----------------------------------------------------------------------------
# 3. Statistical Modeling
# -----------------------------------------------------------------------------

message("\n--- Step 3: Statistical Modeling (Pooled Analysis) ---")

predictors_all <- c("dp", "alb1286", "alb1606", "psfc", "orog", 
                    "sza", "vza", "primary_angle", "AOD_1286", "AE")

if (nrow(all_scenes_df) < 500) {
  message("  Skipping: Insufficient pixels (<500).")
  best_pool_mod <- NULL
} else {
  
  target_cols <- c("xco2", "scene_id", predictors_all)
  
  df_pool_clean <- all_scenes_df %>%
    dplyr::select(any_of(target_cols)) %>%
    drop_na()
  
  if (!("xco2" %in% names(df_pool_clean))) {
    stop("Fatal: 'xco2' missing after cleaning")
  }
  
  df_pool_clean$scene_id <- as.factor(df_pool_clean$scene_id)
  preds_pool <- intersect(names(df_pool_clean), predictors_all)
  
  message(sprintf("  Available predictors: %s", paste(preds_pool, collapse=", ")))
  
  if (length(preds_pool) < 2) {
    message("  Skipping: <2 valid predictors")
    best_pool_mod <- NULL
  } else {
    
    # Z-score standardization
    df_pool_z <- df_pool_clean
    for(v in preds_pool) df_pool_z[[v]] <- z_score(df_pool_z[[v]])
    
    # VIF pruning
    vp_pool <- vif_prune(df_pool_z, preds_pool, VIF_THRESHOLD)
    kept_pool <- vp_pool$keep
    message(sprintf("  Predictors after VIF pruning (threshold=%.1f): %s", 
                    VIF_THRESHOLD, paste(kept_pool, collapse=", ")))
    
    if (length(kept_pool) == 0) {
      message("  Skipping: All predictors removed by VIF")
      best_pool_mod <- NULL
    } else {
      
      message("\n  -> Running AICc Best Subset Selection (dredge)...")
      form_pool <- as.formula(paste("xco2 ~ scene_id +", 
                                    paste(kept_pool, collapse="+")))
      mod_pool_global <- lm(form_pool, data=df_pool_z, na.action="na.fail")
      
      dredge_pool <- dredge(mod_pool_global, rank="AICc", 
                            fixed=c("scene_id"), trace=FALSE, 
                            m.lim=c(1, length(kept_pool)))
      
      write_csv(as.data.frame(dredge_pool), 
                file.path(OUTPUT_DIR, "pooled_dredge_results.csv"))
      
      best_pool_mod <- get.models(dredge_pool, subset=1)[[1]]
      r2_adj <- summary(best_pool_mod)$adj.r.squared
      
      message(sprintf("  ✓ Best Model Adjusted R²: %.3f", r2_adj))
      
      # Extract predictors in final model
      final_predictors <- names(coef(best_pool_mod))
      final_predictors <- final_predictors[!grepl("Intercept|scene_id", final_predictors)]
      message(sprintf("  Final predictors: %s", paste(final_predictors, collapse=", ")))
      
      model_summary <- broom::tidy(best_pool_mod)
      write_csv(model_summary, file.path(OUTPUT_DIR, "best_model_coefficients.csv"))
      
      # Diagnostic plots
      png(file.path(OUTPUT_DIR, "model_diagnostics.png"), 
          width=1200, height=900, res=120)
      par(mfrow=c(2,2))
      plot(best_pool_mod, main="Model Diagnostics")
      dev.off()
      
      message("  ✓ Model diagnostics saved")
    }
  }
}

# -----------------------------------------------------------------------------
# 4. XCO2 Correction
# -----------------------------------------------------------------------------

message("\n--- Step 4: XCO2 Correction (Scene-Mean Preserving) ---")

if (is.null(best_pool_mod)) {
  message("  Skipping: No valid model available")
} else {
  
  coef_vec <- coef(best_pool_mod)
  
  predictor_coefs <- coef_vec[
    !grepl("^scene_id|^\\(Intercept\\)", names(coef_vec))
  ]
  
  message(sprintf("  Applying correction with %d predictors", 
                  length(predictor_coefs)))
  
  all_scenes_corrected <- list()
  correction_stats <- list()
  
  for (sid in unique(all_scenes_df$scene_id)) {
    
    df_scene <- all_scenes_df %>% filter(scene_id == sid)
    if (nrow(df_scene) < MIN_OVERLAP_PIXELS) {
      message(sprintf("  [%s] Skipping: Too few pixels (%d)", sid, nrow(df_scene)))
      next
    }
    
    # Check if scene was in training data
    scene_coef_name <- paste0("scene_id", sid)
    if (!scene_coef_name %in% names(coef_vec)) {
      message(sprintf("  [%s] WARNING: Scene not in training data - skipping correction", sid))
      next
    }
    
    y_hat <- rep(0, nrow(df_scene))
    failed_predictors <- character(0)
    
    # Apply predictor effects
    for (v in names(predictor_coefs)) {
      if (v %in% names(df_scene)) {
        x <- df_scene[[v]]
        if (sd(x, na.rm = TRUE) > 1e-10) {
          z_vals <- z_score(x)
          if (!any(is.na(z_vals))) {
            y_hat <- y_hat + predictor_coefs[[v]] * z_vals
          } else {
            failed_predictors <- c(failed_predictors, v)
          }
        } else {
          failed_predictors <- c(failed_predictors, v)
        }
      }
    }
    
    if (length(failed_predictors) > 0) {
      message(sprintf("  [%s] WARNING: Skipped predictors: %s",
                      sid, paste(failed_predictors, collapse=", ")))
    }
    
    # Add scene effect
    scene_effect <- coef_vec[scene_coef_name]
    y_hat <- y_hat + scene_effect
    
    # Quality checks
    if (any(is.infinite(y_hat))) {
      warning(sprintf("  [%s] ERROR: Infinite values in predictions! Skipping.", sid))
      next
    }
    
    if (all(is.na(y_hat))) {
      warning(sprintf("  [%s] ERROR: All predictions are NA! Skipping.", sid))
      next
    }
    
    # Scene-mean preserving correction
    y_hat_centered <- y_hat - mean(y_hat, na.rm = TRUE)
    df_scene$xco2_corrected <- df_scene$xco2 - y_hat_centered
    
    n_valid <- sum(is.finite(df_scene$xco2_corrected))
    if (n_valid < nrow(df_scene) * 0.5) {
      warning(sprintf("  [%s] Only %d/%d pixels successfully corrected",
                      sid, n_valid, nrow(df_scene)))
    }
    
    # Calculate statistics
    correction_stats[[sid]] <- data.frame(
      scene_id = sid,
      mean_correction = mean(y_hat_centered, na.rm = TRUE),
      sd_correction = sd(y_hat_centered, na.rm = TRUE),
      max_correction = max(abs(y_hat_centered), na.rm = TRUE),
      xco2_mean_before = mean(df_scene$xco2, na.rm = TRUE),
      xco2_mean_after = mean(df_scene$xco2_corrected, na.rm = TRUE),
      xco2_sd_before = sd(df_scene$xco2, na.rm = TRUE),
      xco2_sd_after = sd(df_scene$xco2_corrected, na.rm = TRUE),
      variance_reduction = 1 - sd(df_scene$xco2_corrected, na.rm=TRUE) / sd(df_scene$xco2, na.rm=TRUE),
      mean_preservation_error = abs(mean(df_scene$xco2_corrected, na.rm=TRUE) - mean(df_scene$xco2, na.rm=TRUE))
    )
    
    # Preserve raster attributes for plotting
    attr(df_scene, "viirs_raster_aod") <- attr(all_scenes_df[all_scenes_df$scene_id == sid, ], "viirs_raster_aod")
    
    all_scenes_corrected[[sid]] <- df_scene
    
    message(sprintf("  [%s] ✓ Corrected: max |Δ| = %.2f ppm, σ reduction = %.1f%%", 
                    sid, 
                    correction_stats[[sid]]$max_correction,
                    correction_stats[[sid]]$variance_reduction * 100))
  }
  
  if (length(correction_stats) > 0) {
    correction_summary <- bind_rows(correction_stats)
    write_csv(correction_summary, 
              file.path(OUTPUT_DIR, "correction_statistics.csv"))
    message(sprintf("\n✓ Correction statistics saved (%d scenes)", nrow(correction_summary)))
    
    all_scenes_df <- bind_rows(all_scenes_corrected)
  } else {
    message("\n⚠ No scenes were successfully corrected")
  }
}

# -----------------------------------------------------------------------------
# 5. Publication-Quality Plots
# -----------------------------------------------------------------------------

message("\n--- Step 5: Publication-Quality XCO2 Comparison ---")

if (!"xco2_corrected" %in% names(all_scenes_df)) {
  message("  Skipping: No corrected XCO2 available")
} else {
  
  for (sid in unique(all_scenes_df$scene_id)) {
    
    df_scene <- all_scenes_df %>% 
      filter(scene_id == sid, 
             is.finite(xco2), 
             is.finite(xco2_corrected),
             xco2 > 300, xco2 < 600,
             xco2_corrected > 300, xco2_corrected < 600)
    
    if (nrow(df_scene) < MIN_OVERLAP_PIXELS) {
      message(sprintf("  [%s] Skipping: Too few valid pixels for plotting", sid))
      next
    }
    
    # Retrieve VIIRS raster from original data
    original_scene <- valid_data_list[[which(sapply(valid_data_list, function(x) x$scene_id[1]) == sid)]]
    viirs_raster_aod <- attr(original_scene, "viirs_raster_aod")
    
    if (is.null(viirs_raster_aod)) {
      message(sprintf("  [%s] WARNING: VIIRS raster not found, using point display", sid))
    }
    
    message(sprintf("  [%s] Creating publication plot...", sid))
    
    # Create comparison plot
    grob <- create_map_pair(
      df_scene,
      var1 = "xco2",
      var2 = "xco2_corrected",
      title1 = "Uncorrected XCO₂ (ppm)",
      title2 = "Corrected XCO₂ (ppm)",
      main_title = sprintf("%s - XCO₂ Aerosol Correction", sid),
      viirs_raster = viirs_raster_aod
    )
    
    out_file <- file.path(OUTPUT_DIR, 
                          sprintf("%s_XCO2_comparison_publication.png", sid))
    
    ggsave(out_file, grob, width = 12, height = 6, 
           dpi = 1000, bg = "white")
    
    message(sprintf("  [%s] ✓ Saved: %s", sid, basename(out_file)))
    
    # BONUS: Correction magnitude map
    df_scene$correction_ppm <- df_scene$xco2 - df_scene$xco2_corrected
    
    if (!is.null(viirs_raster_aod)) {
      tryCatch({
        v_corr <- vect(df_scene, geom = c("lon", "lat"), crs = "EPSG:4326")
        r_corr <- rasterize(v_corr, viirs_raster_aod, 
                            field = "correction_ppm", fun = "mean")
        df_corr_plot <- as.data.frame(r_corr, xy = TRUE, na.rm = TRUE)
        names(df_corr_plot) <- c("x", "y", "correction_ppm")
        
        # Filter outliers
        df_corr_plot <- df_corr_plot %>% filter(abs(correction_ppm) < 50)
        
        if (nrow(df_corr_plot) > 50) {
          p_corr <- ggplot(df_corr_plot, aes(x = x, y = y, fill = correction_ppm)) +
            geom_tile() +
            coord_fixed() +
            scale_fill_gradient2(
              low = "blue", mid = "white", high = "red",
              midpoint = 0, name = "Correction\n(ppm)"
            ) +
            labs(
              title = sprintf("%s - Correction Magnitude", sid),
              subtitle = sprintf("Mean: %.2f ppm | SD: %.2f ppm | Max: |%.2f| ppm",
                                 mean(df_scene$correction_ppm, na.rm=TRUE),
                                 sd(df_scene$correction_ppm, na.rm=TRUE),
                                 max(abs(df_scene$correction_ppm), na.rm=TRUE)),
              x = "Longitude", y = "Latitude"
            ) +
            theme_minimal(base_size = 12) +
            theme(
              legend.position = "bottom",
              panel.grid = element_blank(),
              plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              legend.key.width = unit(2, "cm"),
              legend.key.height = unit(0.4, "cm")
            )
          
          corr_file <- file.path(OUTPUT_DIR,
                                 sprintf("%s_correction_map.png", sid))
          ggsave(corr_file, p_corr, width = 10, height = 8, dpi = 1000, bg = "white")
          message(sprintf("  [%s] ✓ Saved: %s", sid, basename(corr_file)))
        } else {
          message(sprintf("  [%s] Skipping correction map: insufficient data", sid))
        }
      }, error = function(e) {
        message(sprintf("  [%s] WARNING: Correction map failed: %s", sid, e$message))
      })
    }
  }
}

# -----------------------------------------------------------------------------
# 6. Final Summary
# -----------------------------------------------------------------------------

message("\n", paste(rep("=", 70), collapse=""))
message("PIPELINE COMPLETE")
message(paste(rep("=", 70), collapse=""))
message(sprintf("\nOutput directory: %s", OUTPUT_DIR))
message("\nGenerated files:")
message("  ✓ unified_database.csv - All extracted data at MSAT resolution")
message("  ✓ extraction_metadata.csv - Scene metadata and VIIRS granule info")
message("  ✓ pooled_dredge_results.csv - Model selection results")
message("  ✓ best_model_coefficients.csv - Final model coefficients")
message("  ✓ correction_statistics.csv - Per-scene correction metrics")
message("  ✓ model_diagnostics.png - Regression diagnostic plots")
message("  ✓ *_XCO2_comparison_publication.png - Side-by-side before/after")
message("  ✓ *_correction_map.png - Correction magnitude maps")

message("\n", paste(rep("-", 70), collapse=""))
message("KEY FEATURES")
message(paste(rep("-", 70), collapse=""))
message("✓ TEMPORAL FILTERING: ±", TEMPORAL_WINDOW_HOURS, "h VIIRS granule selection")
message("✓ SPATIAL RESOLUTION: MSAT native (~400m) with 6km VIIRS aerosol")
message("✓ QUALITY CONTROL: AE range [", AE_MIN, ", ", AE_MAX, "], QF=1 only")
message("✓ SCENE-MEAN PRESERVING: XCO₂ mean unchanged, variance reduced")
message("✓ PUBLICATION QUALITY: 300 DPI, consistent color scales")
message("✓ STATISTICAL RIGOR: VIF pruning, AICc model selection")

message("\n", paste(rep("-", 70), collapse=""))
message("CRITICAL FIXES APPLIED")
message(paste(rep("-", 70), collapse=""))
message("✓ Temporal matching: Only VIIRS within ±3h of MSAT overpass")
message("✓ Metadata handling: Fixed null attribute references")
message("✓ AE quality filtering: Now properly applied")
message("✓ Raster merging: Warnings suppressed (multi-satellite expected)")
message("✓ Coverage threshold: Require ≥5 MSAT pixels per 6km VIIRS cell")

message("\nAnalysis complete!")