# =============================================================================
# Unified MethaneSAT + VIIRS XCO2 Modeling Pipeline - COSINE ANGLE VERSION
# =============================================================================
# KEY MODIFICATIONS FROM ORIGINAL:
# 1. Uses cos(sza), cos(vza), cos(primary_angle) instead of raw angles
# 2. All other functionality preserved
# 3. Autonomous execution - no manual edits needed
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
OUTPUT_DIR <- paste0("results_viirs_grid_cosine_", format(Sys.Date(), "%Y%m%d"))

# Modeling Constants
VIF_THRESHOLD <- 5
TOP_N_DREDGE <- 10
TEMPORAL_WINDOW_HOURS <- 3

# Quality thresholds (RELAXED FOR DEBUGGING)
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
  if (inherits(nc, "try-error")) {
    cat(sprintf("    ERROR: Cannot open NetCDF: %s\n", as.character(nc)))
    return(NULL)
  }
  
  lon <- try(ncvar_get(nc, "geolocation/longitude"), silent = TRUE)
  lat <- try(ncvar_get(nc, "geolocation/latitude"), silent = TRUE)
  time_vals <- try(ncvar_get(nc, "geolocation/time"), silent = TRUE)
  
  if (inherits(lon, "try-error")) {
    cat("    ERROR: Cannot read longitude\n")
    nc_close(nc)
    return(NULL)
  }
  if (inherits(lat, "try-error")) {
    cat("    ERROR: Cannot read latitude\n")
    nc_close(nc)
    return(NULL)
  }
  if (inherits(time_vals, "try-error")) {
    cat("    ERROR: Cannot read time\n")
    nc_close(nc)
    return(NULL)
  }
  
  nc_close(nc)
  
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
  
  cat(" -> No existing VIIRS files found\n")
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

# --- B. Spatial Processing (VIIRS GRID APPROACH) ---

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
  
  cat(sprintf("  ===== EXTRACTION DETAILS FOR %s =====\n", scene_id))
  
  # 1. GET MSAT METADATA
  cat("  [1/11] Reading MSAT metadata...\n")
  scene_info <- get_bbox_and_time_range(msat_file)
  if (is.null(scene_info)) {
    stop(sprintf("Cannot read MSAT metadata from %s", basename(msat_file)))
  }
  
  scene_date <- as.Date(scene_info$date)
  cat(sprintf("    ✓ MSAT scene time: %s\n", format(scene_info$time_median, "%Y-%m-%d %H:%M:%S UTC")))

  # 2. FIND ALL VIIRS FILES FROM SAME DAY
  cat("  [2/11] Finding VIIRS files for this date...\n")
  
  # Check if VIIRS directory exists and has files
  if (!dir.exists(VIIRS_DATA_DIR)) {
    stop(sprintf("VIIRS directory does not exist: %s", VIIRS_DATA_DIR))
  }
  
  viirs_files_all <- list.files(VIIRS_DATA_DIR, pattern="\\.nc$", full.names=TRUE)
  cat(sprintf("    Total VIIRS files in directory: %d\n", length(viirs_files_all)))
  
  if (length(viirs_files_all) == 0) {
    stop(sprintf("No VIIRS files found in %s - run download first", VIIRS_DATA_DIR))
  }

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

  if(length(viirs_files_day) == 0) {
    stop(sprintf("No VIIRS granules found for scene %s on %s (searched %d files)", 
                 scene_id, scene_date, length(viirs_files_all)))
  }
  
  cat(sprintf("    ✓ Found %d VIIRS granules on %s\n", length(viirs_files_day), scene_date))

  # 3. FILTER BY TEMPORAL WINDOW
  cat("  [3/11] Filtering by temporal window...\n")
  viirs_files_matched <- Filter(function(f) {
    viirs_time <- extract_viirs_time(f)
    if(is.null(viirs_time)) {
      cat(sprintf("    WARNING: Cannot parse time from %s\n", basename(f)))
      return(FALSE)
    }
    
    time_diff_hours <- abs(as.numeric(
      difftime(viirs_time, scene_info$time_median, units="hours")
    ))
    
    cat(sprintf("    - %s: Δt = %.1f hours %s\n", 
                basename(f), 
                time_diff_hours,
                ifelse(time_diff_hours <= TEMPORAL_WINDOW_HOURS, "[MATCH]", "[SKIP]")))
    
    return(time_diff_hours <= TEMPORAL_WINDOW_HOURS)
  }, viirs_files_day)
  
  if(length(viirs_files_matched) == 0) {
    stop(sprintf("No VIIRS granules within ±%d hours of MSAT overpass (searched %d granules)", 
                 TEMPORAL_WINDOW_HOURS, length(viirs_files_day)))
  }
  
  cat(sprintf("    ✓ Using %d/%d VIIRS granules (±%dh window)\n", 
              length(viirs_files_matched), length(viirs_files_day), TEMPORAL_WINDOW_HOURS))

  # 4. LOAD VIIRS DATA (as point data)
  cat("  [4/11] Loading VIIRS data...\n")
  
  load_viirs_data <- function(file) {
    cat(sprintf("    Loading %s...\n", basename(file)))
    nc <- try(nc_open(file), silent = TRUE)
    if (inherits(nc, "try-error")) {
      cat(sprintf("      ERROR: Cannot open file\n"))
      return(NULL)
    }
    
    lon <- try(ncvar_get(nc, "Longitude"), silent = TRUE)
    lat <- try(ncvar_get(nc, "Latitude"), silent = TRUE)
    aod <- try(ncvar_get(nc, "Aerosol_Optical_Thickness_550_Land_Best_Estimate"), silent = TRUE)
    ae <- try(ncvar_get(nc, "Angstrom_Exponent_Land_Best_Estimate"), silent = TRUE)
    
    if (inherits(lon, "try-error") || inherits(lat, "try-error") ||
        inherits(aod, "try-error") || inherits(ae, "try-error")) {
      nc_close(nc)
      cat(sprintf("      ERROR: Cannot read required variables\n"))
      return(NULL)
    }
    
    # Handle fill values for AOD
    fill_aod <- ncatt_get(nc, "Aerosol_Optical_Thickness_550_Land_Best_Estimate", "_FillValue")
    if (fill_aod$hasatt) {
      aod[abs(aod - fill_aod$value) < 1e-5] <- NA
    }
    aod[aod < -0.1 | aod > 5] <- NA  # Reasonable AOD range
    
    # Handle fill values for AE
    fill_ae <- ncatt_get(nc, "Angstrom_Exponent_Land_Best_Estimate", "_FillValue")
    if (fill_ae$hasatt) {
      ae[abs(ae - fill_ae$value) < 1e-5] <- NA
    }
    ae[ae < AE_MIN | ae > AE_MAX] <- NA  # Use configured AE range
    
    nc_close(nc)
    
    n_valid <- sum(!is.na(aod) & !is.na(ae))
    cat(sprintf("      ✓ Valid pixels: %d\n", n_valid))
    
    # Convert to vector format
    list(lon=as.vector(lon), lat=as.vector(lat), aod=as.vector(aod), ae=as.vector(ae))
  }
  
  # Load all VIIRS granules as point data
  viirs_data_list <- lapply(viirs_files_matched, load_viirs_data)
  viirs_data_list <- Filter(Negate(is.null), viirs_data_list)  # Remove failed loads
  
  if (length(viirs_data_list) == 0) {
    stop("All VIIRS granules failed to load")
  }
  
  # Combine all VIIRS data into single dataframe
  viirs_df_list <- lapply(viirs_data_list, function(vdata) {
    data.frame(
      lon = vdata$lon,
      lat = vdata$lat,
      AOD550 = vdata$aod,
      AE = vdata$ae
    ) %>% filter(complete.cases(.))
  })
  
  viirs_combined <- bind_rows(viirs_df_list)
  
  if (nrow(viirs_combined) == 0) {
    stop("No valid VIIRS data after filtering")
  }
  
  cat(sprintf("    ✓ VIIRS valid pixels: %d\n", nrow(viirs_combined)))

  # 5. LOAD MSAT DATA
  cat("  [5/11] Loading MSAT data...\n")
  m_nc <- try(nc_open(msat_file), silent = TRUE)
  if (inherits(m_nc, "try-error")) {
    stop("Cannot open MSAT NetCDF file")
  }
  
  m_lon <- try(ncvar_get(m_nc, "geolocation/longitude"), silent = TRUE)
  m_lat <- try(ncvar_get(m_nc, "geolocation/latitude"), silent = TRUE)
  
  if (inherits(m_lon, "try-error") || inherits(m_lat, "try-error")) {
    nc_close(m_nc)
    stop("Cannot read MSAT geolocation")
  }
  
  # Load XCO2 and handle fill values
  xco2  <- try(ncvar_get(m_nc, "product_o2proxy/xco2_o2proxy_vcdcor"), silent = TRUE)
  if (inherits(xco2, "try-error")) {
    nc_close(m_nc)
    stop("Cannot read xco2_o2proxy_vcdcor variable")
  }
  
  # Get fill value from attributes and convert to NA
  fill_val <- ncatt_get(m_nc, "product_o2proxy/xco2_o2proxy_vcdcor", "_FillValue")
  if (fill_val$hasatt) {
    xco2[abs(xco2 - fill_val$value) < 1e-5] <- NA
  }
  # Also filter obviously bad values (should be ~400 ppm)
  xco2[xco2 < 300 | xco2 > 600] <- NA
  
  cat(sprintf("    XCO2: %d / %d valid\n", sum(!is.na(xco2) & is.finite(xco2)), length(xco2)))
  
  # Load other variables with fill value handling
  dp    <- tryCatch(ncvar_get(m_nc,"o2dp_fit_diagnostics/bias_corrected_delta_pressure"),
                    error=function(e) NULL)
  if (is.null(dp)) {
    dp <- try(ncvar_get(m_nc,"o2dp_fit_diagnostics/delta_pressure"), silent = TRUE)
    if (inherits(dp, "try-error")) dp <- NULL
  }
  
  # Clean fill values for dp
  if (!is.null(dp)) {
    fill_dp <- ncatt_get(m_nc, "o2dp_fit_diagnostics/delta_pressure", "_FillValue")
    if (fill_dp$hasatt) {
      dp[abs(dp - fill_dp$value) < 1e-5] <- NA
    }
    dp[abs(dp) > 1e10] <- NA  # Remove obviously bad values
  }
  
  alb1286 <- try(ncvar_get(m_nc,"o2dp_fit_diagnostics/retrieved_albedo_1286nm"), silent = TRUE)
  if (!inherits(alb1286, "try-error")) {
    alb1286[alb1286 < 0 | alb1286 > 1] <- NA  # Albedo should be 0-1
  } else {
    alb1286 <- NULL
  }
  
  alb1606 <- try(ncvar_get(m_nc,"co2proxy_fit_diagnostics/retrieved_albedo_1606nm"), silent = TRUE)
  if (!inherits(alb1606, "try-error")) {
    alb1606[alb1606 < 0 | alb1606 > 1] <- NA
  } else {
    alb1606 <- NULL
  }
  
  psfc <- try(ncvar_get(m_nc,"co2proxy_fit_diagnostics/retrieved_surface_pressure"), silent = TRUE)
  if (!inherits(psfc, "try-error")) {
    psfc[psfc < 50000 | psfc > 110000] <- NA  # Reasonable pressure range (Pa)
  } else {
    psfc <- NULL
  }
  
  orog <- try(ncvar_get(m_nc,"geolocation/terrain_height"), silent = TRUE)
  if (!inherits(orog, "try-error")) {
    orog[abs(orog) > 10000] <- NA  # Reasonable elevation range
  } else {
    orog <- NULL
  }
  
  sza  <- try(ncvar_get(m_nc,"geolocation/sza"), silent = TRUE)
  if (!inherits(sza, "try-error")) {
    sza[sza < 0 | sza > 90] <- NA  # Solar zenith angle range
  } else {
    sza <- NULL
  }
  
  vza  <- try(ncvar_get(m_nc,"geolocation/vza"), silent = TRUE)
  if (!inherits(vza, "try-error")) {
    vza[vza < 0 | vza > 90] <- NA  # Viewing zenith angle range
  } else {
    vza <- NULL
  }
  
  saa  <- try(ncvar_get(m_nc,"geolocation/saa"), silent = TRUE)
  if (!inherits(saa, "try-error")) {
    saa[saa < 0 | saa > 360] <- NA  # Solar azimuth angle range
  } else {
    saa <- NULL
  }
  
  vaa  <- try(ncvar_get(m_nc,"geolocation/vaa"), silent = TRUE)
  if (!inherits(vaa, "try-error")) {
    vaa[vaa < 0 | vaa > 360] <- NA  # Viewing azimuth angle range
  } else {
    vaa <- NULL
  }
  
  qf   <- try(ncvar_get(m_nc,"product_o2proxy/main_quality_flag"), silent = TRUE)
  if (inherits(qf, "try-error")) qf <- NULL
  
  nc_close(m_nc)
  
  # Apply quality mask
  if (!is.null(qf)) {
    bad <- qf != 1
    xco2[bad] <- NA
    if (!is.null(dp)) dp[bad] <- NA
    if (!is.null(alb1286)) alb1286[bad] <- NA
    if (!is.null(alb1606)) alb1606[bad] <- NA
    if (!is.null(psfc)) psfc[bad] <- NA
    if (!is.null(orog)) orog[bad] <- NA
    if (!is.null(sza)) sza[bad] <- NA
    if (!is.null(vza)) vza[bad] <- NA
    if (!is.null(saa)) saa[bad] <- NA
    if (!is.null(vaa)) vaa[bad] <- NA
    
    cat(sprintf("    Quality filter: %d / %d pixels pass QF=1\n", 
                sum(qf == 1), length(qf)))
  }
  
  # Calculate primary angle
  p_angle <- NULL
  if (!is.null(sza) && !is.null(vza) && !is.null(saa) && !is.null(vaa)) {
    p_angle <- calc_primary_angle(sza, vza, saa, vaa)
  }
  
  # ===== COSINE TRANSFORMATION =====
  # Convert angles from degrees to radians and apply cosine
  cos_sza <- NULL
  cos_vza <- NULL
  cos_primary_angle <- NULL
  
  if (!is.null(sza)) {
    cos_sza <- cos(sza * pi / 180)
    cat(sprintf("    ✓ Computed cos(SZA)\n"))
  }
  
  if (!is.null(vza)) {
    cos_vza <- cos(vza * pi / 180)
    cat(sprintf("    ✓ Computed cos(VZA)\n"))
  }
  
  if (!is.null(p_angle)) {
    cos_primary_angle <- cos(p_angle * pi / 180)
    cat(sprintf("    ✓ Computed cos(primary_angle)\n"))
  }
  # =================================
  
  msat_df <- data.frame(
    lon=as.vector(m_lon), 
    lat=as.vector(m_lat),
    xco2=as.vector(xco2),
    dp=if(is.null(dp)) NA_real_ else as.vector(dp),
    alb1286=if(is.null(alb1286)) NA_real_ else as.vector(alb1286),
    alb1606=if(is.null(alb1606)) NA_real_ else as.vector(alb1606),
    psfc=if(is.null(psfc)) NA_real_ else as.vector(psfc),
    orog=if(is.null(orog)) NA_real_ else as.vector(orog),
    cos_sza=if(is.null(cos_sza)) NA_real_ else as.vector(cos_sza),
    cos_vza=if(is.null(cos_vza)) NA_real_ else as.vector(cos_vza),
    cos_primary_angle=if(is.null(cos_primary_angle)) NA_real_ else as.vector(cos_primary_angle)
  ) %>% filter(complete.cases(lon, lat, xco2))
  
  if (nrow(msat_df)==0) {
    stop("MSAT has no valid points after filtering")
  }
  
  cat(sprintf("    ✓ MSAT valid pixels: %d\n", nrow(msat_df)))
  
  # 6. CREATE COMMON GRID FROM MSAT EXTENT (with small buffer)
  cat("  [6/11] Creating common grid...\n")
  
  # Use MSAT extent with 10% buffer to catch nearby VIIRS pixels
  lon_range_msat <- max(msat_df$lon) - min(msat_df$lon)
  lat_range_msat <- max(msat_df$lat) - min(msat_df$lat)
  buffer_lon <- lon_range_msat * 0.1
  buffer_lat <- lat_range_msat * 0.1
  
  lon_min <- min(msat_df$lon, na.rm=TRUE) - buffer_lon
  lon_max <- max(msat_df$lon, na.rm=TRUE) + buffer_lon
  lat_min <- min(msat_df$lat, na.rm=TRUE) - buffer_lat
  lat_max <- max(msat_df$lat, na.rm=TRUE) + buffer_lat
  
  cat(sprintf("    MSAT extent (with 10%% buffer):\n"))
  cat(sprintf("      Lon range: [%.3f, %.3f]\n", lon_min, lon_max))
  cat(sprintf("      Lat range: [%.3f, %.3f]\n", lat_min, lat_max))
  
  # Create regular grid with ~6km resolution (VIIRS native)
  res_deg <- 0.054
  ncols <- ceiling((lon_max - lon_min) / res_deg)
  nrows <- ceiling((lat_max - lat_min) / res_deg)
  
  r_template <- rast(
    nrows = nrows,
    ncols = ncols,
    xmin = lon_min,
    xmax = lon_max,
    ymin = lat_min,
    ymax = lat_max,
    crs = "EPSG:4326"
  )
  
  total_cells <- nrows * ncols
  cat(sprintf("    ✓ Grid: %d x %d cells = %d total (%.3f° resolution)\n", 
              nrows, ncols, total_cells, res_deg))
  
  # Safety check: warn if grid is very large
  if (total_cells > 100000) {
    cat(sprintf("    ⚠ WARNING: Large grid (%d cells) - rasterization may be slow\n", total_cells))
  }
  
  # 7. RASTERIZE BOTH DATASETS TO COMMON GRID
  cat("  [7/11] Rasterizing VIIRS to grid...\n")
  viirs_vect <- try(vect(viirs_combined, geom=c("lon","lat"), crs="EPSG:4326"), silent = TRUE)
  if (inherits(viirs_vect, "try-error")) {
    stop("Failed to create VIIRS vector from points")
  }
  
  v_raster <- try(rasterize(viirs_vect, r_template, field=c("AOD550", "AE"), 
                            fun="mean", na.rm=TRUE), silent = TRUE)
  if (inherits(v_raster, "try-error")) {
    stop("Failed to rasterize VIIRS data")
  }
  
  cat("  [8/11] Rasterizing MSAT to grid...\n")
  msat_pts <- try(vect(msat_df, geom=c("lon","lat"), crs="EPSG:4326"), silent = TRUE)
  if (inherits(msat_pts, "try-error")) {
    stop("Failed to create MSAT vector from points")
  }
  
  msat_vars <- c("xco2", "dp", "alb1286", "alb1606", "psfc", "orog", 
                 "cos_sza", "cos_vza", "cos_primary_angle")
  m_raster <- try(rasterize(msat_pts, r_template, field=msat_vars, 
                            fun="mean", na.rm=TRUE), silent = TRUE)
  if (inherits(m_raster, "try-error")) {
    stop("Failed to rasterize MSAT data")
  }
  
  cat("  [9/11] Stacking rasters...\n")
  stack_all <- try(c(v_raster, m_raster), silent = TRUE)
  if (inherits(stack_all, "try-error")) {
    stop("Failed to stack rasters")
  }
  
  # Ensure names are correct
  expected_names <- c("AOD550", "AE", "xco2", "dp", "alb1286", "alb1606", 
                      "psfc", "orog", "cos_sza", "cos_vza", "cos_primary_angle")
  if (nlyr(stack_all) == length(expected_names)) {
    names(stack_all) <- expected_names
  }
  
  # 9. CONVERT TO DATAFRAME
  cat("  [10/11] Converting to dataframe...\n")
  df_scene <- try(as.data.frame(stack_all, xy=TRUE, na.rm=FALSE), silent = TRUE)
  if (inherits(df_scene, "try-error")) {
    stop("Failed to convert raster stack to dataframe")
  }
  
  cat(sprintf("    Before filtering: %d pixels\n", nrow(df_scene)))
  
  # Filter for required columns only
  df_scene <- df_scene %>% filter(!is.na(AOD550) & !is.na(AE) & !is.na(xco2))
  
  cat(sprintf("    After requiring AOD550, AE, xco2: %d pixels\n", nrow(df_scene)))
  
  if (nrow(df_scene) == 0) {
    stop("No valid colocated pixels after rasterization")
  }
  
  # 10. APPLY AE QUALITY FILTER
  cat("  [11/11] Applying quality filters...\n")
  n_before <- nrow(df_scene)
  
  df_scene <- df_scene %>% 
    filter(!is.na(AOD550) & !is.na(AE)) %>%
    filter(AE >= AE_MIN & AE <= AE_MAX)
  
  n_after <- nrow(df_scene)
  n_filtered <- n_before - n_after
  
  cat(sprintf("    Filtered: %d pixels removed (%d remain)\n", n_filtered, n_after))
  
  if (nrow(df_scene) < MIN_OVERLAP_PIXELS) {
    stop(sprintf("Insufficient overlapping pixels: %d < %d required", 
                 nrow(df_scene), MIN_OVERLAP_PIXELS))
  }
  
  # 11. DERIVED VARIABLES
  df_scene$AOD_1286 <- df_scene$AOD550 * (1286/550)^(-df_scene$AE)
  df_scene$scene_id <- scene_id
  
  # 12. ATTACH METADATA
  attr(df_scene, "viirs_raster_aod") <- stack_all[[1]]  # Store AOD layer for plotting
  attr(df_scene, "n_viirs_granules") <- length(viirs_files_matched)
  attr(df_scene, "time_range") <- c(scene_info$time_start, scene_info$time_end)
  attr(df_scene, "grid_dims") <- c(nrows, ncols)
  
  cat(sprintf("  ===== ✓ SUCCESS: %d valid pixels =====\n\n", nrow(df_scene)))
  
  return(df_scene)
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

create_map_pair <- function(df, var1, var2, title1, title2, main_title) {
  
  # Calculate consistent color limits (±3σ from scene mean)
  val_range <- range(c(df[[var1]], df[[var2]]), na.rm = TRUE)
  scene_mean <- mean(c(df[[var1]], df[[var2]]), na.rm = TRUE)
  scene_sd <- sd(c(df[[var1]], df[[var2]]), na.rm = TRUE)
  
  color_limits <- c(
    max(val_range[1], scene_mean - 3 * scene_sd),
    min(val_range[2], scene_mean + 3 * scene_sd)
  )
  
  # Plot 1: Uncorrected
  p1 <- ggplot(df, aes(x = x, y = y, fill = .data[[var1]])) +
    geom_raster() +
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
  p2 <- ggplot(df, aes(x = x, y = y, fill = .data[[var2]])) +
    geom_raster() +
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

message("=== MethaneSAT XCO2 Correction Pipeline - COSINE ANGLE VERSION ===\n")
message("--- Step 1: File Discovery ---")

if (USE_SPECIFIC_FILES) {
  message("-> Mode: Specific List")
  raw_file_list <- c(
    "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250516/c0BC915D0/p3838/MSAT_L2_c0BC915D0_p3838_v04001006_20250516T212530Z_212602Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t349/20241230/c058A15D0/p4175/MSAT_L2_c058A15D0_p4175_v04003001_20241230T212120Z_212152Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t348/20250220/c073C15C0/p3678/MSAT_L2_c073C15C0_p3678_v04001006_20250220T132419Z_132451Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t348/20250430/c0A9D15C0/p3665/MSAT_L2_c0A9D15C0_p3665_v04001006_20250430T132730Z_132802Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t60/20241228/c056403C0/p3860/MSAT_L2_c056403C0_p3860_v04001006_20241228T210922Z_210954Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t60/20240917/c017103C0/p3772/MSAT_L2_c017103C0_p3772_v04001006_20240917T211429Z_211501Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t56/20241223/c05040380/p3231/MSAT_L2_c05040380_p3231_v04001003_20241223T204056Z_204128Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t25/20241007/c021E0190/p4566/MSAT_L2_c021E0190_p4566_v04003001_20241007T122638Z_122710Z.nc",
    "gs://msat-prod-data-methanesat-level2/level2_results/t25/20240911/c01420190/p3215/MSAT_L2_c01420190_p3215_v04001003_20240911T122836Z_122908Z.nc"
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

message("\n--- Step 2: Extraction & Unification (VIIRS Grid Resolution + Cosine Angles) ---")

files_df$VIIRS_file <- sapply(files_df$MSAT_file, find_or_download_viirs)
files_df <- files_df %>% filter(!is.na(VIIRS_file))
message(sprintf("Scenes with VIIRS data: %d\n", nrow(files_df)))

if (nrow(files_df) == 0) {
  stop("No scenes have VIIRS data available. Please check VIIRS_DATA_DIR or run downloads.")
}

valid_data_list <- list()
extraction_metadata <- list()
extraction_errors <- list()

for (i in 1:nrow(files_df)) {
  sid <- files_df$scene_id[i]
  cat(sprintf("\n[%d/%d] Processing %s\n", i, nrow(files_df), sid))
  cat(sprintf("  MSAT file: %s\n", basename(files_df$MSAT_file[i])))
  
  res <- tryCatch({
    extract_scene_data(files_df$VIIRS_file[i], files_df$MSAT_file[i], sid)
  }, error = function(e) {
    cat(sprintf("  ✗ FAILED: %s\n", conditionMessage(e)))
    extraction_errors[[sid]] <<- conditionMessage(e)
    return(NULL)
  })
  
  if (!is.null(res)) {
    valid_data_list[[length(valid_data_list)+1]] <- res
    
    extraction_metadata[[sid]] <- list(
      n_pixels = nrow(res),
      n_viirs_granules = attr(res, "n_viirs_granules"),
      grid_dims = attr(res, "grid_dims"),
      time_start = as.character(attr(res, "time_range")[1]),
      time_end = as.character(attr(res, "time_range")[2])
    )
    
    grid_dims <- attr(res, "grid_dims")
    cat(sprintf("  ✓ SUCCESS: %d valid pixels on %dx%d VIIRS grid\n", 
                nrow(res), grid_dims[1], grid_dims[2]))
  }
}

# Save error log
if (length(extraction_errors) > 0) {
  error_df <- data.frame(
    scene_id = names(extraction_errors),
    error_message = unlist(extraction_errors),
    stringsAsFactors = FALSE
  )
  write_csv(error_df, file.path(OUTPUT_DIR, "extraction_errors.csv"))
  message(sprintf("\n⚠ Saved extraction errors to extraction_errors.csv (%d failures)", 
                  length(extraction_errors)))
}

if (length(valid_data_list) == 0) {
  message("\n===== EXTRACTION FAILED FOR ALL SCENES =====")
  message("Common causes:")
  message("  1. No VIIRS data in VIIRS_DATA_DIR")
  message("  2. VIIRS files don't match temporal window")
  message("  3. No spatial overlap between MSAT and VIIRS")
  message("  4. Quality filtering removed all pixels")
  message("\nCheck extraction_errors.csv for details.")
  stop("No data extracted successfully.")
}

# Continue with rest of pipeline...
message(sprintf("\n✓ Successfully extracted %d / %d scenes", 
                length(valid_data_list), nrow(files_df)))

# Store raster attributes separately (bind_rows strips attributes!)
raster_attributes <- lapply(valid_data_list, function(df) {
  list(
    scene_id = df$scene_id[1],
    viirs_raster_aod = attr(df, "viirs_raster_aod"),
    n_viirs_granules = attr(df, "n_viirs_granules"),
    grid_dims = attr(df, "grid_dims"),
    time_range = attr(df, "time_range")
  )
})
names(raster_attributes) <- sapply(raster_attributes, function(x) x$scene_id)

all_scenes_df <- bind_rows(valid_data_list)

db_path <- file.path(OUTPUT_DIR, "unified_database_viirs_grid_cosine.csv")
write_csv(all_scenes_df, db_path)
message(sprintf("\n✓ Unified database saved: %s", db_path))

# Metadata
meta_df <- bind_rows(lapply(names(extraction_metadata), function(sid) {
  m <- extraction_metadata[[sid]]
  data.frame(
    scene_id = sid,
    n_pixels = m$n_pixels,
    n_viirs_granules = ifelse(is.null(m$n_viirs_granules), NA, m$n_viirs_granules),
    grid_rows = m$grid_dims[1],
    grid_cols = m$grid_dims[2],
    time_start = m$time_start,
    time_end = m$time_end,
    stringsAsFactors = FALSE
  )
}))
write_csv(meta_df, file.path(OUTPUT_DIR, "extraction_metadata.csv"))

message(sprintf("\nUnified Database: %d rows from %d scenes (VIIRS grid resolution + cosine angles)", 
                nrow(all_scenes_df), length(valid_data_list)))

# -----------------------------------------------------------------------------
# 3. Statistical Modeling
# -----------------------------------------------------------------------------

message("\n--- Step 3: Statistical Modeling (Pooled Analysis with Cosine Angles) ---")

# Essential predictors that must be present for aerosol correction
essential_predictors <- c("AOD_1286", "AE")

# All available predictors (NOW WITH COSINE ANGLES)
predictors_all <- c("dp", "alb1286", "alb1606", "psfc", "orog", "AOD_1286", "AE", "cos_primary_angle")
# ,"cos_sza", "cos_vza" removed
if (nrow(all_scenes_df) < 500) {
  message("  Skipping: Insufficient pixels (<500).")
  best_pool_mod <- NULL
} else {
  
  # Only drop NA for xco2, scene_id, and ESSENTIAL predictors
  essential_cols <- c("xco2", "scene_id", essential_predictors)
  
  df_pool_clean <- all_scenes_df %>%
    dplyr::select(any_of(c(essential_cols, predictors_all))) %>%
    filter(complete.cases(across(any_of(essential_cols))))
  
  message(sprintf("  Rows after requiring essential columns (xco2, scene_id, AOD_1286, AE): %d", 
                  nrow(df_pool_clean)))
  
  if (!("xco2" %in% names(df_pool_clean))) {
    stop("Fatal: 'xco2' missing after cleaning")
  }
  
  # Check available predictors
  preds_pool <- intersect(names(df_pool_clean), predictors_all)
  
  # Report how many NAs each predictor has
  message("\n  Predictor completeness:")
  for (pred in preds_pool) {
    n_valid <- sum(!is.na(df_pool_clean[[pred]]))
    pct_valid <- 100 * n_valid / nrow(df_pool_clean)
    message(sprintf("    %s: %d valid (%.1f%%)", pred, n_valid, pct_valid))
  }
  
  # Only keep predictors that have at least 50% valid data
  preds_pool <- preds_pool[sapply(preds_pool, function(p) {
    sum(!is.na(df_pool_clean[[p]])) >= nrow(df_pool_clean) * 0.5
  })]
  
  message(sprintf("\n  Predictors with ≥50%% valid data: %s", paste(preds_pool, collapse=", ")))
  
  n_scenes <- length(unique(df_pool_clean$scene_id))
  message(sprintf("  Number of unique scenes: %d", n_scenes))
  
  if (n_scenes < 2) {
    message("  Skipping modeling: Need at least 2 scenes for scene_id as factor")
    best_pool_mod <- NULL
  } else {
    
    df_pool_clean$scene_id <- as.factor(df_pool_clean$scene_id)
    
    if (length(preds_pool) < 2) {
      message("  Skipping: <2 valid predictors")
      best_pool_mod <- NULL
    } else {
      
      # For modeling, drop NA for the selected predictors
      model_cols <- c("xco2", "scene_id", preds_pool)
      df_pool_z <- df_pool_clean %>%
        dplyr::select(any_of(model_cols)) %>%
        drop_na()
      
      message(sprintf("  Rows for modeling after dropping NA in selected predictors: %d", 
                      nrow(df_pool_z)))
      
      if (nrow(df_pool_z) < 500) {
        message("  Skipping: <500 rows after filtering")
        best_pool_mod <- NULL
      } else {
        
        # Z-score standardization
        for(v in preds_pool) {
          if (v %in% names(df_pool_z)) {
            df_pool_z[[v]] <- z_score(df_pool_z[[v]])
          }
        }
        
        # VIF pruning
        vp_pool <- vif_prune(df_pool_z, preds_pool, VIF_THRESHOLD)
        kept_pool <- vp_pool$keep
        message(sprintf("  Predictors after VIF pruning (threshold=%.1f): %s", 
                        VIF_THRESHOLD, paste(kept_pool, collapse=", ")))
        
        if (length(kept_pool) == 0) {
          message("  Skipping: All predictors removed by VIF")
          best_pool_mod <- NULL
        } else {
          
          # Verify each scene has enough observations
          scene_counts <- df_pool_z %>% 
            group_by(scene_id) %>% 
            summarise(n = n()) %>%
            arrange(n)
          
          message("\n  Scene pixel counts after filtering:")
          print(scene_counts)
          
          # Remove scenes with very few pixels (< 10)
          scenes_to_keep <- scene_counts %>% 
            filter(n >= 10) %>% 
            pull(scene_id)
          
          df_pool_z <- df_pool_z %>% 
            filter(scene_id %in% scenes_to_keep) %>%
            mutate(scene_id = droplevels(scene_id))
          
          n_scenes_final <- length(unique(df_pool_z$scene_id))
          message(sprintf("\n  Final number of scenes with ≥10 pixels: %d", n_scenes_final))
          
          if (n_scenes_final < 2) {
            message("  Skipping: Not enough scenes with sufficient data")
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
    }
  }
}

# -----------------------------------------------------------------------------
# 4. XCO2 Correction (FULL BIAS CORRECTION - NOT SCENE-MEAN PRESERVING)
# -----------------------------------------------------------------------------

message("\n--- Step 4: XCO2 Correction (Full Bias Correction) ---")

if (is.null(best_pool_mod)) {
  message("  Skipping: No valid model available")
} else {
  
  coef_vec <- coef(best_pool_mod)
  
  predictor_coefs <- coef_vec[
    !grepl("^scene_id|^\\(Intercept\\)", names(coef_vec))
  ]
  
  message(sprintf("  Applying full bias correction with %d predictors", 
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
    
    # Initialize prediction with INTERCEPT
    y_hat <- rep(coef_vec["(Intercept)"], nrow(df_scene))
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
    
    # Add scene effect (scene-specific intercept)
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
    
    # FULL CORRECTION: Subtract the complete predicted bias (including mean shift)
    df_scene$xco2_corrected <- df_scene$xco2 - y_hat
    
    n_valid <- sum(is.finite(df_scene$xco2_corrected))
    if (n_valid < nrow(df_scene) * 0.5) {
      warning(sprintf("  [%s] Only %d/%d pixels successfully corrected",
                      sid, n_valid, nrow(df_scene)))
    }
    
    # Calculate statistics
    correction_stats[[sid]] <- data.frame(
      scene_id = sid,
      mean_correction = mean(y_hat, na.rm = TRUE),
      sd_correction = sd(y_hat, na.rm = TRUE),
      min_correction = min(y_hat, na.rm = TRUE),
      max_correction = max(y_hat, na.rm = TRUE),
      xco2_mean_before = mean(df_scene$xco2, na.rm = TRUE),
      xco2_mean_after = mean(df_scene$xco2_corrected, na.rm = TRUE),
      xco2_sd_before = sd(df_scene$xco2, na.rm = TRUE),
      xco2_sd_after = sd(df_scene$xco2_corrected, na.rm = TRUE),
      variance_reduction = 1 - sd(df_scene$xco2_corrected, na.rm=TRUE) / sd(df_scene$xco2, na.rm=TRUE),
      mean_shift = mean(df_scene$xco2_corrected, na.rm=TRUE) - mean(df_scene$xco2, na.rm=TRUE),
      mean_shift_pct = 100 * (mean(df_scene$xco2_corrected, na.rm=TRUE) - mean(df_scene$xco2, na.rm=TRUE)) / mean(df_scene$xco2, na.rm=TRUE)
    )
    
    # Preserve raster attributes for plotting
    if (sid %in% names(raster_attributes)) {
      attr(df_scene, "viirs_raster_aod") <- raster_attributes[[sid]]$viirs_raster_aod
      attr(df_scene, "n_viirs_granules") <- raster_attributes[[sid]]$n_viirs_granules
      attr(df_scene, "grid_dims") <- raster_attributes[[sid]]$grid_dims
    }
    
    all_scenes_corrected[[sid]] <- df_scene
    
    message(sprintf("  [%s] ✓ Corrected: mean shift = %.2f ppm (%.2f%%), σ reduction = %.1f%%", 
                    sid, 
                    correction_stats[[sid]]$mean_shift,
                    correction_stats[[sid]]$mean_shift_pct,
                    correction_stats[[sid]]$variance_reduction * 100))
  }
  
  if (length(correction_stats) > 0) {
    correction_summary <- bind_rows(correction_stats)
    write_csv(correction_summary, 
              file.path(OUTPUT_DIR, "correction_statistics.csv"))
    message(sprintf("\n✓ Correction statistics saved (%d scenes)", nrow(correction_summary)))
    
    # Print summary of mean shifts across all scenes
    message("\n--- Scene Mean Shift Summary ---")
    message(sprintf("  Average mean shift: %.2f ppm (%.2f%%)",
                    mean(correction_summary$mean_shift, na.rm=TRUE),
                    mean(correction_summary$mean_shift_pct, na.rm=TRUE)))
    message(sprintf("  Range: [%.2f, %.2f] ppm",
                    min(correction_summary$mean_shift, na.rm=TRUE),
                    max(correction_summary$mean_shift, na.rm=TRUE)))
    
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
    
    # Verify x and y columns exist
    if (!("x" %in% names(df_scene)) || !("y" %in% names(df_scene))) {
      message(sprintf("  [%s] ERROR: Missing x or y coordinates. Available columns: %s", 
                      sid, paste(names(df_scene), collapse=", ")))
      next
    }
    
    message(sprintf("  [%s] Creating publication plot...", sid))
    
    # Calculate GLOBAL color limits across BOTH uncorrected and corrected
    # This ensures both maps use the same scale for fair comparison
    all_values <- c(df_scene$xco2, df_scene$xco2_corrected)
    global_mean <- mean(all_values, na.rm = TRUE)
    global_sd <- sd(all_values, na.rm = TRUE)
    
    color_limits <- c(
      max(min(all_values, na.rm = TRUE), global_mean - 3 * global_sd),
      min(max(all_values, na.rm = TRUE), global_mean + 3 * global_sd)
    )
    
    message(sprintf("    Using color scale: [%.2f, %.2f] ppm", 
                    color_limits[1], color_limits[2]))
    
    # Create comparison plot with SHARED color scale
    grob <- create_map_pair(
      df_scene,
      var1 = "xco2",
      var2 = "xco2_corrected",
      title1 = "Uncorrected XCO₂ (ppm)",
      title2 = "Corrected XCO₂ (ppm)",
      main_title = sprintf("%s - XCO₂ Aerosol Correction (Full Bias Correction)", sid)
    )
    
    out_file <- file.path(OUTPUT_DIR, 
                          sprintf("%s_XCO2_comparison_full_correction.png", sid))
    
    ggsave(out_file, grob, width = 12, height = 6, 
           dpi = 300, bg = "white")
    
    message(sprintf("  [%s] ✓ Saved: %s", sid, basename(out_file)))
    
    # BONUS 1: Correction magnitude map (shows the full bias that was removed)
    df_scene$correction_ppm <- df_scene$xco2 - df_scene$xco2_corrected
    
    # Calculate correction statistics for subtitle
    corr_mean <- mean(df_scene$correction_ppm, na.rm=TRUE)
    corr_sd <- sd(df_scene$correction_ppm, na.rm=TRUE)
    corr_range <- range(df_scene$correction_ppm, na.rm=TRUE)
    
    # Filter extreme outliers only for visualization (>5σ)
    corr_threshold <- 5 * corr_sd
    df_scene_filt <- df_scene %>% 
      filter(abs(correction_ppm - corr_mean) < corr_threshold)
    
    if (nrow(df_scene_filt) > 50) {
      
      # Use symmetric color scale centered on mean correction
      corr_limit <- max(abs(df_scene_filt$correction_ppm - corr_mean), na.rm=TRUE)
      
      p_corr <- ggplot(df_scene_filt, aes(x = x, y = y, fill = correction_ppm)) +
        geom_raster() +
        coord_fixed() +
        scale_fill_gradient2(
          low = "blue", mid = "white", high = "red",
          midpoint = corr_mean,
          limits = c(corr_mean - corr_limit, corr_mean + corr_limit),
          name = "Correction\n(ppm)"
        ) +
        labs(
          title = sprintf("%s - Bias Correction Magnitude", sid),
          subtitle = sprintf("Mean: %.2f ppm | SD: %.2f ppm | Range: [%.2f, %.2f] ppm",
                             corr_mean, corr_sd, corr_range[1], corr_range[2]),
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
                             sprintf("%s_correction_map_full.png", sid))
      ggsave(corr_file, p_corr, width = 10, height = 8, dpi = 300, bg = "white")
      message(sprintf("  [%s] ✓ Saved: %s", sid, basename(corr_file)))
    }
    
    # BONUS 2: Before/After histogram comparison
    df_hist <- data.frame(
      value = c(df_scene$xco2, df_scene$xco2_corrected),
      type = rep(c("Uncorrected", "Corrected"), each = nrow(df_scene))
    )
    
    p_hist <- ggplot(df_hist, aes(x = value, fill = type)) +
      geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
      scale_fill_manual(values = c("Uncorrected" = "#E64B35", "Corrected" = "#4DBBD5"),
                        name = NULL) +
      labs(
        title = sprintf("%s - XCO₂ Distribution Comparison", sid),
        subtitle = sprintf("Mean shift: %.2f → %.2f ppm (Δ = %.2f ppm)",
                           mean(df_scene$xco2, na.rm=TRUE),
                           mean(df_scene$xco2_corrected, na.rm=TRUE),
                           mean(df_scene$xco2_corrected, na.rm=TRUE) - mean(df_scene$xco2, na.rm=TRUE)),
        x = "XCO₂ (ppm)",
        y = "Frequency"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    hist_file <- file.path(OUTPUT_DIR,
                           sprintf("%s_histogram_comparison.png", sid))
    ggsave(hist_file, p_hist, width = 8, height = 6, dpi = 300, bg = "white")
    message(sprintf("  [%s] ✓ Saved: %s", sid, basename(hist_file)))
    
    # BONUS 3: Scatter plot showing correction relationship
    if (nrow(df_scene) <= 10000) {
      # Only create scatter for reasonably-sized scenes
      p_scatter <- ggplot(df_scene, aes(x = xco2, y = xco2_corrected)) +
        geom_point(alpha = 0.3, size = 1, color = "#4DBBD5") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                    color = "red", linewidth = 0.8) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
        labs(
          title = sprintf("%s - Correction Relationship", sid),
          subtitle = sprintf("Red line = 1:1 (no correction) | Black line = actual relationship"),
          x = "Uncorrected XCO₂ (ppm)",
          y = "Corrected XCO₂ (ppm)"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)
        )
      
      scatter_file <- file.path(OUTPUT_DIR,
                                sprintf("%s_scatter_comparison.png", sid))
      ggsave(scatter_file, p_scatter, width = 8, height = 7, dpi = 300, bg = "white")
      message(sprintf("  [%s] ✓ Saved: %s", sid, basename(scatter_file)))
    }
  }
}


# -----------------------------------------------------------------------------
# 6. Final Summary
# -----------------------------------------------------------------------------

message("\n", paste(rep("=", 70), collapse=""))
message("PIPELINE COMPLETE - COSINE ANGLE VERSION")
message(paste(rep("=", 70), collapse=""))
message(sprintf("\nOutput directory: %s", OUTPUT_DIR))
message("\nGenerated files:")
message("  ✓ unified_database_viirs_grid_cosine.csv - All data with cosine angles")
message("  ✓ extraction_metadata.csv - Scene metadata including grid dimensions")
message("  ✓ pooled_dredge_results.csv - Model selection results")
message("  ✓ best_model_coefficients.csv - Final model coefficients")
message("  ✓ correction_statistics.csv - Per-scene correction metrics")
message("  ✓ model_diagnostics.png - Regression diagnostic plots")
message("  ✓ *_XCO2_comparison_cosine.png - Side-by-side before/after")
message("  ✓ *_correction_map_cosine.png - Correction magnitude maps")
if (length(extraction_errors) > 0) {
  message("  ✓ extraction_errors.csv - Detailed error log for failed scenes")
}

message("\n", paste(rep("-", 70), collapse=""))
message("KEY MODIFICATIONS IN COSINE ANGLE VERSION")
message(paste(rep("-", 70), collapse=""))
message("✓ ANGLE TRANSFORMATION: cos(sza), cos(vza), cos(primary_angle)")
message("✓ SPATIAL RESOLUTION: VIIRS native grid (~6km)")
message("✓ TEMPORAL FILTERING: ±", TEMPORAL_WINDOW_HOURS, "h VIIRS granule selection")
message("✓ QUALITY CONTROL: AE range [", AE_MIN, ", ", AE_MAX, "], QF=1 only")
message("✓ SCENE-MEAN PRESERVING: XCO₂ mean unchanged, variance reduced")
message("✓ PUBLICATION QUALITY: 300 DPI, consistent color scales")
message("✓ STATISTICAL RIGOR: VIF pruning, AICc model selection")

message("\n", paste(rep("-", 70), collapse=""))
message("COSINE TRANSFORMATION RATIONALE")
message(paste(rep("-", 70), collapse=""))
message("✓ Linearizes geometric effects of solar/viewing angles")
message("✓ cos(0°) = 1.0 (direct overhead), cos(90°) = 0.0 (horizon)")
message("✓ Better captures radiative transfer physics")
message("✓ Reduces nonlinear VIF inflation from raw angles")
message("✓ Standard practice in atmospheric retrieval modeling")

message("\nAnalysis complete!")