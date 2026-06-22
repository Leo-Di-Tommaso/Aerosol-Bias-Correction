# =============================================================================
# Unified MethaneSAT + VIIRS XCO2 Modeling Pipeline
# =============================================================================
# This script:
# 1. Discovers MethaneSAT files (specific list or all recursive).
# 2. Downloads/Finds matching VIIRS aerosol data.
# 3. Extracts, Rasterizes, and Colocates data into a single UNIFIED DATABASE.
# 4. Runs Unified Statistical Modeling (AICc, PLS, accounting for scene effects).
# TO RUN THIS SCRIPT: source ~/earthdata-env/bin/activate


# -----------------------------------------------------------------------------
# 0. Libraries and Constants
# -----------------------------------------------------------------------------
library(terra)
library(ncdf4)
library(tidyr)
library(readr)
library(writexl)
library(MASS) # stepAIC, rlm
library(broom)  # tidy
library(purrr)
library(jsonlite)
library(stringr)
library(caret)
library(MuMIn)  # dredge, AICc
library(pls)  # plsr
library(car)  # vif
library(speedglm)
library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)
library(dplyr)


# --- CONFIGURATION ---


# >>> USER TOGGLE: Run specific files only? <<<
USE_SPECIFIC_FILES <- TRUE # Set FALSE to search MSAT_ROOT_DIR recursively


# Directories
# NOTE: These paths are specific to your environment and likely rely on
# Google Cloud Storage (GCS) mounts or local directories.
MSAT_ROOT_DIR <- "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results"
VIIRS_DATA_DIR <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/viirs_data2"
# Using the hardcoded path for the Python script as provided in the snippet
PYTHON_SCRIPT_PATH <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/download_viirs_deep_blue2.py"
OUTPUT_DIR <- paste0("results_unified_model_", format(Sys.Date(), "%Y%m%d"))


# Modeling Constants
VIF_THRESHOLD <- 5
TOP_N_DREDGE <- 10


# Create directories
if (!dir.exists(VIIRS_DATA_DIR)) dir.create(VIIRS_DATA_DIR, recursive = TRUE)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. Helper Functions
# -----------------------------------------------------------------------------

# --- A. VIIRS Acquisition ---
# Extracts bounding box and *median* time/date from an MSAT NetCDF file for VIIRS search.
get_bbox_and_date <- function(msat_file) {
  nc <- try(nc_open(msat_file), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)
 
  lon <- try(ncvar_get(nc, "geolocation/longitude"), silent = TRUE)
  lat <- try(ncvar_get(nc, "geolocation/latitude"), silent = TRUE)
  time_vals <- try(ncvar_get(nc, "geolocation/time"), silent = TRUE)
  nc_close(nc)
 
  if (inherits(lon, "try-error") || inherits(lat, "try-error") || inherits(time_vals, "try-error")) return(NULL)


  origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
  # UPDATE: Use median time for better VIIRS matchup, as requested
  mid_time <- origin + median(time_vals, na.rm = TRUE)
 
  list(
    bbox = c(min(lon, na.rm=TRUE), min(lat, na.rm=TRUE), max(lon, na.rm=TRUE), max(lat, na.rm=TRUE)),
    date = as.character(as.Date(mid_time))
  )
}


find_or_download_viirs <- function(msat_file) {

  script_path <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/download_viirs_deep_blue2.py"

  scene_info <- get_bbox_and_date(msat_file)
  if (is.null(scene_info)) return(NA_character_)

  # --- Fix longitude wrap for CMR (critical!) ---
  bbox <- scene_info$bbox
  lon1 <- ifelse(bbox[1] > 180, bbox[1] - 360, bbox[1])
  lon2 <- ifelse(bbox[3] > 180, bbox[3] - 360, bbox[3])
  bbox_clean <- c(lon1, bbox[2], lon2, bbox[4])

  # --- Build search pattern ---
  date_obj <- as.Date(scene_info$date)
  year_doy <- format(date_obj, "%Y%j")
  search_pattern <- sprintf(
    "AERDB_L2_VIIRS_(SNPP|NOAA20|NOAA21).*A%s.*\\.(nc|nc4)$",
    year_doy
  )

  cat(sprintf(" -> Searching for local VIIRS file (Date %s)\n", scene_info$date))

  # --- Search local directory ---
  local_files <- list.files(
    VIIRS_DATA_DIR,
    pattern = search_pattern,
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(local_files) > 0) {
    cat(sprintf(" -> Found existing VIIRS file: %s\n", basename(local_files[1])))
    return(local_files[1])
  }

  # =====================================================================
  # 2. DOWNLOAD IF MISSING
  # =====================================================================

  cat(" -> No local file found. Attempting to download via Python script...\n")

  if (!file.exists(script_path)) {
    warning(sprintf("Downloader script not found at: %s", script_path))
    return(NA_character_)
  }

  # Build argument vector EXACTLY as Python expects:
  # python3 script --bbox W S E N --date yyyy-mm-dd --outdir dir --platforms NOAA20,NOAA21
  cmd_args <- c(
    script_path,
    "--bbox", as.character(bbox_clean[1]),
    as.character(bbox_clean[2]),
    as.character(bbox_clean[3]),
    as.character(bbox_clean[4]),
    "--date", scene_info$date,
    "--outdir", VIIRS_DATA_DIR,
    "--platforms", "NOAA20,NOAA21"
  )

  # Run python
  res <- tryCatch(
    system2("python3", args = cmd_args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      warning("Download script execution failed: ", conditionMessage(e))
      return(NULL)
    }
  )

  # If system2 returned NULL or signaled an error
  if (is.null(res)) {
    warning("Download script returned NULL (execution failed).")
    return(NA_character_)
  }

  cat(" -> Python Script Output:\n")
  cat(paste(res, collapse = "\n"), "\n")

  # =====================================================================
  # 3. CHECK FOR DOWNLOADED FILES
  # =====================================================================

  cat(" -> Rechecking for VIIRS file...\n")

  downloaded <- list.files(
    VIIRS_DATA_DIR,
    pattern = search_pattern,
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(downloaded) > 0) {
    cat(sprintf(" -> Successfully found downloaded file: %s\n",
                basename(downloaded[1])))
    return(downloaded[1])
  }

  cat(" -> VIIRS file could not be found or downloaded.\n")
  return(NA_character_)
}


# --- B. Spatial Processing ---
# Calculates the scattering/phase angle (often called Primary Angle in atmospheric remote sensing).
calc_primary_angle <- function(sza, vza, saa, vaa) {
  dtor <- pi / 180
  rtod <- 180 / pi
  # Formula for the cosine of the scattering angle (phase angle)
  cos_phase <- cos(sza*dtor) * cos(vza*dtor) +
               sin(sza*dtor) * sin(vza*dtor) * cos((saa - vaa)*dtor)
 # Clamp cosine values to [-1, 1] to avoid NaN from floating point errors
 # and convert back to degrees.
  acos(pmax(-1, pmin(1, cos_phase))) * rtod
}


extract_scene_data <- function(viirs_file, msat_file, scene_id) {
  
  # ---------------------------------------------
  # 1. LOAD & MOSAIC VIIRS GRANULES FROM SAME DAY
  # ---------------------------------------------

  scene_info <- get_bbox_and_date(msat_file)
  scene_date <- as.Date(scene_info$date)

  # List all VIIRS granules in VIIRS_DATA_DIR
  viirs_files_all <- list.files(VIIRS_DATA_DIR, pattern="\\.nc$", full.names=TRUE)

  # Extract date from filename robustly: AYYYYDDD
  viirs_files_day <- Filter(function(f) {
    fn <- basename(f)
    match <- regmatches(fn, regexpr("A(\\d{7})", fn))
    if(length(match) == 0) return(FALSE)
    yyyydoy <- substr(match, 2, 8)  # "2025120"
    year <- as.integer(substr(yyyydoy,1,4))
    doy  <- as.integer(substr(yyyydoy,5,7))
    granule_date <- as.Date(doy - 1, origin = paste0(year,"-01-01"))
    !is.na(granule_date) && granule_date == scene_date
  }, viirs_files_all)

  if(length(viirs_files_day) == 0) stop(sprintf("No VIIRS granules found for scene %s on %s", scene_id, scene_date))

  # Load and mosaic AOD and AE rasters
  load_viirs_raster <- function(file, varname) {
    nc <- nc_open(file)
    lon <- ncvar_get(nc, "Longitude")
    lat <- ncvar_get(nc, "Latitude")
    vals <- ncvar_get(nc, varname)
    nc_close(nc)
    rast(nrows=nrow(vals), ncols=ncol(vals),
        xmin=min(lon), xmax=max(lon),
        ymin=min(lat), ymax=max(lat),
        vals=as.vector(vals))
  }

  v_raster_aod <- do.call(terra::merge, lapply(viirs_files_day, load_viirs_raster, varname="Aerosol_Optical_Thickness_550_Land_Best_Estimate"))
  v_raster_ae  <- do.call(terra::merge, lapply(viirs_files_day, load_viirs_raster, varname="Angstrom_Exponent_Land_Best_Estimate"))

  crs(v_raster_aod) <- "EPSG:4326"
  crs(v_raster_ae)  <- "EPSG:4326"



  # ---------------------------------------------
  # 2. LOAD MSAT POINTS
  # ---------------------------------------------
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
  
  msat_pts <- vect(msat_df, geom=c("lon","lat"), crs="EPSG:4326")
  
  # ---------------------------------------------
  # 3. EXTRACT VIIRS VALUES AT MSAT POINTS
  # ---------------------------------------------
  ## DEBUG
  print(terra::ext(v_raster_aod))  # VIIRS raster extent
  print(range(msat_df$lon, na.rm=TRUE))  # MSAT lon range
  print(range(msat_df$lat, na.rm=TRUE))  # MSAT lat range

  aod_vals <- terra::extract(v_raster_aod, msat_pts)
  ae_vals  <- terra::extract(v_raster_ae, msat_pts)

  
  msat_df$AOD550 <- aod_vals[,2]  # terra returns ID in first column
  msat_df$AE     <- ae_vals[,2]
  
  # Filter points with valid VIIRS data
  msat_df <- msat_df %>% filter(!is.na(AOD550) & !is.na(AE))
  
  if (nrow(msat_df) < 50)
    stop("Insufficient overlapping pixels")
  
  # ---------------------------------------------
  # 4. DERIVED VARIABLES
  # ---------------------------------------------
  msat_df$AOD_1286 <- msat_df$AOD550 * (1286/550)^(-msat_df$AE)
  msat_df$primary_angle <- calc_primary_angle(msat_df$sza, msat_df$vza, msat_df$saa, msat_df$vaa)
  msat_df$scene_id <- scene_id
  
  return(msat_df)
}



# --- C. Stats Helpers ---
# Standard Z-score normalization
z_score <- function(x) as.numeric(scale(x))


# Calculates the maximum Variance Inflation Factor (VIF) for a set of predictors.
max_vif_num <- function(df, xvars) {
  if (length(xvars) < 2) return(1.0)
  form <- as.formula(paste("xco2 ~", paste(xvars, collapse="+")))
  # Use a subsample for faster VIF calculation on large datasets
  dat_fit <- if(nrow(df) > 5000) df[sample(nrow(df), 5000), ] else df
  tryCatch({
    mod <- lm(form, data = dat_fit)
    vifs <- car::vif(mod)
    # VIF output can be a vector (if one predictor) or matrix
    if(is.null(dim(vifs))) max(vifs) else max(vifs[,1])
  }, error = function(e) NA_real_)
}


# Iteratively removes the variable with the highest VIF until all remaining VIFs are below the threshold.
vif_prune <- function(df, xvars, threshold=5) {
  keep <- xvars
  repeat {
    if (length(keep) < 2) break
    curr_max <- max_vif_num(df, keep)
    if (is.na(curr_max) || curr_max <= threshold) break


    form <- as.formula(paste("xco2 ~", paste(keep, collapse="+")))
    # Use a subsample for calculation
    dat_fit <- if(nrow(df) > 5000) df[sample(nrow(df), 5000), ] else df
    mod <- lm(form, data=dat_fit)
    v <- car::vif(mod)
    val_to_check <- if(is.null(dim(v))) v else v[,1]
    # Remove the variable with the highest VIF
    keep <- setdiff(keep, names(which.max(val_to_check)))
  }
  list(keep=keep)
}


# -----------------------------------------------------------------------------
# 2. Main: Discovery & Extraction
# -----------------------------------------------------------------------------
message("--- Step 1: File Discovery ---")
  
if (USE_SPECIFIC_FILES) {
  message("-> Mode: Specific List")
  all_nc_files <- c(
    # List of specific MethaneSAT files to process
  
  #   "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t348/20250430/c0A9D15C0/p3665/MSAT_L2_c0A9D15C0_p3665_v04001006_20250430T132730Z_132802Z.nc",
  #   "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t108/20250120/c06B506C0/p3188/MSAT_L2_c06B506C0_p3188_v04001003_20250120T071415Z_071447Z.nc",
  #   "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t2/20240905/c010C0020/p3262/MSAT_L2_c010C0020_p3262_v04001003_20240905T194216Z_194248Z.nc",
  #   "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t60/20241228/c056403C0/p3860/MSAT_L2_c056403C0_p3860_v04001006_20241228T210922Z_210954Z.nc",
  #   "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t56/20250320/c09060380/p3221/MSAT_L2_c09060380_p3221_v04001003_20250320T210203Z_210235Z.nc",
  #   "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t147/20250606/c0CF30930/p3651/MSAT_L2_c0CF30930_p3651_v04001006_20250606T115534Z_115606Z.nc",
  #   "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t25/20241122/c046E0190/p3212/MSAT_L2_c046E0190_p3212_v04001003_20241122T124404Z_124436Z.nc",
  #   "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t25/20241220/c04D90190/p3232/MSAT_L2_c04D90190_p3232_v04001003_20241220T122155Z_122227Z.nc"
  # )

    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t349/20250516/c0BC915D0/p3838/MSAT_L2_c0BC915D0_p3838_v04001006_20250516T212530Z_212602Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t349/20250119/c069C15D0/p4210/MSAT_L2_c069C15D0_p4210_v04003001_20250119T213734Z_213806Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t348/20250220/c073C15C0/p3678/MSAT_L2_c073C15C0_p3678_v04001006_20250220T132419Z_132451Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t348/20250430/c0A9D15C0/p3665/MSAT_L2_c0A9D15C0_p3665_v04001006_20250430T132730Z_132802Z.nc")

  all_nc_files <- all_nc_files[file.exists(all_nc_files)]
} else {
  message("-> Mode: Recursive Search")
  all_nc_files <- list.files(MSAT_ROOT_DIR, recursive = TRUE, pattern = "\\.nc$", full.names = TRUE)
}


if (length(all_nc_files) == 0) stop("No valid MSAT files found.")


# ID Extraction: Creates a unique scene identifier from the file path
extract_scene_n <- function(path) {
  t_m <- str_extract(path, "t[0-9]+")
  p_m <- str_extract(path, "p[0-9]+")
  if (!is.na(t_m) && !is.na(p_m)) return(paste0(t_m, p_m))
  if (!is.na(t_m)) return(t_m)
  return(basename(dirname(path)))
}


# Filter to unique scenes based on extracted ID
files_df <- data.frame(MSAT_file = all_nc_files, stringsAsFactors=FALSE) %>%
  mutate(scene_id = sapply(MSAT_file, extract_scene_n)) %>%
  distinct(scene_id, .keep_all=TRUE)


message(sprintf("Found %d unique scenes.", nrow(files_df)))


message("\n--- Step 2: Extraction & Unification (Serial Mode) ---")


# A. Download VIIRS (Pre-pass)
files_df$VIIRS_file <- sapply(files_df$MSAT_file, find_or_download_viirs)
files_df <- files_df %>% filter(!is.na(VIIRS_file))
message(sprintf("Scenes with available VIIRS data: %d", nrow(files_df)))


# B. Extract Data
valid_data_list <- list()


for (i in 1:nrow(files_df)) {
  sid <- files_df$scene_id[i]
  cat(sprintf("[%d/%d] Processing %s ... ", i, nrow(files_df), sid))
 
  res <- tryCatch({
    extract_scene_data(files_df$VIIRS_file[i], files_df$MSAT_file[i], sid)
  }, error = function(e) {
    # Custom error handling for user feedback
    if (grepl("Insufficient overlapping pixels", e$message)) {
      cat(sprintf("FAILED (Insufficient Pixels). Error: %s\n", e$message))
    } else {
      cat(sprintf("FAILED (Extraction Error). Error: %s\n", e$message))
    }
    return(NULL)
  })
 
  if (!is.null(res)) {
    valid_data_list[[length(valid_data_list)+1]] <- res
    cat(sprintf("SUCCESS (%d pixels)\n", nrow(res)))
  }
}


# C. Create Unified Database
if (length(valid_data_list) == 0) stop("No data extracted successfully.")


all_scenes_df <- bind_rows(valid_data_list)


# Save Unified DB
db_path <- file.path(OUTPUT_DIR, "unified_database.csv")
write_csv(all_scenes_df, db_path)
message(sprintf("\nUnified Database created with %d rows. Saved to: %s", nrow(all_scenes_df), db_path))




# -----------------------------------------------------------------------------
# 3. Statistical Modeling (Unified/Pooled)
# -----------------------------------------------------------------------------
message("\n--- Step 3: Statistical Modeling (Unified/Pooled) ---")


# List of all potential predictors (MethaneSAT and VIIRS derived)
predictors_all <- c("dp", "alb1286", "alb1606", "psfc", "orog", "sza", "vza", "primary_angle", "AOD_1286", "AE")


if (nrow(all_scenes_df) < 500) {
  message("  Skipping pooled analysis: Insufficient total pixels (<500).")
} else {
 
  # Select only necessary columns and remove any remaining NA values
  target_cols <- c("xco2", "scene_id", predictors_all)
 
  df_pool_clean <- all_scenes_df %>%
    dplyr::select(any_of(target_cols)) %>%
    drop_na()
 
  # Final check: Dependent variable must exist
  if (!("xco2" %in% names(df_pool_clean))) {
      stop("Fatal Error: The dependent variable 'xco2' is missing after data cleaning.")
  }
 
  # Prepare data for modeling
  df_pool_clean$scene_id <- as.factor(df_pool_clean$scene_id)
  preds_pool <- intersect(names(df_pool_clean), predictors_all)
 
  if (length(preds_pool) < 2) {
    message("  Skipping pooled analysis: Less than 2 valid predictors.")
  } else {
   
    # 1. Z-Score predictors for stable regression coefficients
    df_pool_z <- df_pool_clean
    for(v in preds_pool) df_pool_z[[v]] <- z_score(df_pool_z[[v]])
   
    # 2. VIF Pruning (addressing multicollinearity)
    vp_pool <- vif_prune(df_pool_z, preds_pool, VIF_THRESHOLD)
    kept_pool <- vp_pool$keep
    message("  Pooled Predictors (Post-VIF): ", paste(kept_pool, collapse=", "))
   
    if (length(kept_pool) == 0) {
      message("  Skipping pooled analysis: All predictors removed by VIF pruning.")
    } else {
     
      # # --- A.1 AICc Best Subset (Dredge) with Scene ID as a fixed factor ---
      message("\n  -> Running AICc Best Subset (Dredge)...")
      # Global model including all kept predictors and fixed 'scene_id' factor
      form_pool <- as.formula(paste("xco2 ~ scene_id +", paste(kept_pool, collapse="+")))
      mod_pool_global <- lm(form_pool, data=df_pool_z, na.action="na.fail")
     
      # Dredge: Explores all combinations of predictors, keeping scene_id always present
      dredge_pool <- dredge(mod_pool_global, rank="AICc", fixed=c("scene_id"), trace=TRUE, m.lim=c(1, length(kept_pool)))
     
      # Save top models
      write_csv(as.data.frame(dredge_pool), file.path(OUTPUT_DIR, "pooled_dredge_results.csv"))
     
      best_pool_mod <- get.models(dredge_pool, subset=1)[[1]]
      r2_adj <- summary(best_pool_mod)$adj.r.squared
      message(sprintf("  Pooled MLR Adjusted R-squared: %.3f", r2_adj))
    }
  }
}

message("\nDone. Check output folder: ", OUTPUT_DIR)

      #       # --- A.2 Fit Unified Regression Using speedglm (Big-Data Safe) --- FOR HUGE DATASETS (the whole level2_results)
      # message("\n  -> Running speedglm Unified Regression...")

      # # Form using kept predictors only (scene_id is excluded: factors not supported efficiently in speedglm)
      # form_speed <- as.formula(
      #   paste("xco2 ~", paste(kept_pool, collapse = " + "))
      # )

      # # Fit model using speedglm (handles 40M rows safely)
      # mod_speed <- speedglm(
      #   formula = form_speed,
      #   data = df_pool_z
      # )

      # message("  speedglm model fit complete.")
      # print(summary(mod_speed))

      # # Extract Adjusted R² manually (speedglm does not compute it)
      # # R² = 1 - SSE/SST
      # preds <- mod_speed$fitted.values
      # SSE <- sum((df_pool_z$xco2 - preds)^2)
      # SST <- sum((df_pool_z$xco2 - mean(df_pool_z$xco2))^2)
      # r2_speed <- 1 - SSE/SST

      # message(sprintf("  speedglm R-squared: %.3f", r2_speed))
      
      # Save coefficient table
      # coef_table <- data.frame(
      #   term = names(mod_speed$coefficients),
      #   estimate = mod_speed$coefficients
      # )
      # write_csv(coef_table, file.path(OUTPUT_DIR, "speedglm_pooled_coefficients.csv"))



#       # --- B. PLS Regression (Alternative modeling technique) ---
#       message("\n  -> Running PLS Regression...")
#       # PLS is typically run on continuous predictors only
#       pls_data <- df_pool_z %>% dplyr::select(xco2, any_of(kept_pool))
#       # Cross-validated PLS to determine optimal number of components
#       pls_res <- plsr(xco2 ~ ., data=pls_data, validation="CV")
#       ncomp_opt <- which.min(pls_res$validation$PRESS)
#       rmsep_opt <- min(RMSEP(pls_res)$val[1,,])
     
#       message(sprintf("  PLS Optimal Components: %d, RMSEP: %.3f", ncomp_opt, rmsep_opt))
     
#       # Save main results summary
#       results_summary <- tibble(
#         analysis = "Pooled_Unified",
#         n_pixels = nrow(df_pool_z),
#         adj_r2_mlr = r2_adj,
#         pls_ncomp = ncomp_opt,
#         pls_rmsep = rmsep_opt,
#         kept_vars_vif = paste(kept_pool, collapse=", ")
#       )
#       write_csv(results_summary, file.path(OUTPUT_DIR, "unified_model_summary.csv"))
# =============================================================================
# Step 4: XCO2 Correction Using POOLED Best AICc Model
# =============================================================================

message("\n--- Step 4: XCO2 Correction Using Pooled AICc Model ---")

# Extract best pooled model coefficients ONCE
coef_vec <- coef(best_pool_mod)

# Separate predictor coefficients (exclude intercept + scene effects)
predictor_coefs <- coef_vec[
  !grepl("^scene_id", names(coef_vec)) &
  names(coef_vec) != "(Intercept)"
]

for (sid in unique(all_scenes_df$scene_id)) {

  df_scene <- all_scenes_df %>% filter(scene_id == sid)
  if (nrow(df_scene) < 50) next

  # Compute predictor-driven bias
  y_hat <- rep(0, nrow(df_scene))
  for (v in names(predictor_coefs)) {
    if (v %in% names(df_scene)) {
      y_hat <- y_hat + predictor_coefs[[v]] * scale(df_scene[[v]])
    }
  }

  # Scene-mean preserving correction
  df_scene <- df_scene %>%
    mutate(
      y_hat_m = y_hat,
      xco2_corrected = xco2 - (y_hat_m - mean(y_hat_m, na.rm = TRUE))
    )

  message(sprintf("[%s] XCO2 correction complete (%d pixels).",
                  sid, nrow(df_scene)))
}

# ---------------------------------------------------------------------------
# 5. Raster-based XCO2 Comparison (Uncorrected vs Corrected)
# ---------------------------------------------------------------------------
create_map_pair <- function(df, var1, var2, title1, title2, main_title) {

  r1 <- ggplot(df, aes(x, y, fill = .data[[var1]])) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis(option = "C") +
    labs(title = title1, fill = "ppm") +
    theme_minimal()

  r2 <- ggplot(df, aes(x, y, fill = .data[[var2]])) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis(option = "C") +
    labs(title = title2, fill = "ppm") +
    theme_minimal()

  gridExtra::grid.arrange(r1, r2, ncol = 2,
                          top = grid::textGrob(main_title,
                                gp = grid::gpar(fontsize = 14, fontface = "bold")))
}

plot_xco2_comparison <- function(df_scene, xco2_corrected, scene_id, out_dir) {

  stopifnot(nrow(df_scene) == length(xco2_corrected))

  # -----------------------------------------------------------------------
  # 1. Attach corrected XCO2 to scene dataframe
  # -----------------------------------------------------------------------
  df_scene$xco2_corrected <- xco2_corrected

  # Drop rows with missing values (important for rasterization)
  df_scene <- df_scene %>%
    filter(is.finite(xco2), is.finite(xco2_corrected))

  if (nrow(df_scene) < 50) {
    warning(sprintf("[%s] Too few valid pixels for plotting", scene_id))
    return(NULL)
  }

  # -----------------------------------------------------------------------
  # 2. Convert to SpatVector
  # -----------------------------------------------------------------------
  v_scene <- vect(
    df_scene,
    geom = c("lon", "lat"),
    crs  = "EPSG:4326"
  )

  # -----------------------------------------------------------------------
  # 3. Build raster template from scene extent
  #    (matches your VIIRS/MSAT grid logic)
  # -----------------------------------------------------------------------
  r_template <- rast(
    xmin = min(df_scene$lon), xmax = max(df_scene$lon),
    ymin = min(df_scene$lat), ymax = max(df_scene$lat),
    ncols = 500, nrows = 500,   # reasonable visualization resolution
    crs = "EPSG:4326"
  )


  # -----------------------------------------------------------------------
  # 4. Rasterize XCO2 fields
  # -----------------------------------------------------------------------
  r_xco2_raw <- rasterize(
    v_scene, r_template,
    field = "xco2",
    fun = "mean", na.rm = TRUE
  )

  r_xco2_corr <- rasterize(
    v_scene, r_template,
    field = "xco2_corrected",
    fun = "mean", na.rm = TRUE
  )

  # -----------------------------------------------------------------------
  # 5. Convert rasters back to dataframe for ggplot
  # -----------------------------------------------------------------------
  df_plot <- as.data.frame(
    c(r_xco2_raw, r_xco2_corr),
    xy = TRUE, na.rm = TRUE
  )

  names(df_plot) <- c("x", "y", "xco2_raw", "xco2_corrected")

  # -----------------------------------------------------------------------
  # 6. Plot side-by-side maps (same style as your pipeline)
  # -----------------------------------------------------------------------
  p <- create_map_pair(
    df_plot,
    var1 = "xco2_raw",
    var2 = "xco2_corrected",
    title1 = "Uncorrected XCO2 (ppm)",
    title2 = "Corrected XCO2 (ppm)",
    main_title = paste(scene_id, "- XCO2 Correction Comparison")
  )

  out_file <- file.path(out_dir, paste0(scene_id, "_XCO2_raw_vs_corrected.png"))
  ggsave(out_file, p, width = 10, height = 5, bg = "white")

  message(sprintf("[%s] Step 5 complete: raster XCO2 comparison saved.", scene_id))
}
message("\n--- Step 5: Raster-based XCO2 Comparison ---")

for (sid in unique(all_scenes_df$scene_id)) {

  df_scene <- all_scenes_df %>% filter(scene_id == sid)
  if (nrow(df_scene) < 50) next

  xco2_corrected <- df_scene$xco2 - df_scene$dp  # matches Step 4 logic

  plot_xco2_comparison(
    df_scene       = df_scene,
    xco2_corrected = xco2_corrected,
    scene_id       = sid,
    out_dir        = OUTPUT_DIR
  )
}
