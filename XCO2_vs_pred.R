# IMPORTANT: to run the following files (this and "download_viirs_deep_blue.py"), you must first activate the earthdata environment with the following command:
# source ~/earthdata-env/bin/activate

# -----------------------------------------------------------------------------
# Step 0: Load Libraries and Define Constants
# -----------------------------------------------------------------------------
# Load all necessary libraries
library(terra)
library(ncdf4)
library(dplyr)
library(tidyr)
library(viridis)
library(ggplot2)
library(MASS)
library(broom)
library(purrr)      # Used for iterating and combining data
library(jsonlite)
library(stringr)    # Used for robust path parsing
library(caret)      # Added for multicollinearity check
library(parallel)   # Added for parallel processing

# --- Constants ---
# Root directory for MethaneSAT Level 2 data
MSAT_ROOT_DIR <- "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results"

# Root directory where VIIRS data is stored or will be downloaded.
# PLEASE UPDATE THIS PATH if it's different on your system.
VIIRS_DATA_DIR <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/viirs_data"

# Directory to save plots and results
PLOT_DIR <- "by_aerosol_data"
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE) ## DEBUG: Added recursive = TRUE for safety.

# -----------------------------------------------------------------------------
# Step 1: Discover MethaneSAT Files and Build Scene List
# -----------------------------------------------------------------------------
message("--- Step 1: Discovering MethaneSAT .nc files ---")

# --- USER TOGGLE: Set to TRUE to run specific files, FALSE to run all files ---
USE_SPECIFIC_FILES <- TRUE
# ---

if (USE_SPECIFIC_FILES) {
  # Option 1: Use a predefined list of specific MethaneSAT files.
  message("-> Using a predefined list of specific scenes.")

  all_nc_files <- c(
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t26/20241124/c049501A0/p748/MSAT_L2_c049501A0_p748_v03002000_20241124T143538Z_143610Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t380/20250319/c08F617C0/p1524/MSAT_L2_c08F617C0_p1524_v03005001_20250319T080703Z_080735Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/t159/2024/09/12/014E09F0/202409171319_po-2103/level2-results/MSAT_L2_post_20240912T221258_20240912T221330_post.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t25/20241119/c04350190/p898/MSAT_L2_c04350190_p898_v03002000_20241119T122526Z_122558Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t381/20250426/c0A4217D0/p1618/MSAT_L2_c0A4217D0_p1618_v03005001_20250426T144316Z_144348Z.nc"
  )

  # Verify that the specified files exist before proceeding
  existing_files <- file.exists(all_nc_files)
  if (any(!existing_files)) {
    warning("The following specified MethaneSAT files were not found and will be skipped:")
    print(all_nc_files[!existing_files])
  }
  all_nc_files <- all_nc_files[existing_files]

} else {
  # Option 2: Find all .nc files recursively under the root directory.
  message("-> Searching for all .nc files in ", MSAT_ROOT_DIR)

  all_nc_files <- list.files(
    path = MSAT_ROOT_DIR,
    recursive = TRUE,
    pattern = "\\.nc$",
    full.names = TRUE
  )
}


if (length(all_nc_files) == 0) {
  stop("No MethaneSAT .nc files found to process.")
}

# Create a dataframe from the file paths
files <- data.frame(MSAT_file_nc = all_nc_files, stringsAsFactors = FALSE)

# Function to extract a scene identifier (e.g., "t26p748") from a file path
extract_scene_n <- function(path) {
  # Look for t### and p### patterns in the file path
  target_match <- str_extract(path, "t[0-9]+")
  plot_match <- str_extract(path, "p[0-9]+")

  if (!is.na(target_match) && !is.na(plot_match)) {
    return(paste0(target_match, plot_match))
  } else if (!is.na(target_match)) {
    return(target_match)
  } else {
    # Fallback to a directory name if patterns are not found
    parts <- strsplit(path, "/")[[1]]
    if (length(parts) >= 9) {
        return(parts[9])
    }
    return(basename(dirname(path)))
  }
}

# Populate the scene_n column and remove any duplicate scenes
files <- files %>%
  mutate(scene_n = sapply(MSAT_file_nc, extract_scene_n)) %>%
  distinct(scene_n, .keep_all = TRUE)

message(sprintf("Found %d unique MethaneSAT scenes to process.", nrow(files)))

# -----------------------------------------------------------------------------
# Step 2: Acquire VIIRS Data for All Scenes
# -----------------------------------------------------------------------------
message("\n--- Step 2: Acquiring VIIRS data for all scenes ---")

# Ensure the VIIRS download directory exists
if (!dir.exists(VIIRS_DATA_DIR)) {
  message("Creating VIIRS data directory: ", VIIRS_DATA_DIR)
  dir.create(VIIRS_DATA_DIR, recursive = TRUE)
}

# Helper function to get bounding box and date from a MethaneSAT file
get_bbox_and_date <- function(msat_file) {
  nc <- nc_open(msat_file)

  lon <- ncvar_get(nc, "geolocation/longitude")
  lat <- ncvar_get(nc, "geolocation/latitude")
  time_vals <- ncvar_get(nc, "geolocation/time")

  nc_close(nc)

  origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
  start_time <- origin + min(time_vals, na.rm = TRUE)

  list(
    bbox = c(min(lon, na.rm = TRUE), min(lat, na.rm = TRUE), max(lon, na.rm = TRUE), max(lat, na.rm = TRUE)),
    date = as.character(as.Date(start_time))
  )
}

# Function to find a local VIIRS file or download it if not found
find_or_download_viirs <- function(msat_file) {
  scene_info <- get_bbox_and_date(msat_file)

  # 1. Search for existing local file
  date_obj <- as.Date(scene_info$date)
  year_doy <- format(date_obj, "%Y%j")
  search_pattern <- paste0("AERDB_L2_VIIRS_SNPP.A", year_doy, ".*\\.nc$")

  cat("  -> Searching for local VIIRS file with pattern:", search_pattern, "\n")
  local_files <- list.files(VIIRS_DATA_DIR, pattern = search_pattern, full.names = TRUE, recursive = TRUE)

  if (length(local_files) > 0) {
    cat("  -> Found local file:", basename(local_files[1]), "\n")
    return(local_files[1]) # Return the first match
  }

  # 2. If not found, attempt to download
  cat("  -> No local file found. Attempting to download...\n")

  cmd_args <- c(
    "download_viirs_deep_blue.py",
    "--date", scene_info$date,
    "--outdir", VIIRS_DATA_DIR,
    "--bbox", as.character(scene_info$bbox)
  )

  # Execute the python script
  res <- tryCatch(
    system2("python", args = cmd_args, stdout = TRUE, stderr = TRUE),
    error = function(e) e
  )

  # Check for execution errors or if the script returned no files
  # A non-zero status attribute often indicates an error from the script itself.
  if (inherits(res, "error") || !is.null(attr(res, "status"))) {
    cat("  -> FAILED to download. Error:\n", paste(res, collapse="\n"), "\n")
    return(NA_character_)
  }

  # Parse the JSON output from the python script
  viirs_files <- tryCatch(fromJSON(paste(res, collapse = "")), error = function(e) NA)

  if (any(is.na(viirs_files)) || length(viirs_files) == 0) {
    cat("  -> Download failed or no VIIRS files found for this date/location.\n")
    return(NA_character_)
  }

  cat("  -> Successfully downloaded:", basename(viirs_files[[1]]), "\n")
  return(viirs_files[[1]])
}

# Apply the function to get all required VIIRS file paths
files$VIIRS_file_nc <- purrr::map_chr(files$MSAT_file_nc, find_or_download_viirs)

# Filter out scenes where VIIRS data could not be found/downloaded
n_before <- nrow(files)
files <- files %>% filter(!is.na(VIIRS_file_nc) & file.exists(VIIRS_file_nc)) ## DEBUG: Added file.exists check
if (nrow(files) < n_before) {
  message(sprintf("Warning: Filtered out %d scenes due to missing VIIRS data.", n_before - nrow(files)))
}

message("--- VIIRS data acquisition complete ---\n")

# -----------------------------------------------------------------------------
# Step 3: Function to Extract and Process Data for a Single Scene
# -----------------------------------------------------------------------------
# This function reads VIIRS and MSAT files, processes them, and returns a
# single dataframe of colocated pixels for one scene.
extract_scene_data <- function(viirs_file_nc, MSAT_file_nc, scene_n) {
  if (is.na(viirs_file_nc) || !file.exists(viirs_file_nc)) {
    warning(paste("Skipping", scene_n, "— VIIRS file missing or NA"))
    return(NULL)
  }

  if (is.na(MSAT_file_nc) || !file.exists(MSAT_file_nc)) {
    warning(paste("Skipping", scene_n, "— MSAT file missing or NA"))
    return(NULL)
  }
  # --- Process VIIRS Data ---
  viirs_nc <- nc_open(viirs_file_nc)
  lon_viirs <- ncvar_get(viirs_nc, "Longitude")
  lat_viirs <- ncvar_get(viirs_nc, "Latitude")
  viirs_data_df <- data.frame(
    lon = as.vector(lon_viirs),
    lat = as.vector(lat_viirs),
    AOD550 = as.vector(ncvar_get(viirs_nc, "Aerosol_Optical_Thickness_550_Land_Best_Estimate")),
    AEROSOL_TYPE = as.vector(ncvar_get(viirs_nc, "Aerosol_Type_Land")),
    AE_QA = as.vector(ncvar_get(viirs_nc, "Angstrom_Exponent_Land_Best_Estimate")),
    AE = as.vector(ncvar_get(viirs_nc, "Angstrom_Exponent_Land"))
  )
  nc_close(viirs_nc)

  viirs_data_df <- viirs_data_df[complete.cases(viirs_data_df), ]
  if (nrow(viirs_data_df) == 0) return(NULL) ## DEBUG: early exit if no valid VIIRS data
  viirs_vect <- vect(viirs_data_df, geom = c("lon", "lat"), crs = "EPSG:4326")
  r_template_viirs <- rast(nrows = dim(lon_viirs)[1], ncols = dim(lon_viirs)[2],
                           xmin = min(lon_viirs), xmax = max(lon_viirs),
                           ymin = min(lat_viirs), ymax = max(lat_viirs),
                           crs = "EPSG:4326")

  viirs_vars_to_rasterize <- c("AOD550", "AEROSOL_TYPE", "AE_QA", "AE")
  viirs_rast <- rasterize(viirs_vect, r_template_viirs,
                          field = viirs_vars_to_rasterize,
                          fun = "mean")
  names(viirs_rast) <- viirs_vars_to_rasterize

  # --- Process MSAT Data ---
  MSAT_nc <- nc_open(MSAT_file_nc)
  o2_main_quality_flag <- ncvar_get(MSAT_nc, "product_o2proxy/main_quality_flag")
  xco2_o2proxy <- ncvar_get(MSAT_nc, "product_o2proxy/xco2_o2proxy_vcdcor")
  albedo_1286 <- ncvar_get(MSAT_nc, "o2dp_fit_diagnostics/retrieved_albedo_1286nm")
  dp <- ncvar_get(MSAT_nc, "o2dp_fit_diagnostics/bias_corrected_delta_pressure")
  airglow_vcd <- ncvar_get(MSAT_nc, "o2dp_fit_diagnostics/airglow_vcd")
  co2_prior <- tryCatch(ncvar_get(MSAT_nc, "apriori_data/co2_profile"), error = function(e) NA)

  # screen for clouds and low quality retrievals
  mask <- o2_main_quality_flag == 1
  dp[mask] <- NA
  xco2_o2proxy[mask] <- NA
  albedo_1286[mask] <- NA # error could be here


  msat_data_df <- data.frame(
    lon = as.vector(ncvar_get(MSAT_nc, "geolocation/longitude")),
    lat = as.vector(ncvar_get(MSAT_nc, "geolocation/latitude")),
    XCO2_O2PROXY = as.vector(xco2_o2proxy),
    ALBEDO_1286 = as.vector(albedo_1286),
    DP = as.vector(dp),
    AIRGLOW_VCD = as.vector(airglow_vcd),
    XCO2_PRIOR = as.vector(co2_prior)
  )
  nc_close(MSAT_nc)

  msat_data_df <- msat_data_df[complete.cases(msat_data_df), ]
  if (nrow(msat_data_df) == 0) return(NULL) ## DEBUG: early exit if no valid MSAT data
  msat_vect <- vect(msat_data_df, geom = c("lon", "lat"), crs = "EPSG:4326")

  msat_vars_to_rasterize <- c("XCO2_O2PROXY", "ALBEDO_1286",
                              "DP", "AIRGLOW_VCD", "XCO2_PRIOR")
  msat_rast <- rasterize(msat_vect, r_template_viirs,
                         field = msat_vars_to_rasterize,
                         fun = "mean")
  names(msat_rast) <- msat_vars_to_rasterize

  # --- Colocate and Combine ---
  ext_intersect <- terra::intersect(ext(viirs_rast), ext(msat_rast))
  if (is.null(ext_intersect)) {
    warning(paste("No spatial overlap for scene:", scene_n))
    return(NULL)
  }

  stack_all <- c(crop(viirs_rast, ext_intersect), crop(msat_rast, ext_intersect))
  df_scene <- as.data.frame(stack_all, xy = TRUE, na.rm = TRUE)

  if (nrow(df_scene) == 0) {
    warning(paste("No valid overlapping pixels found for scene:", scene_n))
    return(NULL)
  }

  df_scene$scene <- scene_n

  return(df_scene)
}

# -----------------------------------------------------------------------------
# Step 4: Loop Through All Files, Extract Data, and Combine
# -----------------------------------------------------------------------------
message("--- Step 4: Starting Data Extraction for all scenes ---")

# Use parallel processing for scene extraction
num_cores <- detectCores() - 1  # Use all but one core
cl <- makeCluster(num_cores)
clusterExport(cl, varlist = c("extract_scene_data", "files"))
clusterEvalQ(cl, {
  library(terra)
  library(ncdf4)
})

all_pixels_list <- parLapply(cl, 1:nrow(files), function(i) {
  cat("Processing scene:", files$scene_n[i], "\n")
  tryCatch({
    extract_scene_data(
      viirs_file_nc = files$VIIRS_file_nc[i],
      MSAT_file_nc = files$MSAT_file_nc[i],
      scene_n = files$scene_n[i]
    )
  }, error = function(e) {
    cat("ERROR processing scene", files$scene_n[i], ":", e$message, "\n")
    return(NULL)  # Return NULL on error to skip the file
  })
})
stopCluster(cl)

# Combine results into a single dataframe
all_pixels_df <- bind_rows(all_pixels_list) ## DEBUG: Switched to bind_rows for safety and consistency.
message("--- Data Extraction Complete ---")


# -----------------------------------------------------------------------------
# Step 5: Post-Process the Combined DataFrame
# -----------------------------------------------------------------------------
## FIX: Added an is.null() check. The error "argument is of length zero" occurs
## when all_pixels_df is NULL, which happens if no data is extracted for any scene.
## This makes the check more robust.
if (is.null(all_pixels_df) || nrow(all_pixels_df) == 0) {
  stop("No data was extracted. Please check file paths, scene overlap, or filters.")
}

# Ensure necessary libraries are loaded in the main session
library(dplyr)
library(MASS)
library(broom)
library(tidyr)
library(caret)

# Compute AOD at 1286nm using the Angstrom power law
all_pixels_df <- all_pixels_df %>%
  mutate(AOD_1286 = AOD550 * (1286 / 550)^(-AE_QA))

# Define aerosol labels for clarity and add as a new column
aerosol_labels <- c(
  `-999` = "No retrieval", `0` = "Dust", `1` = "Smoke",
  `2` = "High altitude smoke", `3` = "Pyrocumulonimbus",
  `4` = "Non-smoke fine mode", `5` = "Mixed", `6` = "Background"
)
all_pixels_df <- all_pixels_df %>%
  mutate(Aerosol_Label = dplyr::recode(as.character(round(AEROSOL_TYPE)), !!!aerosol_labels, .default = "Unknown"))

# -----------------------------------------------------------------------------
# Step 6: Run Regression Analysis for Each Predictor Group (Sequentially)
# -----------------------------------------------------------------------------
unique_groups <- unique(all_pixels_df$Aerosol_Label)

message("\n--- Step 6: Starting Analysis per Predictor Group ---")

# Use lapply for sequential processing of regression analysis
analysis_results_list <- lapply(unique_groups, function(current_group) {
  cat("\nAnalyzing Group:", current_group, "\n")

  group_df <- all_pixels_df %>%
    filter(Aerosol_Label == current_group) %>%
    dplyr::select(XCO2_O2PROXY, AOD_1286, AE, ALBEDO_1286, DP, AIRGLOW_VCD) %>%
    drop_na()

  if (nrow(group_df) < 50) {
    cat("Skipping '", current_group, "' due to insufficient data (", nrow(group_df), " valid pixels).\n", sep = "")
    return(NULL)
  }

  predictors <- group_df %>% dplyr::select(AOD_1286, AE, ALBEDO_1286)
  
  # Remove any predictors that are constant
  non_constant_vars <- predictors[, sapply(predictors, function(x) length(unique(x)) > 1), drop = FALSE]

  # DEBUG: Added a check for near-zero variance predictors. This is a more robust
  # way to remove variables that are almost constant, which can also cause singularity.
  if (ncol(non_constant_vars) > 1) {
    nzv_indices <- nearZeroVar(non_constant_vars)
    if (length(nzv_indices) > 0) {
      cat("   -> For group '", current_group, "', removing near-zero variance predictors: ", paste(names(non_constant_vars)[nzv_indices], collapse=", "), "\n", sep="")
      non_constant_vars <- non_constant_vars[, -nzv_indices, drop = FALSE]
    }
  }

  # Check for multicollinearity among the remaining predictors
  if (ncol(non_constant_vars) > 1) {
      correlation_matrix <- cor(non_constant_vars)
      highly_correlated <- findCorrelation(correlation_matrix, cutoff = 0.90, names = TRUE)
      
      if(length(highly_correlated) > 0) {
          cat("   -> For group '", current_group, "', removing highly correlated predictors: ", paste(highly_correlated, collapse=", "), "\n", sep="")
          non_constant_vars <- non_constant_vars %>%
              dplyr::select(-all_of(highly_correlated))
      }
  }

  if (ncol(non_constant_vars) == 0) {
    cat("Skipping '", current_group, "' due to all predictors being constant or collinear.\n")
    return(NULL)
  }
  
  model_formula <- as.formula(
    paste("XCO2_O2PROXY ~", paste(names(non_constant_vars), collapse = " + "))
  )

  model <- tryCatch({
    rlm(model_formula, data = group_df, maxit = 100)
  }, error = function(e) {
    cat("   -> ERROR: Could not fit model for group '", current_group, "'.\n", sep = "")
    cat("       Reason:", e$message, "\n")
    return(NULL)
  })

  if (is.null(model)) {
    return(NULL)
  }

  pseudo_r_squared <- cor(fitted(model), group_df$XCO2_O2PROXY)^2

  # Combine results into a single dataframe for this group
  tidy_model <- tidy(model)
  tidy_model$group <- current_group
  tidy_model$r_squared <- pseudo_r_squared
  tidy_model$n_pixels <- nrow(group_df)
  
  return(tidy_model)
})

# -----------------------------------------------------------------------------
# Step 7: Consolidate and Save Final Results
# -----------------------------------------------------------------------------
final_summary_df <- bind_rows(analysis_results_list)

if (!is.null(final_summary_df) && nrow(final_summary_df) > 0) {
  final_summary_df <- final_summary_df %>%
    dplyr::select(group, n_pixels, r_squared, term, estimate, std.error, statistic)

  # Create the directory if it doesn't exist
  output_dir <- file.path(PLOT_DIR, "XCO2")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_csv_path <- file.path(output_dir, "DP_predictor_regression_summary.csv")
  write.csv(final_summary_df, output_csv_path, row.names = FALSE)

  message("\nAnalysis summary saved to: ", output_csv_path)
  print(final_summary_df)
} else {
  message("\nNo analysis was performed, likely due to lack of data or filtering.")
}
