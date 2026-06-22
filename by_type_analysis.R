# -----------------------------------------------------------------------------
# Step 0: Load Libraries and Define File Paths
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
library(purrr) # Used for iterating and combining data
library(utils) # For combn function
library(car) # For VIF (Variance Inflation Factor) analysis
# Note: Using custom Moran's I calculation instead of spdep to avoid dependency issues

# ===== INSTRUCTION 1: KNN-BASED MORAN'S I CALCULATION =====
# Replaces full distance matrix approach with k-nearest neighbors.
# Handles large datasets (>2000 pixels) via subsampling; permutation-based p-value.
calculate_morans_i_knn <- function(residuals, coords, k = 8, max_n = 2000) {
  "Calculate Moran's I using k-nearest neighbor weights and permutation test
  
  Args:
    residuals: vector of residuals
    coords: nx2 matrix of (x, y) coordinates
    k: number of nearest neighbors for weight construction (default 8)
    max_n: max observations for computation; if n > max_n, subsample with warning
  
  Returns:
    List with I statistic, expected value, p-value, method, k, n_used
  "
  n_orig <- length(residuals)
  
  # Remove NAs first
  valid_idx <- !is.na(residuals)
  residuals <- residuals[valid_idx]
  coords <- coords[valid_idx, ]
  n_after_na <- length(residuals)
  
  if (n_after_na < 3) {
    return(list(I = NA, expected = NA, p = NA, method = "KNN-permutation", k = k, 
                n_used = n_after_na, note = "too few valid residuals"))
  }
  
  # Subsample if necessary
  subsampled <- FALSE
  if (n_after_na > max_n) {
    set.seed(42)
    idx <- sample(n_after_na, max_n)
    residuals <- residuals[idx]
    coords <- coords[idx, ]
    n <- max_n
    subsampled <- TRUE
    cat(sprintf("  [Moran's I] Subsampled from %d to %d obs for computation\n", n_after_na, max_n))
  } else {
    n <- n_after_na
  }
  
  # Build k-nearest neighbor weight matrix
  # Each observation gets weight 1 for its k nearest neighbors, 0 otherwise
  dists <- as.matrix(dist(coords))
  diag(dists) <- Inf
  W <- matrix(0, n, n)
  for (i in 1:n) {
    nn_idx <- order(dists[i, ])[1:min(k, n - 1)]
    W[i, nn_idx] <- 1
  }
  
  # Row-standardize weights
  row_sums <- rowSums(W)
  row_sums[row_sums == 0] <- 1  # Avoid division by zero
  W <- W / row_sums
  W0 <- sum(W)
  
  # Center residuals
  y <- residuals - mean(residuals)
  
  # Calculate Moran's I
  if (sum(y^2) == 0) {
    return(list(I = NA, expected = NA, p = NA, method = "KNN-permutation", k = k, 
                n_used = n, note = "zero variance in residuals"))
  }
  
  I <- (n / W0) * (sum(y * (W %*% y)) / sum(y^2))
  E_I <- -1 / (n - 1)
  
  # Permutation-based p-value (999 permutations, more robust than analytical)
  set.seed(42 + n)  # n varies by stratum, giving reproducible but distinct sequences
  perm_I <- replicate(999, {
    y_perm <- sample(y)
    (n / W0) * (sum(y_perm * (W %*% y_perm)) / sum(y_perm^2))
  })
  p_value <- mean(abs(perm_I) >= abs(I))
  
  return(list(I = I, expected = E_I, p = p_value, method = "KNN-permutation", 
              k = k, n_used = n, subsampled = subsampled, permutation_seed = 42 + n))
}

# ===== INSTRUCTION 3: CROSS-SCENE CONSISTENCY TEST =====
# Test whether DP ~ predictor slopes are consistent across scenes via F-test
test_scene_consistency <- function(best_model, best_combo, data_scaled, weights = NULL) {
  "Fit interaction model and test slope consistency across scenes
  
  Args:
    best_model: the best linear model
    best_combo: character vector of predictor names
    data_scaled: analysis dataframe with scene factor
    weights: optional observation weights
  
  Returns:
    List with F-stat, p-value, slopes_consistent flag, interpretation note
  "
  if (length(best_combo) == 0) {
    return(list(f_stat = NA, p_value = NA, slopes_consistent = NA, 
                note = "No predictors to test"))
  }
  
  tryCatch({
    # Build interaction formula: predictors * scene
    interaction_terms <- paste(paste(best_combo, collapse = " + "), 
                              "* as.factor(scene)", sep = " ")
    interaction_formula <- as.formula(paste("DP ~", interaction_terms))
    
    interaction_model <- if (!is.null(weights)) {
      lm(interaction_formula, data = data_scaled, weights = weights)
    } else {
      lm(interaction_formula, data = data_scaled)
    }
    
    # F-test: compare models via ANOVA
    anova_result <- anova(best_model, interaction_model)
    f_stat <- anova_result$F[2]
    p_val  <- anova_result$`Pr(>F)`[2]
    
    slopes_consistent <- !is.na(p_val) && p_val >= 0.05
    
    note <- ifelse(slopes_consistent,
                   "Slopes consistent across scenes (interaction NS)",
                   "WARNING: Slopes differ significantly across scenes")
    
    return(list(f_stat = f_stat, p_value = p_val, slopes_consistent = slopes_consistent, 
                note = note))
  }, error = function(e) {
    return(list(f_stat = NA, p_value = NA, slopes_consistent = NA, 
                note = paste("interaction test failed:", e$message)))
  })
}

# ===== INSTRUCTION 2: SPATIAL BLOCK CROSS-VALIDATION =====
# Out-of-sample validation accounting for spatial autocorrelation
spatial_block_cv <- function(model_formula, data_scaled, coords, weights = NULL, n_blocks = 5) {
  "Spatial block cross-validation: fit on train blocks, predict on held-out block
  
  Args:
    model_formula: model formula
    data_scaled: analysis dataframe
    coords: nx2 matrix of coordinates
    weights: optional observation weights (must be same length as nrow(data_scaled))
    n_blocks: number of spatial blocks (default 5)
  
  Returns:
    List with CV RMSE (mean, sd), CV R², number of blocks used
  "
  if (nrow(data_scaled) < 30 || nrow(coords) != nrow(data_scaled)) {
    return(list(cv_rmse_mean = NA, cv_rmse_sd = NA, cv_r2_mean = NA, 
                n_blocks_used = 0, note = "insufficient data for CV"))
  }
  
  # Create spatial blocks using quantile-based grid
  x_breaks <- quantile(coords[, 1], probs = seq(0, 1, length.out = n_blocks + 1), na.rm = TRUE)
  y_breaks <- quantile(coords[, 2], probs = seq(0, 1, length.out = n_blocks + 1), na.rm = TRUE)
  
  x_breaks <- unique(x_breaks)
  y_breaks <- unique(y_breaks)
  
  block_x <- cut(coords[, 1], breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  block_y <- cut(coords[, 2], breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  block_id <- interaction(block_x, block_y)
  
  unique_blocks <- unique(block_id)
  cv_rmse <- numeric(length(unique_blocks))
  cv_r2 <- numeric(length(unique_blocks))
  n_valid <- 0
  
  for (b in seq_along(unique_blocks)) {
    test_idx <- which(block_id == unique_blocks[b])
    train_idx <- which(block_id != unique_blocks[b])
    
    if (length(train_idx) < 20 || length(test_idx) < 5) next
    
    train_data <- data_scaled[train_idx, ]
    test_data  <- data_scaled[test_idx, ]
    train_weights <- if (!is.null(weights)) weights[train_idx] else NULL
    
    fold_model <- tryCatch({
      if (!is.null(train_weights)) {
        lm(model_formula, data = train_data, weights = train_weights)
      } else {
        lm(model_formula, data = train_data)
      }
    }, error = function(e) NULL)
    
    if (is.null(fold_model)) next
    
    preds <- tryCatch(
      predict(fold_model, newdata = test_data),
      error = function(e) NULL
    )
    if (is.null(preds)) next
    
    actuals <- test_data$DP
    valid <- !is.na(preds) & !is.na(actuals)
    if (sum(valid) < 5) next
    
    cv_rmse[b] <- sqrt(mean((preds[valid] - actuals[valid])^2))
    ss_res <- sum((actuals[valid] - preds[valid])^2)
    ss_tot <- sum((actuals[valid] - mean(actuals[valid]))^2)
    cv_r2[b] <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA)
    n_valid <- n_valid + 1
  }
  
  if (n_valid == 0) {
    return(list(cv_rmse_mean = NA, cv_rmse_sd = NA, cv_r2_mean = NA, 
                n_blocks_used = 0, note = "all folds failed"))
  }
  
  return(list(
    cv_rmse_mean = mean(cv_rmse[1:n_valid], na.rm = TRUE),
    cv_rmse_sd   = sd(cv_rmse[1:n_valid], na.rm = TRUE),
    cv_r2_mean   = mean(cv_r2[1:n_valid], na.rm = TRUE),
    n_blocks_used = n_valid
  ))
}

# ===== INSTRUCTION 5: NONLINEARITY TESTING =====
# Test whether adding quadratic terms improves model fit
test_nonlinearity <- function(best_model, best_combo, data_scaled, weights = NULL, single_scene_flag = FALSE) {
  "Fit quadratic model and compare with linear via AIC
  
  Args:
    best_model: the best linear model
    best_combo: character vector of predictor names (must be continuous)
    data_scaled: analysis dataframe
    weights: optional observation weights
    single_scene_flag: if TRUE, omit as.factor(scene) from quadratic formula
  
  Returns:
    List with delta_aic, aic values, nonlinear_preferred flag
  "
  if (length(best_combo) == 0) {
    return(list(delta_aic = NA, aic_linear = NA, aic_quadratic = NA, 
                nonlinear_preferred = FALSE))
  }
  
  tryCatch({
    # Build linear and quadratic formulas
    linear_terms <- paste(best_combo, collapse = " + ")
    quad_terms <- paste(paste0("I(", best_combo, "^2)"), collapse = " + ")
    
    # Conditionally include scene fixed effect (only if multi-scene stratum)
    if (single_scene_flag) {
      quad_formula <- as.formula(paste("DP ~", linear_terms, "+", quad_terms))
    } else {
      quad_formula <- as.formula(paste("DP ~", linear_terms, "+", quad_terms, "+ as.factor(scene)"))
    }
    
    quad_model <- if (!is.null(weights)) {
      lm(quad_formula, data = data_scaled, weights = weights)
    } else {
      lm(quad_formula, data = data_scaled)
    }
    
    aic_lin <- AIC(best_model)
    aic_quad <- AIC(quad_model)
    delta_aic <- aic_lin - aic_quad
    
    # Positive delta_aic means quadratic model is preferred (lower AIC)
    nonlinear_preferred <- delta_aic > 2  # Threshold: AIC reduction > 2
    
    return(list(
      delta_aic = delta_aic,
      aic_linear = aic_lin,
      aic_quadratic = aic_quad,
      nonlinear_preferred = nonlinear_preferred
    ))
  }, error = function(e) {
    return(list(delta_aic = NA, aic_linear = NA, aic_quadratic = NA, 
                nonlinear_preferred = FALSE, note = paste("test failed:", e$message)))
  })
}

# Define file paths and scene identifiers
# NOTE: Ensure you have at least one scene uncommented for the script to run.
files <- data.frame(
  scene_name = c(
    "MAURITANIA, 1",
    "BANGLADESH, DHAKA",
    "UKRAINE, KHARKIV",
    "USA, CALIFORNIA, TULARE",
    "LIBYA, 4",
    "SPAIN, MADRID",
    # "MAURITANIA, 2",
    # "MAURITANIA, 2.1",
    # "MAURITANIA, 2.2",
    "CHINA, HEFEI"
  ),
  viirs_file_nc = c(
    "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/diff_bands/AERDB_L2_VIIRS_SNPP_2-20250707_203516/AERDB_L2_VIIRS_SNPP.A2024329.1254.002.2024330011126.nc",
    "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/diff_bands/AERDB_L2_VIIRS_SNPP_2-20250709_161443/AERDB_L2_VIIRS_SNPP.A2025078.0654.002.2025078193026.nc",
    "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/diff_bands/AERDB_L2_VIIRS_SNPP_2-20250715_012523/AERDB_L2_VIIRS_SNPP.A2025158.0848.002.2025161020511.nc",
    "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/diff_bands/AERDB_L2_VIIRS_SNPP_2-20250708_200258/AERDB_L2_VIIRS_SNPP.A2024256.2054.002.2024257091442.nc",
    "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/diff_bands/AERDB_L2_VIIRS_SNPP_2-20250723_205854/AERDB_L2_VIIRS_SNPP.A2024324.1106.002.2024324232733.nc",
    "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/diff_bands/AERDB_L2_VIIRS_SNPP_2-20250725_144124/AERDB_L2_VIIRS_SNPP.A2025116.1330.002.2025117015659.nc",
    # "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/diff_bands/AERDB_L2_VIIRS_SNPP_2-20250725_145734/AERDB_L2_VIIRS_SNPP.A2025152.1354.002.2025153020439.nc",
    # "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/diff_bands/AERDB_L2_VIIRS_SNPP_2-20250725_150428/AERDB_L2_VIIRS_SNPP.A2025130.1406.002.2025131022022.nc",
    # "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/diff_bands/AERDB_L2_VIIRS_SNPP_2-20250725_201702/AERDB_L2_VIIRS_SNPP.A2024298.1418.002.2024299023741.nc",
    "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/diff_bands/AERDB_L2_VIIRS_SNPP_2-20250725_202930/AERDB_L2_VIIRS_SNPP.A2025079.0454.002.2025079174029.nc"
  ),
  MSAT_file_nc = c(
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t26/20241124/c049501A0/p748/MSAT_L2_c049501A0_p748_v03002000_20241124T143538Z_143610Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t380/20250319/c08F617C0/p1524/MSAT_L2_c08F617C0_p1524_v03005001_20250319T080703Z_080735Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t122/20250607/c0D0707A0/p2610/MSAT_L2_c0D0707A0_p2610_v03005001_20250607T120614Z_120646Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/t159/2024/09/12/014E09F0/202409171319_po-2103/level2-results/MSAT_L2_post_20240912T221258_20240912T221330_post.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t25/20241119/c04350190/p898/MSAT_L2_c04350190_p898_v03002000_20241119T122526Z_122558Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t381/20250426/c0A4217D0/p1618/MSAT_L2_c0A4217D0_p1618_v03005001_20250426T144316Z_144348Z.nc",
    # "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t27/20250601/c0C7F01B0/p2309/MSAT_L2_c0C7F01B0_p2309_v03005001_20250601T145352Z_145424Z.nc",
    # "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t27/20250510/c0B3601B0/p1734/MSAT_L2_c0B3601B0_p1734_v03005001_20250510T143014Z_143046Z.nc",
    # "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t27/20241024/c032201B0/p668/MSAT_L2_c032201B0_p668_v03002000_20241024T142657Z_142729Z.nc",
    "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t361/20250320/c09041690/p1403/MSAT_L2_c09041690_p1403_v03004000_20250320T063412Z_063444Z.nc"
  ),
  scene_n = c(
    "t26p748",
    "t380p1524",
    "t122p2610",
    "t159",
    "t25p898",
    "t381p1618",
    # "t27p2309",
    # "t27p1734",
    # "t27p668",
    "t361p1403"
  )
)

# Directory to save plots and results
plot_dir <- "by_aerosol_data"
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# -----------------------------------------------------------------------------
# Step 1: Function to Extract and Process Data for a Single Scene
# -----------------------------------------------------------------------------
# This function reads VIIRS and MSAT files, processes them, and returns a
# single dataframe of colocated pixels for one scene.
extract_scene_data <- function(viirs_file_nc, MSAT_file_nc, scene_n) {
  
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
  viirs_vect <- vect(viirs_data_df, geom = c("lon", "lat"), crs = "EPSG:4326")
  r_template_viirs <- rast(nrows = dim(lon_viirs)[1], ncols = dim(lon_viirs)[2],
                           xmin = min(lon_viirs), xmax = max(lon_viirs),
                           ymin = min(lat_viirs), ymax = max(lat_viirs),
                           crs = "EPSG:4326")
  viirs_rast <- rasterize(viirs_vect, r_template_viirs,
                          field = c("AOD550", "AEROSOL_TYPE", "AE_QA", "AE"),
                          fun = "mean")

  # --- Process MSAT Data ---
  MSAT_nc <- nc_open(MSAT_file_nc)
  o2_main_quality_flag <- ncvar_get(MSAT_nc, "product_o2proxy/main_quality_flag")
  
  # Extract core MSAT variables (always available)
  dp_raw <- ncvar_get(MSAT_nc, "o2dp_fit_diagnostics/bias_corrected_delta_pressure")
  albedo_1286_raw <- ncvar_get(MSAT_nc, "o2dp_fit_diagnostics/retrieved_albedo_1286nm")
  
  # Extract SZA (with robust error handling)
  sza_raw <- tryCatch(
    ncvar_get(MSAT_nc, "geolocation/sza"),
    error = function(e) {
      cat("  [Warning] SZA not found, using NA\n")
      rep(NA_real_, length(dp_raw))
    }
  )
  
  # Try to extract XCO2 from o2 proxy (bias-corrected version preferred)
  xco2_raw <- tryCatch(
    ncvar_get(MSAT_nc, "product_o2proxy/xco2_o2proxy_vcdcor"),
    error = function(e) {
      tryCatch(
        ncvar_get(MSAT_nc, "product_o2proxy/xco2_o2proxy_raw"),
        error = function(e2) {
          tryCatch(
            ncvar_get(MSAT_nc, "product_o2proxy/xco2"),
            error = function(e3) {
              cat("  [Warning] XCO2 not found, using NA\n")
              rep(NA_real_, length(dp_raw))
            }
          )
        }
      )
    }
  )
  
  # Try to extract second albedo (try multiple wavelengths)
  albedo_secondary_raw <- tryCatch(
    ncvar_get(MSAT_nc, "o2dp_fit_diagnostics/retrieved_albedo_1606nm"),
    error = function(e) {
      tryCatch(
        ncvar_get(MSAT_nc, "o2dp_fit_diagnostics/retrieved_albedo_1294nm"),
        error = function(e2) {
          cat("  [Warning] Secondary albedo not found, using NA\n")
          rep(NA_real_, length(dp_raw))
        }
      )
    }
  )
  
  # Apply cloud filtering
  dp_raw[o2_main_quality_flag == 1] <- NA
  albedo_1286_raw[o2_main_quality_flag == 1] <- NA
  xco2_raw[o2_main_quality_flag == 1] <- NA
  sza_raw[o2_main_quality_flag == 1] <- NA
  albedo_secondary_raw[o2_main_quality_flag == 1] <- NA

  msat_data_df <- data.frame(
    lon = as.vector(ncvar_get(MSAT_nc, "geolocation/longitude")),
    lat = as.vector(ncvar_get(MSAT_nc, "geolocation/latitude")),
    DP = as.vector(dp_raw),
    ALBEDO_1286 = as.vector(albedo_1286_raw),
    SZA = as.vector(sza_raw),
    XCO2 = as.vector(xco2_raw)
  )
  nc_close(MSAT_nc)
  
  # Filter for complete core data (lon, lat, DP, ALBEDO_1286)
  # SZA and XCO2 are optional; if mostly NA, we'll skip them
  complete_core <- complete.cases(msat_data_df[, c("lon", "lat", "DP", "ALBEDO_1286")])
  msat_data_df <- msat_data_df[complete_core, ]
  
  if (nrow(msat_data_df) == 0) {
    warning(paste("No complete core MSAT data for scene:", scene_n))
    return(NULL)
  }
  
  # Build raster with available variables
  available_fields <- c("DP", "ALBEDO_1286")
  
  # Only include SZA if it has >30% non-NA values
  if (sum(!is.na(msat_data_df$SZA)) > 0.3 * nrow(msat_data_df)) {
    available_fields <- c(available_fields, "SZA")
  } else {
    msat_data_df$SZA <- NULL
  }
  
  # Only include XCO2 if it has >30% non-NA values
  if (!is.null(msat_data_df$XCO2) && sum(!is.na(msat_data_df$XCO2)) > 0.3 * nrow(msat_data_df)) {
    available_fields <- c(available_fields, "XCO2")
  } else {
    msat_data_df$XCO2 <- NULL
  }
  
  msat_vect <- vect(msat_data_df, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  # ===== INSTRUCTION 7: Add pixel-count weighting for colocation uncertainty =====
  # Track how many MSAT pixels (~100m) contribute to each VIIRS cell (~750m)
  msat_rast_vals  <- rasterize(msat_vect, r_template_viirs,
                                field = available_fields,
                                fun = "mean")
  msat_rast_count <- rasterize(msat_vect, r_template_viirs,
                                field = available_fields[1],
                                fun = "count")
  names(msat_rast_count) <- "MSAT_pixel_count"
  
  # Stack count with values
  msat_rast <- c(msat_rast_vals, msat_rast_count)
  
  # Add count to available fields so it gets included in output
  available_fields_with_count <- c(available_fields, "MSAT_pixel_count")

  # --- Colocate and Combine ---
  # FIX: Explicitly use terra::intersect to avoid masking by dplyr::intersect
  ext_intersect <- terra::intersect(ext(viirs_rast), ext(msat_rast))
  if (is.null(ext_intersect)) {
    warning(paste("No spatial overlap for scene:", scene_n))
    return(NULL) # Return NULL if there's no overlap
  }
  
  stack_all <- c(crop(viirs_rast, ext_intersect), crop(msat_rast, ext_intersect))
  
  # Build names dynamically based on which variables were included
  viirs_names <- c("AOD550", "AEROSOL_TYPE", "AE_QA", "AE")
  msat_names <- available_fields_with_count  # Now includes MSAT_pixel_count
  all_names <- c(viirs_names, msat_names)
  
  names(stack_all) <- all_names
  
  df_scene <- as.data.frame(stack_all, xy = TRUE, na.rm = TRUE)
  
  # Add scene identifier
  df_scene$scene <- scene_n
  
  return(df_scene)
}

# -----------------------------------------------------------------------------
# Step 2: Loop Through All Files, Extract Data, and Combine into One DataFrame
# -----------------------------------------------------------------------------
# Use purrr::map_dfr to iterate through each row of the 'files' dataframe,
# call the extraction function, and row-bind the results into a single df.
cat("--- Starting Data Extraction ---\n")
all_pixels_df <- map_dfr(1:nrow(files), ~{
  cat("Processing scene:", files$scene_n[.x], "\n")
  # Use a try-catch block to handle potential errors in file processing
  tryCatch({
    extract_scene_data(
      viirs_file_nc = files$viirs_file_nc[.x],
      MSAT_file_nc = files$MSAT_file_nc[.x],
      scene_n = files$scene_n[.x]
    )
  }, error = function(e) {
    cat("ERROR processing scene", files$scene_n[.x], ":", e$message, "\n")
    return(NULL) # Return NULL on error to skip the file
  })
})
cat("--- Data Extraction Complete ---\n")

# -----------------------------------------------------------------------------
# Step 3: Post-Process the Combined DataFrame
# -----------------------------------------------------------------------------
# Check if the dataframe is empty
if (nrow(all_pixels_df) == 0) {
  stop("No data was extracted. Please check file paths and scene overlap.")
}

# Compute AOD at 1286nm using the Angstrom power law
all_pixels_df <- all_pixels_df %>%
  mutate(
    AOD_1286 = AOD550 * (1286 / 550)^(-AE_QA),
    COS_SZA = cos(SZA * pi / 180)  # Convert SZA to radians and compute cosine
  )

# Define aerosol labels for clarity
aerosol_labels <- c(
  `-999` = "No retrieval", `0` = "Dust", `1` = "Smoke",
  `2` = "High altitude smoke", `3` = "Pyrocumulonimbus",
  `4` = "Non-smoke fine mode", `5` = "Mixed", `6` = "Background"
)

# Add a human-readable aerosol label column
all_pixels_df <- all_pixels_df %>%
  mutate(
    # Round AEROSOL_TYPE in case of rasterization averaging, then convert to character for recoding
    Aerosol_Label = dplyr::recode(as.character(round(AEROSOL_TYPE)), !!!aerosol_labels)
  )

# Function to categorize aerosols by Angstrom Exponent
categorize_by_ae <- function(ae_values) {
  cut(ae_values, 
      breaks = c(-Inf, 0.5, 1, 2, Inf),
      labels = c("Very Large", "Large", "Medium", "Small"),
      right = FALSE)
}

# Add AE categories
all_pixels_df <- all_pixels_df %>%
  mutate(
    AE_Category = categorize_by_ae(AE)
  )

# Display the structure and summary of the final combined dataframe
message("Structure of the final combined dataframe (all_pixels_df):")
str(all_pixels_df)
message("\nSummary of key variables:")
summary_vars <- c("DP", "AOD_1286", "ALBEDO_1286", "AE", "COS_SZA", "Aerosol_Label", "AE_Category")
# Add optional variables if they exist
if ("SZA" %in% names(all_pixels_df)) summary_vars <- c(summary_vars, "SZA")
if ("XCO2" %in% names(all_pixels_df)) summary_vars <- c(summary_vars, "XCO2")
summary(all_pixels_df[, summary_vars[summary_vars %in% names(all_pixels_df)]])

# ===== INSTRUCTION 4: GLOBAL OUTLIER EXCLUSION (Physics-Justified) =====
# Remove extreme DP values that likely indicate retrieval failures or cloud contamination
# Use 3 × IQR from scene-specific median (robust to outliers, comparable across scenes)
cat("\n--- Removing DP Outliers (3 × IQR threshold) ---\n")
n_pixels_before <- nrow(all_pixels_df)

all_pixels_df <- all_pixels_df %>%
  group_by(scene) %>%
  mutate(
    DP_IQR = IQR(DP, na.rm = TRUE),
    DP_median = median(DP, na.rm = TRUE),
    DP_outlier = abs(DP - DP_median) > 3 * DP_IQR
  ) %>%
  ungroup() %>%
  filter(!DP_outlier) %>%
  dplyr::select(-DP_IQR, -DP_median, -DP_outlier)

n_outliers_removed <- n_pixels_before - nrow(all_pixels_df)
cat(sprintf("Removed %d outlier pixels (%.2f%% of original %d pixels)\n", 
            n_outliers_removed, 100 * n_outliers_removed / n_pixels_before, n_pixels_before))
cat("Rationale: Extreme DP values (> 3 × IQR from scene median) likely indicate\n")
cat("           retrieval failures or cloud contamination, not genuine aerosol effects.\n")

# -----------------------------------------------------------------------------
# Step 4: Run Regression Analysis for Each Aerosol Type
# -----------------------------------------------------------------------------
# Get the unique aerosol types present in the dataset
unique_aerosol_labels <- unique(all_pixels_df$Aerosol_Label)

# Prepare a list to store the results of each analysis
analysis_results <- list()
model_comparison_results <- list()

cat("\n--- Starting Rigorous Analysis per Aerosol Type with AE Stratification ---\n")

for (aerosol in unique_aerosol_labels) {
  
  cat("\n========================================\n")
  cat("Analyzing Aerosol Type:", aerosol, "\n")
  cat("========================================\n")
  
  # Filter data for the current aerosol type
  aerosol_base_df <- all_pixels_df %>%
    filter(Aerosol_Label == aerosol)
  
  if (nrow(aerosol_base_df) < 50) {
    cat("Skipping '", aerosol, "' due to insufficient data (", nrow(aerosol_base_df), " valid pixels).\n", sep = "")
    next
  }
  
  # Get unique AE categories for this aerosol type
  ae_categories <- unique(aerosol_base_df$AE_Category)
  
  for (ae_cat in ae_categories) {
    
    cat("\n--- AE Category:", as.character(ae_cat), "---\n")
    
    # Filter data for current aerosol type and AE category
    aerosol_ae_df <- aerosol_base_df %>%
      filter(AE_Category == ae_cat)
    
    if (nrow(aerosol_ae_df) < 30) {
      cat("Skipping due to insufficient data (", nrow(aerosol_ae_df), " pixels).\n", sep = "")
      next
    }
    
    # NOTE: Background aerosol filtering removed (per INSTRUCTION 4)
    # Global outlier exclusion (3*IQR) already applied in Step 3 above
    # Do NOT drop_na() yet — defer complete-case filtering to after predictors are finalized
    
    # Select only numeric columns for analysis
    numeric_cols <- sapply(aerosol_ae_df, is.numeric)
    analysis_df <- aerosol_ae_df[, numeric_cols]
    
    # Remove target variable and non-predictor columns
    # ===== INSTRUCTION 5: Exclude AE from predictor pool =====
    # Rationale: AE is used to define AE_Category bins. Within a bin, AE has
    # artificially restricted range and its coefficient reflects within-bin
    # residual variation, not the full physical relationship. Exclude to avoid circularity.
    # ===== INSTRUCTION 7: Exclude MSAT_pixel_count from predictors (use as weights) =====
    cols_to_exclude <- c("x", "y", "AEROSOL_TYPE", "AE_QA", "AE", "MSAT_pixel_count")
    predictors_initial <- setdiff(names(analysis_df), c("DP", cols_to_exclude))
    
    # Remove columns with zero variance
    analysis_df <- analysis_df %>% dplyr::select(DP, all_of(predictors_initial))
    
    # Remove columns with low variance: require (1) min 5 unique values AND (2) CV > 0.5% OR IQR > 0
    # This prevents near-constant predictors that cause numerical instability
    variance_check <- sapply(analysis_df[, -1], function(x) {
      n_unique <- length(unique(x))
      if (n_unique < 5) return(FALSE)  # At least 5 unique values
      mean_x <- mean(x, na.rm = TRUE)
      sd_x <- sd(x, na.rm = TRUE)
      cv <- if (mean_x != 0) abs(sd_x / mean_x) else sd_x  # Avoid division by zero
      iqr_x <- IQR(x, na.rm = TRUE)
      # Pass if CV > 0.5% OR IQR > 0
      return(cv > 0.005 || iqr_x > 0)
    })
    predictors <- names(analysis_df)[-1][variance_check]
    
    if (length(predictors) == 0) {
      cat("Skipping: no valid predictors with variation.\n")
      next
    }
    
    # ===== INSTRUCTION 3: Complete-case filtering (deferred to after predictors finalized) =====
    # Filter to complete cases on only the relevant columns needed for analysis
    cols_needed <- c("DP", "x", "y", "scene", predictors)
    if ("MSAT_pixel_count" %in% names(aerosol_ae_df)) {
      cols_needed <- c(cols_needed, "MSAT_pixel_count")
    }
    # Keep only cols that actually exist (defensive check)
    cols_needed <- cols_needed[cols_needed %in% names(aerosol_ae_df)]
    
    n_before_cc <- nrow(aerosol_ae_df)
    aerosol_ae_df <- aerosol_ae_df[complete.cases(aerosol_ae_df[, cols_needed]), ]
    n_after_cc <- nrow(aerosol_ae_df)
    
    if (n_before_cc != n_after_cc) {
      cat(sprintf("  Complete-case filter: removed %d rows with NA in analysis columns.\n",
                  n_before_cc - n_after_cc))
    }
    
    if (nrow(aerosol_ae_df) < 30) {
      cat("Skipping: insufficient data after complete-case filter (", 
          nrow(aerosol_ae_df), " pixels).\n")
      next
    }
    
    # ===== INSTRUCTION 1 & 3: Define coordinates early — now guaranteed aligned =====
    # coords_mat is now derived from the cleaned aerosol_ae_df after complete-case filtering
    coords_mat <- cbind(aerosol_ae_df$x, aerosol_ae_df$y)
    
    # ===== INSTRUCTION 6: Extract pixel-count weights BEFORE analysis_df_scaled is built =====
    # This must happen before analysis_df_scaled construction and before analysis_df excludes MSAT_pixel_count
    if ("MSAT_pixel_count" %in% names(aerosol_ae_df)) {
      obs_weights <- aerosol_ae_df$MSAT_pixel_count
      obs_weights[is.na(obs_weights) | obs_weights <= 0] <- 1
      cat(sprintf("  Pixel-count weights: mean=%.1f, sd=%.1f, range=[%d, %d]\n",
                  mean(obs_weights), sd(obs_weights),
                  as.integer(min(obs_weights)), as.integer(max(obs_weights))))
    } else {
      obs_weights <- rep(1, nrow(aerosol_ae_df))
      cat("  No MSAT_pixel_count found — using uniform weights.\n")
    }
    
    # ===== INSTRUCTION 6: Handle single-scene strata explicitly =====
    n_scenes_in_stratum <- length(unique(aerosol_ae_df$scene))
    single_scene_flag <- FALSE
    
    if (n_scenes_in_stratum == 1) {
      cat(sprintf(
        "  ⚠ SINGLE-SCENE STRATUM: Only scene '%s' for %s / %s.\n",
        unique(aerosol_ae_df$scene), aerosol, as.character(ae_cat)
      ))
      cat("    Scene fixed effect not identifiable. Cross-scene consistency test cannot run.\n")
      cat("    Results are scene-specific and cannot generalize. Fitting without scene effect.\n")
      single_scene_flag <- TRUE
    }
    
    cat("Available predictors:", paste(predictors, collapse = ", "), "\n")
    cat("Number of observations:", nrow(aerosol_ae_df), "\n")
    cat("VIF Threshold: 5.0 (models with any predictor having VIF ≥ 5 are excluded)\n")
    if (!single_scene_flag) {
      cat("Scene included as fixed effect in all models.\n")
    } else {
      cat("[Single scene] Scene fixed effect omitted.\n")
    }
    
    # Function to check for structural collinearity
    # AOD_1286 is computed from AOD550 and AE, so exclude combo: AOD_1286 + AOD550 + AE
    is_structurally_collinear <- function(combo) {
      has_aod_1286 <- "AOD_1286" %in% combo
      has_aod550 <- "AOD550" %in% combo
      has_ae <- "AE" %in% combo
      # Exclude if all three are present
      return(has_aod_1286 && has_aod550 && has_ae)
    }
    
    # ===== INSTRUCTION 4: Fixed VIF computation for mixed predictor types =====
    check_vif <- function(formula_str, data, vif_threshold = 5, weights = NULL) {
      temp_model <- tryCatch(
        if (!is.null(weights)) {
          lm(as.formula(formula_str), data = data, weights = weights)
        } else {
          lm(as.formula(formula_str), data = data)
        },
        error = function(e) return(NULL)
      )
      if (is.null(temp_model)) return(list(pass = FALSE, vif_values = NA, max_vif = NA))
      
      vif_raw <- tryCatch(
        vif(temp_model),
        error = function(e) return(NULL)
      )
      if (is.null(vif_raw)) return(list(pass = TRUE, vif_values = NA, max_vif = NA))
      
      # If vif() returns a matrix (GVIF with factors), use GVIF^(1/(2*Df)) which is
      # comparable to sqrt(VIF) for continuous predictors
      if (is.matrix(vif_raw)) {
        # Column 3 is already GVIF^(1/(2*Df)); square it for VIF-comparable statistic
        adjusted_vif <- vif_raw[, 3]^2
      } else {
        adjusted_vif <- vif_raw
      }
      
      # Exclude scene factor from threshold check — it's a fixed effect, not a predictor
      scene_terms <- grepl("scene", names(adjusted_vif), ignore.case = TRUE)
      predictor_vif <- adjusted_vif[!scene_terms]
      
      if (length(predictor_vif) == 0) {
        return(list(pass = TRUE, vif_values = adjusted_vif, max_vif = NA))
      }
      
      max_vif <- max(predictor_vif, na.rm = TRUE)
      pass_vif <- max_vif < vif_threshold
      
      return(list(pass = pass_vif, vif_values = adjusted_vif, max_vif = max_vif))
    }
    
    # ===== INSTRUCTION 2: Build analysis_df_scaled AFTER complete-case filter (not before) =====
    # Rebuild analysis_df from the NOW-CLEANED aerosol_ae_df to ensure row alignment
    numeric_cols <- sapply(aerosol_ae_df[, names(aerosol_ae_df) != "MSAT_pixel_count"], is.numeric)
    analysis_df <- aerosol_ae_df[, c("DP", names(aerosol_ae_df[, numeric_cols])[names(aerosol_ae_df[, numeric_cols]) %in% predictors])]
    
    # Standardize predictors using 3 standard deviations
    analysis_df_scaled <- analysis_df %>%
      mutate(across(all_of(predictors), 
                    function(x) {
                      mean_x <- mean(x, na.rm = TRUE)
                      sd_x <- sd(x, na.rm = TRUE)
                      if (sd_x == 0) return(rep(0, length(x)))
                      return((x - mean_x) / (3 * sd_x))
                    }))
    
    # Add scene as a factor for fixed effects modeling
    analysis_df_scaled$scene <- as.factor(aerosol_ae_df$scene)
    
    # Verify alignment between obs_weights and analysis_df_scaled
    stopifnot(
      "obs_weights length must match analysis_df_scaled rows" =
        length(obs_weights) == nrow(analysis_df_scaled)
    )
    
    cat("Predictors standardized (divided by 3*SD).\n")
    
    # Generate all permutations of predictors (from 1 to max length)
    max_predictors <- min(length(predictors), 5)  # Limit to 5 to avoid too many models
    
    all_models <- list()
    all_aic_results <- data.frame()
    
    # Fit models for each predictor combination
    for (k in 1:max_predictors) {
      model_combos <- combn(predictors, k, simplify = FALSE)
      
      for (combo in model_combos) {
        # Skip structurally collinear combinations
        if (is_structurally_collinear(combo)) {
          next
        }
        
        # Add scene as fixed effect to all models (unless single-scene stratum)
        if (single_scene_flag) {
          formula_str <- paste("DP ~", paste(combo, collapse = " + "))
        } else {
          formula_str <- paste("DP ~", paste(combo, collapse = " + "), "+ as.factor(scene)")
        }
        model_formula <- as.formula(formula_str)
        
        # Check VIF for collinearity (threshold = 5)
        vif_check <- check_vif(formula_str, data = analysis_df_scaled, vif_threshold = 5, weights = obs_weights)
        
        if (!vif_check$pass) {
          # Skip models with high collinearity
          next
        }
        
        # Fit linear model with AIC-minimizing approach
        tryCatch({
          model <- lm(model_formula, data = analysis_df_scaled, weights = obs_weights)
          model_aic <- AIC(model)
          model_bic <- BIC(model)
          model_r2 <- summary(model)$r.squared
          model_adj_r2 <- summary(model)$adj.r.squared
          
          model_name <- paste(combo, collapse = "+")
          max_vif_val <- ifelse(!is.na(vif_check$max_vif), vif_check$max_vif, NA)
          
          all_aic_results <- rbind(all_aic_results, data.frame(
            aerosol_type = aerosol,
            ae_category = as.character(ae_cat),
            model_spec = model_name,
            n_predictors = length(combo),
            aic = model_aic,
            bic = model_bic,
            r_squared = model_r2,
            adj_r_squared = model_adj_r2,
            n_obs = nrow(analysis_df_scaled),
            max_vif = max_vif_val,
            weighted_model = !all(obs_weights == 1),
            stringsAsFactors = FALSE
          ))
          
          all_models[[model_name]] <- model
          
        }, error = function(e) {
          cat("Error fitting model:", model_name, "-", e$message, "\n")
        })
      }
    }
    
    if (nrow(all_aic_results) > 0) {
      # Sort by AIC
      all_aic_results <- all_aic_results %>% 
        arrange(aic) %>%
        mutate(aic_rank = row_number())
      
      message("\nTop 10 models by AIC:")
      print(head(all_aic_results, 10))
      
      # Store best model results
      best_model_idx <- which.min(all_aic_results$aic)
      best_model_spec <- all_aic_results$model_spec[best_model_idx]
      best_model <- all_models[[best_model_spec]]
      
      # Extract predictor names for diagnostic tests
      best_combo <- unlist(strsplit(best_model_spec, "\\+"))
      best_combo <- trimws(best_combo)
      
      # ===== INSTRUCTION 1: STEP A — Test nonlinearity FIRST, before diagnostics =====
      cat("\n📈 NONLINEARITY TESTING (Quadratic Terms)\n")
      nonlin_result <- test_nonlinearity(best_model, best_combo, analysis_df_scaled, weights = obs_weights, single_scene_flag = single_scene_flag)
      
      cat(sprintf("  Linear AIC: %.2f\n", nonlin_result$aic_linear))
      cat(sprintf("  Quadratic AIC: %.2f\n", nonlin_result$aic_quadratic))
      cat(sprintf("  Delta AIC (lin - quad): %.2f\n", nonlin_result$delta_aic))
      
      model_is_quadratic <- FALSE
      
      if (!is.na(nonlin_result$nonlinear_preferred) && nonlin_result$nonlinear_preferred) {
        cat("  → Refitting with quadratic terms (delta AIC > 2).\n")
        
        tryCatch({
          linear_terms <- paste(best_combo, collapse = " + ")
          quad_terms   <- paste(paste0("I(", best_combo, "^2)"), collapse = " + ")
          
          if (single_scene_flag) {
            quad_formula <- as.formula(paste("DP ~", linear_terms, "+", quad_terms))
          } else {
            quad_formula <- as.formula(paste("DP ~", linear_terms, "+", quad_terms, "+ as.factor(scene)"))
          }
          
          quad_model_final <- lm(quad_formula, data = analysis_df_scaled, weights = obs_weights)
          
          best_model      <- quad_model_final
          best_model_spec <- paste0(best_model_spec, "+QUADRATIC")
          model_is_quadratic <- TRUE
          cat(sprintf("  ✓ Quadratic model active. AIC: %.2f | R²: %.4f\n",
                      AIC(best_model), summary(best_model)$r.squared))
        }, error = function(e) {
          cat("  ⚠ Quadratic refit failed:", e$message, "\n")
        })
      }
      
      # ===== INSTRUCTION 1: STEP B — Build best_tidy AFTER model selection is final =====
      best_tidy <- tidy(best_model) %>%
        mutate(
          aerosol_type = aerosol,
          ae_category = as.character(ae_cat),
          model_spec = best_model_spec,
          aic = AIC(best_model),
          bic = BIC(best_model),
          r_squared = summary(best_model)$r.squared,
          adj_r_squared = summary(best_model)$adj.r.squared,
          n_obs = nrow(analysis_df_scaled),
          max_vif = all_aic_results$max_vif[best_model_idx],
          delta_aic_nonlinear = nonlin_result$delta_aic,
          nonlinear_preferred = nonlin_result$nonlinear_preferred,
          model_is_quadratic  = model_is_quadratic,
          single_scene_flag = single_scene_flag,
          n_scenes_in_stratum = n_scenes_in_stratum,
          weighted_model = !all(obs_weights == 1),
          # Placeholders — filled in by subsequent diagnostic blocks
          morans_i = NA_real_,
          morans_p = NA_real_,
          morans_k = NA_real_,
          morans_n_used = NA_real_,
          scene_consistency_f = NA_real_,
          scene_consistency_p = NA_real_,
          slopes_consistent = NA,
          cv_rmse_mean = NA_real_,
          cv_rmse_sd = NA_real_,
          cv_r2_mean = NA_real_,
          n_blocks_used = NA_integer_,
          cv_warning = NA
        )
      
      # ===== INSTRUCTION 1: STEP C — Generate residual diagnostic plots (on correct model) =====
      tryCatch({
        png_diag <- file.path(plot_dir, paste0("Diagnostics_", gsub(" ", "_", aerosol), "_", as.character(ae_cat), ".png"))
        png(png_diag, width = 1200, height = 1000)
        par(mfrow = c(2, 2), mar = c(5, 5, 3, 2), oma = c(1, 1, 2, 1))
        plot(best_model, main = "")
        mtext(paste("Residual Diagnostics:", aerosol, "-", as.character(ae_cat)), 
              outer = TRUE, line = 0.5, cex = 1.5, font = 2)
        dev.off()
        cat("✓ Diagnostic plots saved to:", png_diag, "\n")
      }, error = function(e) {
        cat("⚠ Warning: Could not generate diagnostic plots -", e$message, "\n")
      })
      
      # ===== INSTRUCTION 1: STEP D — Moran's I test =====
      tryCatch({
        residuals_best <- residuals(best_model)
        coords_best <- coords_mat
        
        if (nrow(coords_best) != length(residuals_best)) {
          warning("Coordinate mismatch with residuals; Moran's I skipped")
          best_tidy$morans_i <- NA
          best_tidy$morans_p <- NA
          best_tidy$morans_k <- NA
          best_tidy$morans_n_used <- NA
        } else if (sum(is.na(residuals_best)) > length(residuals_best) * 0.5) {
          warning("Too many NA residuals; Moran's I skipped")
          best_tidy$morans_i <- NA
          best_tidy$morans_p <- NA
          best_tidy$morans_k <- NA
          best_tidy$morans_n_used <- NA
        } else {
          morans_result <- calculate_morans_i_knn(residuals_best, coords_best, k = 8, max_n = 2000)
          
          cat("\n📊 MORAN'S I TEST FOR SPATIAL AUTOCORRELATION (KNN-based):\n")
          cat(sprintf("  Moran's I statistic: %.6f\n", morans_result$I))
          cat(sprintf("  Expected value (under null): %.6f\n", morans_result$expected))
          if (!is.null(morans_result$subsampled) && morans_result$subsampled) {
            cat(sprintf("  NOTE: Subsampled to %d obs (max_n=2000)\n", morans_result$n_used))
          }
          cat(sprintf("  k-nearest neighbors: %d\n", morans_result$k))
          cat(sprintf("  Obs used in test: %d\n", morans_result$n_used))
          cat(sprintf("  Permutation-based p-value: %.2e\n", morans_result$p))
          
          if (morans_result$p < 0.05) {
            cat("  ⚠️  SIGNIFICANT spatial autocorrelation detected (p < 0.05)!\n")
            cat("     Residuals are spatially correlated; standard errors likely UNDERESTIMATED\n")
            cat("     Consider: spatial models, clustering of obs, GEE, or spatial lags\n")
          } else {
            cat("  ✓ No significant spatial autocorrelation detected (p >= 0.05)\n")
          }
          
          best_tidy$morans_i <- morans_result$I
          best_tidy$morans_p <- morans_result$p
          best_tidy$morans_k <- morans_result$k
          best_tidy$morans_n_used <- morans_result$n_used
        }
      }, error = function(e) {
        cat("WARNING: Could not compute Moran's I -", e$message, "\n")
        best_tidy$morans_i <<- NA
        best_tidy$morans_p <<- NA
        best_tidy$morans_k <<- NA
        best_tidy$morans_n_used <<- NA
      })
      
      # ===== INSTRUCTION 1: STEP E — Build CV formula explicitly (INSTRUCTION 2) =====
      if (single_scene_flag) {
        cv_formula <- as.formula(paste("DP ~", paste(best_combo, collapse = " + ")))
      } else {
        cv_formula <- as.formula(paste("DP ~", paste(best_combo, collapse = " + "), "+ as.factor(scene)"))
      }
      
      if (model_is_quadratic) {
        quad_terms_cv <- paste(paste0("I(", best_combo, "^2)"), collapse = " + ")
        base_terms    <- paste(best_combo, collapse = " + ")
        if (single_scene_flag) {
          cv_formula <- as.formula(paste("DP ~", base_terms, "+", quad_terms_cv))
        } else {
          cv_formula <- as.formula(paste("DP ~", base_terms, "+", quad_terms_cv, "+ as.factor(scene)"))
        }
      }
      
      # ===== INSTRUCTION 4: CROSS-SCENE CONSISTENCY TEST =====
      cat("\n🔍 CROSS-SCENE CONSISTENCY TEST\n")
      
      if (single_scene_flag) {
        cat("  Skipped: single-scene stratum — consistency test requires ≥ 2 scenes.\n")
        scene_consistency <- list(
          f_stat           = NA_real_,
          p_value          = NA_real_,
          slopes_consistent = NA,
          note             = "Not applicable: single-scene stratum"
        )
      } else if (length(unique(analysis_df_scaled$scene)) < 2) {
        cat("  Skipped: only one scene level present after filtering.\n")
        scene_consistency <- list(
          f_stat           = NA_real_,
          p_value          = NA_real_,
          slopes_consistent = NA,
          note             = "Not applicable: one scene level after filtering"
        )
      } else {
        scene_consistency <- test_scene_consistency(best_model, best_combo, analysis_df_scaled, weights = obs_weights)
        cat(sprintf("  F-statistic: %.4f | p-value: %.2e\n",
                    scene_consistency$f_stat, scene_consistency$p_value))
        cat(sprintf("  %s\n", scene_consistency$note))
        
        if (isFALSE(scene_consistency$slopes_consistent)) {
          cat("  ⚠️  WARNING: Pooling assumption may be violated.\n")
        }
      }
      
      best_tidy$scene_consistency_f  <- scene_consistency$f_stat
      best_tidy$scene_consistency_p  <- scene_consistency$p_value
      best_tidy$slopes_consistent    <- scene_consistency$slopes_consistent
      
      # ===== INSTRUCTION 1: STEP F — Spatial block CV (with explicit formula) =====
      cat("\n📊 SPATIAL BLOCK CROSS-VALIDATION\n")
      cv_result <- spatial_block_cv(cv_formula, analysis_df_scaled, 
                                    coords_mat, weights = obs_weights, n_blocks = 5)
      
      cat(sprintf("  CV RMSE (mean ± sd): %.4f ± %.4f\n", 
                  cv_result$cv_rmse_mean, cv_result$cv_rmse_sd))
      cat(sprintf("  CV R²: %.4f\n", cv_result$cv_r2_mean))
      cat(sprintf("  Blocks evaluated: %d\n", cv_result$n_blocks_used))
      
      # Check if CV performance degrades substantially (>15% RMSE increase)
      if (!is.na(cv_result$cv_rmse_mean)) {
        in_sample_rmse <- sqrt(mean(residuals(best_model)^2, na.rm = TRUE))
        rmse_ratio <- cv_result$cv_rmse_mean / in_sample_rmse
        
        cat(sprintf("  In-sample RMSE: %.4f\n", in_sample_rmse))
        cat(sprintf("  RMSE ratio (CV/in-sample): %.2f\n", rmse_ratio))
        
        cv_warning <- rmse_ratio > 1.15
        if (cv_warning) {
          cat("  ⚠️  WARNING: CV performance degrades >15% vs in-sample (possible overfitting or spatial leak)\n")
        }
      } else {
        cv_warning <- NA
      }
      
      best_tidy$cv_rmse_mean <- cv_result$cv_rmse_mean
      best_tidy$cv_rmse_sd <- cv_result$cv_rmse_sd
      best_tidy$cv_r2_mean <- cv_result$cv_r2_mean
      best_tidy$n_blocks_used <- cv_result$n_blocks_used
      best_tidy$cv_warning <- cv_warning
      
      # ===== INSTRUCTION 1: STEP G — Store best_tidy AFTER all diagnostics complete =====
      # All NA placeholders have now been filled with real values; store the complete row
      analysis_results[[paste(aerosol, as.character(ae_cat), sep = "_")]] <- best_tidy
      
      # Store comparison table
      all_aic_results$aerosol_type <- aerosol
      all_aic_results$ae_category <- as.character(ae_cat)
      model_comparison_results[[paste(aerosol, as.character(ae_cat), sep = "_")]] <- all_aic_results
      
      # Create visualization of top models by AIC
      top_models <- head(all_aic_results, 5) %>%
        mutate(model_rank = paste0("#", aic_rank))
      
      gg_aic <- ggplot(top_models, aes(x = reorder(model_spec, aic), y = aic, fill = factor(n_predictors))) +
        geom_col() +
        coord_flip() +
        labs(
          title = paste("Top 5 Models by AIC for", aerosol, "-", as.character(ae_cat)),
          subtitle = paste("n =", nrow(analysis_df_scaled), "pixels"),
          x = "Model Specification", y = "AIC", fill = "# Predictors"
        ) +
        theme_bw()
      
      ggsave(
        filename = file.path(plot_dir, paste0("AIC_Models_", gsub(" ", "_", aerosol), "_", as.character(ae_cat), ".png")),
        plot = gg_aic,
        width = 10, height = 6
      )
      
      cat("✓ Best model:", best_model_spec, "with AIC =", round(all_aic_results$aic[best_model_idx], 2), "\n")
      cat("  R² =", round(all_aic_results$r_squared[best_model_idx], 4), "\n")
    }
  }
}

cat("\n--- Analysis Complete ---\n")

# -----------------------------------------------------------------------------
# =============================================================================
# Step 5: Consolidate and Save Final Results (with INSTRUCTION 6 & 7)
# =============================================================================
# Combine all results from the list into a single dataframe
if (length(analysis_results) > 0) {
  # Best model coefficients
  final_summary_df <- bind_rows(analysis_results)
  
  # ===== INSTRUCTION 2 FIX: Pre-compute DP SD lookup table =====
  # INSTRUCTION 6a: EFFECT SIZES (standardized by DP scale)
  # Pre-compute DP standard deviation per stratum (not inside group_by context)
  dp_sd_lookup <- all_pixels_df %>%
    group_by(Aerosol_Label, AE_Category) %>%
    summarise(dp_sd = sd(DP, na.rm = TRUE), .groups = "drop") %>%
    rename(aerosol_type = Aerosol_Label, ae_category = AE_Category) %>%
    mutate(
      aerosol_type = as.character(aerosol_type),
      ae_category  = as.character(ae_category)
    )
  
  # Join lookup table and compute normalized effect sizes
  final_summary_df <- final_summary_df %>%
    left_join(dp_sd_lookup, by = c("aerosol_type", "ae_category")) %>%
    mutate(
      effect_size_normalized = ifelse(
        !grepl("Intercept|scene", term),
        estimate / dp_sd,
        NA_real_
      )
    )
  
  # ===== INSTRUCTION 6b: EFFECTIVE SAMPLE SIZE CORRECTION =====
  # Account for spatial autocorrelation: n_eff ≈ n_obs * (1 - morans_i) / (1 + morans_i)
  final_summary_df <- final_summary_df %>%
    mutate(
      n_eff_approx = round(n_obs * (1 - morans_i) / (1 + morans_i)),
      n_eff_approx = ifelse(is.na(morans_i) | n_eff_approx < 1, n_obs, n_eff_approx),
      low_eff_n_warning = ifelse(!is.na(n_eff_approx) & n_eff_approx < 30, TRUE, FALSE)
    )
  
  # ===== INSTRUCTION 6c: SCENE COMPOSITION =====
  # For each aerosol/AE stratum, list contributing scenes and count
  # ===== INSTRUCTION 5: Fix factor/character mismatch =====
  scene_comp <- all_pixels_df %>%
    group_by(Aerosol_Label, AE_Category) %>%
    summarise(
      scenes_contributing = paste(unique(scene), collapse = ";"),
      n_scenes = n_distinct(scene),
      .groups = "drop"
    ) %>%
    rename(aerosol_type = Aerosol_Label, ae_category = AE_Category) %>%
    mutate(
      aerosol_type = as.character(aerosol_type),
      ae_category  = as.character(ae_category)
    )
  
  # Add scene composition and flag single-scene models
  final_summary_df <- final_summary_df %>%
    left_join(scene_comp, by = c("aerosol_type", "ae_category")) %>%
    mutate(
      single_scene_warning = ifelse(!is.na(n_scenes) & n_scenes == 1, TRUE, FALSE)
    )
  
  # ===== INSTRUCTION 5: Verify joins succeeded =====
  n_missing_dp_sd <- sum(is.na(final_summary_df$dp_sd))
  n_missing_scenes <- sum(is.na(final_summary_df$scenes_contributing))
  
  if (n_missing_dp_sd > 0) {
    cat(sprintf("⚠ WARNING: dp_sd join failed for %d rows — check factor/character mismatch.\n",
                n_missing_dp_sd))
  }
  if (n_missing_scenes > 0) {
    cat(sprintf("⚠ WARNING: scene_comp join failed for %d rows — check factor/character mismatch.\n",
                n_missing_scenes))
  }
  
  # Select columns in logical order with all new INSTRUCTION 6 fields
  final_summary_df <- final_summary_df %>%
    dplyr::select(
      aerosol_type, ae_category, model_spec, term, estimate, std.error, 
      effect_size_normalized, statistic, p.value, 
      aic, bic, r_squared, adj_r_squared, 
      cv_rmse_mean, cv_rmse_sd, cv_r2_mean, n_blocks_used, cv_warning,
      max_vif, 
      morans_i, morans_p, morans_k, morans_n_used,
      scene_consistency_f, scene_consistency_p, slopes_consistent,
      delta_aic_nonlinear, nonlinear_preferred, model_is_quadratic, weighted_model,
      n_obs, n_eff_approx, low_eff_n_warning,
      scenes_contributing, n_scenes, single_scene_warning
    )
  
  # Save the best models summary
  best_models_file <- file.path(plot_dir, "best_models_summary.csv")
  write.csv(final_summary_df, best_models_file, row.names = FALSE)
  message("\n✓ Best models summary saved to: ", best_models_file)
  cat("  Columns include: effect sizes, CV metrics, spatial autocorrelation, scene consistency, nonlinearity tests\n")
  print(head(final_summary_df, 15))
  
  # Model comparison results (all models per aerosol type and AE category)
  if (length(model_comparison_results) > 0) {
    all_models_comparison <- bind_rows(model_comparison_results) %>%
      dplyr::select(aerosol_type, ae_category, model_spec, n_predictors, aic, bic, 
                    r_squared, adj_r_squared, max_vif, weighted_model, n_obs, aic_rank) %>%
      arrange(aerosol_type, ae_category, aic_rank)
    
    comparison_file <- file.path(plot_dir, "model_comparison_all.csv")
    write.csv(all_models_comparison, comparison_file, row.names = FALSE)
    message("\n✓ Full model comparison saved to: ", comparison_file)
  }
  
  # Create summary statistics table (one row per aerosol/AE stratum)
  summary_stats <- final_summary_df %>%
    filter(!grepl("Intercept|scene", term)) %>%
    group_by(aerosol_type, ae_category, model_spec) %>%
    summarise(
      aic = first(aic),
      bic = first(bic),
      r_squared = first(r_squared),
      adj_r_squared = first(adj_r_squared),
      cv_rmse_mean = first(cv_rmse_mean),
      cv_r2_mean = first(cv_r2_mean),
      max_vif = first(max_vif),
      morans_i = first(morans_i),
      morans_p = first(morans_p),
      morans_k = first(morans_k),
      morans_n_used = first(morans_n_used),
      scene_consistency_f = first(scene_consistency_f),
      scene_consistency_p = first(scene_consistency_p),
      slopes_consistent = first(slopes_consistent),
      delta_aic_nonlinear = first(delta_aic_nonlinear),
      nonlinear_preferred = first(nonlinear_preferred),
      n_obs = first(n_obs),
      n_eff_approx = first(n_eff_approx),
      low_eff_n_warning = first(low_eff_n_warning),
      n_scenes = first(n_scenes),
      single_scene_warning = first(single_scene_warning),
      .groups = "drop"
    ) %>%
    arrange(aerosol_type, ae_category)
  
  summary_file <- file.path(plot_dir, "aerosol_bias_analysis_summary.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  message("\n✓ Analysis summary statistics saved to: ", summary_file)
  print(summary_stats)
  
  # ===== INSTRUCTION 8: Add weighted_model warning =====
  n_weighted   <- sum(final_summary_df$weighted_model, na.rm = TRUE)
  n_unweighted <- sum(!final_summary_df$weighted_model, na.rm = TRUE)
  
  if (n_weighted > 0 && n_unweighted > 0) {
    cat(sprintf(
      "\n⚠ WARNING: Model summary contains %d weighted and %d unweighted cases.\n",
      n_weighted, n_unweighted
    ))
    cat("  AIC values are NOT directly comparable across weighted/unweighted strata.\n")
    cat("  Cross-stratum AIC comparisons should be interpreted with caution.\n")
  }
  
  # ===== INSTRUCTION 7: CROSS-AEROSOL COEFFICIENT COMPARISON PLOT =====
  cat("\n--- INSTRUCTION 7: Generating Cross-Aerosol Comparison Plot ---\n")
  tryCatch({
    # Filter to predictor terms (exclude intercept and scene fixed effects)
    coef_plot_df <- final_summary_df %>%
      filter(!grepl("Intercept|scene", term)) %>%
      mutate(
        ci_low  = estimate - 1.96 * std.error,
        ci_high = estimate + 1.96 * std.error,
        significant = p.value < 0.05,
        label = paste0(aerosol_type, "\n(", ae_category, ")")
      )
    
    if (nrow(coef_plot_df) > 0) {
      gg_coef <- ggplot(coef_plot_df, aes(x = term, y = estimate, color = significant)) +
        geom_point(size = 3, position = position_dodge(width = 0.5)) +
        geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                      width = 0.2, position = position_dodge(width = 0.5)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
        facet_wrap(~ label, scales = "free_y") +
        scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue"),
                           labels = c("NS", "p < 0.05")) +
        coord_flip() +
        labs(
          title = "Standardized Predictor Coefficients by Aerosol Type and AE Category",
          subtitle = "Error bars = 95% CI. Note: SEs may be underestimated where Moran's I is significant.",
          x = "Predictor", y = "Standardized Coefficient (per 3 SD)",
          color = "Significance"
        ) +
        theme_bw(base_size = 12)
      
      ggsave(
        filename = file.path(plot_dir, "cross_aerosol_coefficient_comparison.png"),
        plot = gg_coef, width = 16, height = 10, dpi = 150
      )
      message("✓ Cross-aerosol coefficient comparison plot saved\n")
    }
  }, error = function(e) {
    cat("⚠ Warning: Could not generate cross-aerosol plot -", e$message, "\n")
  })
  
} else {
  message("\nNo analysis was performed, likely due to lack of data.")
}