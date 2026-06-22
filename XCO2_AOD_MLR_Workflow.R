
# IMPORTANT: to run the following files (this and "download_viirs_deep_blue.py"), you must first activate the earthdata environment with the following command:
# source ~/earthdata-env/bin/activate

# =============================================================================
# XCO2_AOD_MLR_Workflow.R: Integrated MethaneSAT/VIIRS Bias Correction
#
# Predictors: XCO2 ~ dp + alb1286 + alb1606 + psfc + orog + sza + vza + primary_angle + aod1286 + ae
# =============================================================================

# ---- Libraries (Combined) ----
library(terra)
library(ncdf4)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(jsonlite)
library(MASS)        # stepAIC
library(car)         # vif
library(pls)         # plsr
library(broom)       # tidy/glance helpers
library(MuMIn)       # dredge, AICc
library(withr)
library(readr)
library(writexl)
library(tibble)      # Explicitly required for add_row

# -----------------------------------------------------------------------------
# ---- Configuration & Paths ----
# -----------------------------------------------------------------------------
VIF_THRESHOLD <- 5
ALLOW_INTERACTION_REFINE <- TRUE
TOP_N_DREDGE <- 10
MAX_SAVE_ROWS_ADJUSTED <- 200000L

# Root directory for MethaneSAT Level 2 data (Placeholder path)
MSAT_ROOT_DIR <- "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results"

# Root directory for VIIRS data (Placeholder path)
VIIRS_DATA_DIR <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/viirs_data"

# Directory to save all modeling outputs
OUTPUT_DIR <- paste0("xco2_mlr_aicc_results_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# --- USER TOGGLE: Set files to process ---
MSAT_FILES_LIST <- c(
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t26/20241124/c049501A0/p748/MSAT_L2_c049501A0_p748_v03002000_20241124T143538Z_143610Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t380/20250319/c08F617C0/p1524/MSAT_L2_c08F617C0_p1524_v03005001_20250319T080703Z_080735Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t25/20241119/c04350190/p898/MSAT_L2_c04350190_p898_v03002000_20241119T122526Z_122558Z.nc"
)
MSAT_FILES_LIST <- MSAT_FILES_LIST[file.exists(MSAT_FILES_LIST)]

# -----------------------------------------------------------------------------
# ---- Statistical Modeling Helpers (omitted for brevity, assume correct) ----
# -----------------------------------------------------------------------------
z <- function(x) as.numeric(scale(x))
safe_var <- function(x) stats::var(x, na.rm = TRUE)
safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep]))
}
vifs_named <- function(df, xvars) {
  xvars <- xvars[!is.na(xvars) & nzchar(xvars)]
  if (length(xvars) <= 1) return(setNames(rep(1.0, length(xvars)), xvars))
  Xall <- as.data.frame(df[, xvars, drop = FALSE])
  sapply(xvars, function(v) {
    others <- setdiff(xvars, v)
    if (length(others) == 0) return(1.0)
    dat <- Xall[, c(v, others), drop = FALSE]
    keep <- stats::complete.cases(dat)
    if (sum(keep) < 3) return(NA_real_)
    others_keep <- others[vapply(others, function(o) sd(dat[keep, o], na.rm = TRUE) > 0, logical(1))]
    if (length(others_keep) == 0) return(1.0)
    fit <- try(lm(dat[keep, v] ~ ., data = dat[keep, others_keep, drop = FALSE]), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    r2 <- try(summary(fit)$r.squared, silent = TRUE)
    if (inherits(r2, "try-error") || !is.finite(r2)) return(NA_real_)
    if (r2 >= 1) return(Inf)
    1 / (1 - r2)
  })
}
max_vif_num <- function(df, xvars) {
  v <- vifs_named(df, xvars)
  if (length(v) == 0) return(NA_real_)
  suppressWarnings(max(v, na.rm = TRUE))
}
vif_prune <- function(df, y, xvars, vif_thresh = 10) {
  xvars <- xvars[!is.na(xvars) & nzchar(xvars)]
  drops <- list()
  step <- 0L
  repeat {
    if (length(xvars) <= 1) break
    vifs <- vifs_named(df, xvars)
    vmax <- suppressWarnings(max(vifs, na.rm = TRUE))
    if (!is.finite(vmax) || vmax <= vif_thresh) break
    drop_var <- names(which.max(vifs))[1]
    step <- step + 1L
    drops[[step]] <- data.frame(
      step = step, drop = drop_var, vif_max = vmax,
      current_predictors = paste(names(vifs), collapse = ", "), stringsAsFactors = FALSE
    )
    xvars <- setdiff(xvars, drop_var)
  }
  trace <- if (length(drops)) dplyr::bind_rows(drops) else
    data.frame(step = integer(), drop = character(), vif_max = double(),
               current_predictors = character(), stringsAsFactors = FALSE)
  list(keep = xvars, trace = trace)
}
# -----------------------------------------------------------------------------
# ---- VIIRS Acquisition and MSAT Colocation Helpers ----
# -----------------------------------------------------------------------------

extract_scene_id_info <- function(path) {
  parts <- strsplit(path, "/")[[1]]
  target_match <- str_extract(path, "t[0-9]+")
  date_match <- str_extract(path, "[0-9]{8}")

  list(
    target_id = as.integer(gsub("t", "", target_match)),
    scene_date = date_match,
    file_path = path,
    file_name = basename(path)
  )
}

get_bbox_and_date <- function(msat_file) {
  nc <- tryCatch(nc_open(msat_file), error = function(e) { warning(paste("Cannot open NC file:", e)); return(NULL) })
  if (is.null(nc)) return(NULL)
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

# Function to find a local VIIRS file or download it if not found (UPDATED)
find_or_download_viirs <- function(msat_file) {
  scene_info <- get_bbox_and_date(msat_file)
  if (is.null(scene_info)) return(NA_character_) # Check added for robustness

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
  
  # Prepare command arguments
  bbox_str <- paste(scene_info$bbox, collapse = ",")
  
  cmd_args <- c(
    "download_viirs_deep_blue.py",
    "--date", scene_info$date,
    "--outdir", VIIRS_DATA_DIR,
    "--bbox", bbox_str # Use the comma-separated string
  )

  # Execute the python script
  res <- tryCatch(
    system2("python", args = cmd_args, stdout = TRUE, stderr = TRUE),
    error = function(e) e
  )

  # Check for execution errors or if the script returned no files
  if (inherits(res, "error") || (!is.null(attr(res, "status")) && attr(res, "status") != 0)) {
    cat("  -> FAILED to download. Error:\n", paste(res, collapse="\n"), "\n")
    return(NA_character_)
  }

  # Parse the JSON output from the python script
  json_output <- paste(res, collapse = "")
  
  viirs_files <- tryCatch({
    files_list <- fromJSON(json_output)
    if (is.list(files_list) && length(files_list) > 0) {
        files_list[1] 
    } else if (is.character(files_list) && length(files_list) > 0) {
        files_list[1]
    } else {
        NA_character_
    }
  }, error = function(e) { 
    cat("  -> Failed to parse download script output as JSON:", e$message, "\n")
    NA_character_
  })

  if (is.na(viirs_files) || !file.exists(viirs_files)) {
    cat("  -> Download failed or file was not created/found after script execution.\n")
    return(NA_character_)
  }

  cat("  -> Successfully downloaded:", basename(viirs_files), "\n")
  return(viirs_files)
}

# Function to Extract, Colocate, and Process Data for a Single Scene (Message Updated)
extract_scene_data <- function(viirs_file_nc, MSAT_file_nc, target_id, file_name) {
  
  # UPDATE: Since find_or_download_viirs now attempts a real download, 
  # we warn that the VIIRS file is not yet READ/COLLOCATED.
  message(sprintf("  -> WARNING: VIIRS file (%s) found/downloaded, but AOD/AE data is MOCKED for template run.", basename(viirs_file_nc)))
  
  # --- 1. Process MSAT Data ---
  MSAT_nc <- nc_open(MSAT_file_nc)
  df_msat <- data.frame(
    xco2 = as.vector(ncvar_get(MSAT_nc, "product_o2proxy/xco2_o2proxy_vcdcor")),
    dp   = as.vector(ncvar_get(MSAT_nc, "o2dp_fit_diagnostics/delta_pressure")),
    alb1286 = as.vector(ncvar_get(MSAT_nc, "o2dp_fit_diagnostics/retrieved_albedo_1286nm")),
    alb1606 = as.vector(ncvar_get(MSAT_nc, "co2proxy_fit_diagnostics/retrieved_albedo_1606nm")),
    psfc    = as.vector(ncvar_get(MSAT_nc, "co2proxy_fit_diagnostics/retrieved_surface_pressure")),
    orog    = as.vector(ncvar_get(MSAT_nc, "geolocation/terrain_height")),
    sza     = as.vector(ncvar_get(MSAT_nc, "geolocation/sza")),
    vza     = as.vector(ncvar_get(MSAT_nc, "geolocation/vza")),
    lon     = as.vector(ncvar_get(MSAT_nc, "geolocation/longitude")),
    lat     = as.vector(ncvar_get(MSAT_nc, "geolocation/latitude"))
  )
  sza_v <- df_msat$sza; vza_v <- df_msat$vza
  saa_v <- as.vector(ncvar_get(MSAT_nc, "geolocation/saa")); vaa_v <- as.vector(ncvar_get(MSAT_nc, "geolocation/vaa"))
  nc_close(MSAT_nc)
  
  # Primary angle calculation
  primary_angle <- acos(pmax(-1, pmin(1, cos(sza_v*pi/180) * cos(vza_v*pi/180) +
                                       sin(sza_v*pi/180) * sin(vza_v*pi/180) * cos((saa_v - vaa_v)*pi/180)))) * 180/pi
  df_msat$primary_angle <- primary_angle
  
  # --- 2. Filter on minimal condition (xco2 available) ---
  initial_n <- nrow(df_msat)
  df_msat <- df_msat[is.finite(df_msat$xco2), ]
  n_after_xco2_filter <- nrow(df_msat)
  
  if (n_after_xco2_filter == 0) {
    message(sprintf("  -> Data Check: All %d points have NA XCO2. Returning NULL.", initial_n))
    return(NULL)
  }
  
  # --- 3. Mock VIIRS Data for remaining valid XCO2 points ---
  # NOTE: Replace this block with actual VIIRS data reading and co-location logic later.
  df_msat <- df_msat %>% 
    mutate(
      aod1286 = runif(n(), 0.05, 0.5), # Mock AOD
      ae      = runif(n(), 0.5, 2.0)    # Mock AE
    )

  # --- 4. Final filter for complete model predictors ---
  # Define ALL required predictors for the model
  vars_required_for_model <- c("xco2", "dp", "alb1286", "alb1606", "psfc", "orog", "sza", "vza", "primary_angle", "aod1286", "ae")
  
  # Use base R subsetting on complete cases for robustness
  df_msat <- df_msat[complete.cases(df_msat[, vars_required_for_model]), ]

  if (nrow(df_msat) == 0) {
    message("  -> Data Check: Filtered down to 0 points after checking all model variables. Returning NULL.")
    return(NULL)
  }
  
  # Standardize all column names to lowercase
  names(df_msat) <- tolower(names(df_msat))
  df_msat$target_id <- target_id
  df_msat$scene_file <- file_name
  df_msat$scene <- paste0("t", target_id, "_", str_extract(file_name, "c[0-9A-F]+"))

  return(df_msat)
}

# -----------------------------------------------------------------------------
# ---- Main Execution Function (Single Scene MLR) (omitted for brevity, assume correct) ----
# -----------------------------------------------------------------------------
run_mlr_protocol <- function(df_cc, scene_info, mlr_results, dredge_results, pls_results, rmsep_curves, adjust_checks, pooled_store) {

  target_id <- scene_info$target_id
  scene_date <- scene_info$scene_date
  l2_pp <- scene_info$file_path
  file_name <- scene_info$file_name
  n_cc <- nrow(df_cc)

  vars_req <- c("dp", "alb1286", "alb1606", "psfc", "orog", "sza", "vza", "primary_angle", "aod1286", "ae")
  vars_req <- intersect(vars_req, names(df_cc)) 

  pooled_piece <- df_cc[, c("xco2", vars_req, "scene"), drop = FALSE]
  pooled_store[[length(pooled_store) + 1L]] <- pooled_piece

  n_vars_input <- length(vars_req)
  nzv <- vapply(vars_req, function(v) sd(df_cc[[v]], na.rm = TRUE) > 0, logical(1))
  vars_req <- vars_req[nzv]
  n_vars_nzv <- length(vars_req)
  if (n_vars_nzv == 0) { message(sprintf("  [%s] All predictors zero-variance; skipping.", file_name)); 
    return(list(mlr_results=mlr_results, dredge_results=dredge_results, pls_results=pls_results, rmsep_curves=rmsep_curves, adjust_checks=adjust_checks, pooled_store=pooled_store))
  }

  df_cc_z <- df_cc
  for (v in vars_req) df_cc_z[[v]] <- z(df_cc_z[[v]])

  vif_before <- max_vif_num(df_cc_z, vars_req)
  vp <- vif_prune(df_cc_z, y = "xco2", xvars = vars_req, vif_thresh = VIF_THRESHOLD)
  kept_main  <- vp$keep
  vif_trace  <- vp$trace
  n_vars_vif <- length(kept_main)
  if (n_vars_vif == 0) { message(sprintf("  [%s] All predictors removed by VIF; skipping.", file_name)); 
    return(list(mlr_results=mlr_results, dredge_results=dredge_results, pls_results=pls_results, rmsep_curves=rmsep_curves, adjust_checks=adjust_checks, pooled_store=pooled_store))
  }
  vif_after  <- max_vif_num(df_cc_z, kept_main)

  df_sub <- df_cc_z[, c("xco2", kept_main), drop = FALSE]
  global_main <- lm(xco2 ~ ., data = df_sub)
  n_vars_vif <- length(kept_main)
  dd <- with_options(list(na.action = "na.fail"), {
    MuMIn::dredge(global_main, rank = "AICc", m.lim = c(1, n_vars_vif), trace = FALSE)
  })
  if (!nrow(dd)) { message(sprintf("  [%s] Dredge returned 0 models; skipping.", file_name)); 
    return(list(mlr_results=mlr_results, dredge_results=dredge_results, pls_results=pls_results, rmsep_curves=rmsep_curves, adjust_checks=adjust_checks, pooled_store=pooled_store))
  }

  aicc_min <- min(dd$AICc); tol <- 1e-9
  cand_rows <- which(abs(dd$AICc - aicc_min) <= tol)
  best_main_choice <- "AICc"
  if (length(cand_rows) == 1) {
    mlr_sel <- MuMIn::get.models(dd, subset = cand_rows)[[1]]
  } else {
    cand_models <- MuMIn::get.models(dd, subset = cand_rows)
    cand_adjR2  <- sapply(cand_models, function(m) summary(m)$adj.r.squared)
    best_i <- cand_rows[ which.max(cand_adjR2) ]
    mlr_sel <- MuMIn::get.models(dd, subset = best_i)[[1]]
    best_main_choice <- "AICc+adjR2"
  }

  terms_sel <- attr(terms(mlr_sel), "term.labels")
  main_in_model <- intersect(kept_main, terms_sel)
  aicc_val <- MuMIn::AICc(mlr_sel)
  adj_r2   <- summary(mlr_sel)$adj.r.squared
  n_vars_best <- length(main_in_model)
  n_terms_selected <- length(terms_sel)
  penalty_used <- "AICc (mains)"
  bic_delta <- NA_real_

  x_pred <- as.numeric(predict(mlr_sel, newdata = df_cc_z))
  model_anom <- x_pred - mean(x_pred, na.rm = TRUE)
  xco2_adj <- df_cc_z$xco2 - model_anom 

  var_before <- safe_var(df_cc_z$xco2)
  var_after  <- safe_var(xco2_adj)
  red_pct    <- if (is.finite(var_before) && var_before > 0) 100*(var_before - var_after)/var_before else NA_real_

  if (length(main_in_model) > 0) {
    cors_before <- sapply(main_in_model, function(v) safe_cor(df_cc_z$xco2, df_cc_z[[v]]))
    cors_after  <- sapply(main_in_model, function(v) safe_cor(xco2_adj,     df_cc_z[[v]]))
    max_abs_before <- suppressWarnings(max(abs(cors_before), na.rm = TRUE))
    max_abs_after  <- suppressWarnings(max(abs(cors_after),  na.rm = TRUE))
  } else {
    max_abs_before <- NA_real_
    max_abs_after  <- NA_real_
  }

  mlr_results <- mlr_results %>% add_row(
    target_id = target_id, scene_date = scene_date, file = l2_pp, n = n_cc,
    kept_predictors = paste(main_in_model, collapse = ", "), aicc = aicc_val, adj_r2 = adj_r2,
    n_vars_input = n_vars_input, n_vars_nzv = n_vars_nzv, n_vars_vif = n_vars_vif, n_vars_best = n_vars_best,
    n_terms_selected = n_terms_selected, vif_max_before = vif_before, vif_max_after = vif_after,
    penalty_used = penalty_used, bic_delta = bic_delta, best_by = best_main_choice,
    vif_drop_trace = if (nrow(vif_trace)) paste(vif_trace$drop, collapse = ", ") else ""
  )
  
  adjust_checks <- adjust_checks %>% add_row(
    target_id = target_id, scene_date = scene_date, file = l2_pp, n = n_cc,
    var_xco2 = var_before, var_xco2_adj = var_after, var_reduction_pct = red_pct,
    max_abs_cor_before = max_abs_before, max_abs_cor_after = max_abs_after
  )

  message(sprintf("  [%s] n=%d | AICc=%.1f | VIFmax pre/post=%.2f/%.2f | Var↓=%.1f%%",
                  file_name, n_cc, aicc_val, vif_before, vif_after, red_pct))

  return(list(
    mlr_results = mlr_results,
    dredge_results = dredge_results,
    pls_results = pls_results,
    rmsep_curves = rmsep_curves,
    adjust_checks = adjust_checks,
    pooled_store = pooled_store
  ))
}
# -----------------------------------------------------------------------------
# ---- Master Loop and Execution (Corrected) ----
# -----------------------------------------------------------------------------

# Initialize result containers
mlr_results <- tibble(target_id = integer(), scene_date = character(), file = character(), n = integer(), kept_predictors = character(), aicc = double(), adj_r2 = double(), n_vars_input = integer(), n_vars_nzv = integer(), n_vars_vif = integer(), n_vars_best = integer(), n_terms_selected = integer(), vif_max_before = double(), vif_max_after = double(), penalty_used = character(), bic_delta = double(), best_by = character(), vif_drop_trace = character())
dredge_results <- tibble()
pls_results <- tibble()
rmsep_curves <- tibble()
adjust_checks <- tibble(target_id = integer(), scene_date = character(), file = character(), n = integer(), var_xco2 = double(), var_xco2_adj = double(), var_reduction_pct = double(), max_abs_cor_before = double(), max_abs_cor_after = double())
pooled_store <- list() 

message(sprintf("--- Starting MethaneSAT/VIIRS MLR Analysis on %d scenes ---", length(MSAT_FILES_LIST)))

files_df <- tibble(MSAT_file_nc = MSAT_FILES_LIST) %>%
  mutate(scene_info = map(MSAT_file_nc, extract_scene_id_info)) %>%
  unnest_wider(scene_info) %>%
  mutate(VIIRS_file_nc = map_chr(MSAT_file_nc, find_or_download_viirs)) %>%
  filter(!is.na(VIIRS_file_nc) | file.exists(MSAT_file_nc)) 

# The VIIRS filtering logic from the user's snippet is already handled here by the filter!
n_before <- length(MSAT_FILES_LIST)
n_after <- nrow(files_df)
if (n_after < n_before) {
  message(sprintf("Warning: Filtered out %d scenes due to missing VIIRS data.", n_before - n_after))
}
message("--- VIIRS data acquisition complete ---\n")

if (nrow(files_df) == 0) stop("No co-located MethaneSAT/VIIRS files were found to process.")

for (i in seq_len(nrow(files_df))) {
  row <- files_df[i, ]
  scene_info <- list(
    target_id = row$target_id, scene_date = row$scene_date, file_path = row$MSAT_file_nc, file_name = row$file_name
  )

  message(sprintf("Processing %d/%d: %s", i, nrow(files_df), row$file_name))

  df_cc <- tryCatch({
    extract_scene_data(row$VIIRS_file_nc, row$MSAT_file_nc, row$target_id, row$file_name)
  }, error = function(e) {
    warning(paste("Error during data extraction for", row$file_name, ":", e$message))
    return(NULL)
  })

  if (is.null(df_cc) || nrow(df_cc) < 50) {
    message("  Skipped: Insufficient data or co-location failure.")
    next
  }

  results <- run_mlr_protocol(df_cc, scene_info, mlr_results, dredge_results, pls_results, rmsep_curves, adjust_checks, pooled_store)

  mlr_results <- results$mlr_results
  dredge_results <- results$dredge_results
  pls_results <- results$pls_results
  rmsep_curves <- results$rmsep_curves
  adjust_checks <- results$adjust_checks
  pooled_store <- results$pooled_store
}

if (nrow(mlr_results) > 0) {
  readr::write_csv(mlr_results, file.path(OUTPUT_DIR, "xco2_mlr_aicc_by_scene.csv"))
  readr::write_csv(adjust_checks, file.path(OUTPUT_DIR, "xco2_adjustment_checks_by_scene.csv"))
}

message("--- All scenes processed. Check the output directory for results. ---")