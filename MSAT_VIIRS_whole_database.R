# =============================================================================
# Predictive pipeline (Ranger Random Forest) for MethaneSAT + VIIRS unified DB
# =============================================================================
# - Expects: `all_scenes_df` in the R session OR a CSV file named "unified_database.csv"
# - Trains a ranger random forest on a large sample
# - Predicts in chunks across the full dataset and writes predictions to disk
# - Saves model, importance, and diagnostics
# =============================================================================

# -----------------------
# Libraries
# -----------------------
library(dplyr)
library(readr)
library(data.table)   # fwrite for fast chunked writes
library(ranger)
library(ggplot2)
library(Metrics)      # rmse, mae
library(glue)

# -----------------------
# Configuration (editable)
# -----------------------
OUTPUT_DIR <- "results_predictive_ranger"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# If all_scenes_df is not in the environment, try to find a CSV
DB_PATH <- NA_character_   # if NA, script will search for "unified_database.csv" files under wd

# Model sampling / training sizes
TRAIN_SAMPLE_N <- 1000000   # 1e6 default; increase up to 3e6 if you have memory/CPU
TEST_FRAC <- 0.2            # fraction of sample used as test set
RANGER_NTREES <- 800
RANGER_MTRY <- NULL         # NULL -> auto (sqrt(p) default used below)
RANGER_MIN_NODE_SIZE <- 5

# Chunk size for full-dataset predictions (rows per chunk)
PRED_CHUNK_ROWS <- 500000

SEED <- 42
set.seed(SEED)

# -----------------------
# Predictors (as confirmed)
# -----------------------
PREDICTORS <- c("dp", "alb1286", "psfc", "sza", "vza", "primary_angle", "AOD_1286", "AE", "orog")
TARGET <- "xco2"

# -----------------------
# 0. Load or check unified dataframe
# -----------------------
if (!exists("all_scenes_df")) {
  message("all_scenes_df not present in environment. Searching for CSV...")
  if (is.na(DB_PATH)) {
    csv_candidates <- list.files(pattern = "unified_database\\.csv$", recursive = TRUE, full.names = TRUE)
    if (length(csv_candidates) == 0) {
      stop("Could not find all_scenes_df in the session and no 'unified_database.csv' found. Either load 'all_scenes_df' or set DB_PATH to a file path.")
    }
    # pick the first found
    DB_PATH <- csv_candidates[1]
    message("Using CSV: ", DB_PATH)
  }
  # load with data.table fread for speed
  all_scenes_df <- fread(DB_PATH, data.table = FALSE)
  message("Loaded dataframe with ", nrow(all_scenes_df), " rows and ", ncol(all_scenes_df), " columns.")
} else {
  message("Using all_scenes_df from environment: ", nrow(all_scenes_df), " rows.")
}

# -----------------------
# 1. Basic checks & cleaning
# -----------------------
# Keep only needed columns (target + predictors) to save memory
cols_needed <- unique(c(TARGET, PREDICTORS, "scene_id"))
miss_cols <- setdiff(c(TARGET, PREDICTORS), names(all_scenes_df))
if (length(miss_cols) > 0) {
  stop("Missing predictor columns in all_scenes_df: ", paste(miss_cols, collapse = ", "))
}

df_work <- all_scenes_df %>% dplyr::select(any_of(c(TARGET, PREDICTORS, "scene_id")))

# Drop rows with NA in target or predictors
rows_before <- nrow(df_work)
df_work <- df_work %>% drop_na(any_of(c(TARGET, PREDICTORS)))
rows_after <- nrow(df_work)
message(glue("Dropped {rows_before - rows_after} rows with NA in target/predictors. Remaining rows: {rows_after}."))

if (nrow(df_work) < 1000) stop("Not enough rows after cleaning to train a model.")

# -----------------------
# 2. Build sample for training (random sample)
# -----------------------
n_total <- nrow(df_work)
sample_n <- min(TRAIN_SAMPLE_N, n_total)
message("Drawing training+test sample of ", sample_n, " rows (seed=", SEED, ").")
sample_idx <- sample.int(n_total, sample_n)
df_sample <- df_work[sample_idx, , drop = FALSE]

# Train/test split
n_test <- floor(TEST_FRAC * nrow(df_sample))
n_train <- nrow(df_sample) - n_test
train_idx <- seq_len(n_train)
# shuffle sample for randomness
df_sample <- df_sample[sample(nrow(df_sample)), ]
df_train <- df_sample[1:n_train, ]
df_test  <- df_sample[(n_train+1):nrow(df_sample), ]

message(glue("Training rows: {nrow(df_train)}; Test rows: {nrow(df_test)}"))

# Optional: scale / transform predictors here if you want (Ranger handles raw variables well).
# Example: log-transform AOD if heavily skewed (uncomment if desired)
# df_train$AOD_1286 <- log1p(df_train$AOD_1286); df_test$AOD_1286 <- log1p(df_test$AOD_1286)

# -----------------------
# 3. Ranger model training
# -----------------------
p <- length(PREDICTORS)
if (is.null(RANGER_MTRY)) RANGER_MTRY <- floor(sqrt(p))  # typical default for regression forests

form_rf <- as.formula(paste(TARGET, "~", paste(PREDICTORS, collapse = " + ")))
message("Training ranger model: ", deparse(form_rf))

rf_fit <- ranger::ranger(
  formula = form_rf,
  data = df_train,
  num.trees = RANGER_NTREES,
  mtry = RANGER_MTRY,
  min.node.size = RANGER_MIN_NODE_SIZE,
  importance = "permutation",
  respect.unordered.factors = "order",
  num.threads = parallel::detectCores(),
  verbose = TRUE
)

# Save model
saveRDS(rf_fit, file = file.path(OUTPUT_DIR, "ranger_model.rds"))
message("Saved ranger model to ", file.path(OUTPUT_DIR, "ranger_model.rds"))

# -----------------------
# 4. Evaluate on test set
# -----------------------
pred_test <- predict(rf_fit, data = df_test)$predictions
rmse_test <- rmse(df_test[[TARGET]], pred_test)
mae_test  <- mae(df_test[[TARGET]], pred_test)
r2_test   <- 1 - sum((df_test[[TARGET]] - pred_test)^2) / sum((df_test[[TARGET]] - mean(df_test[[TARGET]]))^2)

message(glue("Test RMSE: {round(rmse_test,3)}; MAE: {round(mae_test,3)}; R2: {round(r2_test,3)}"))

# Save diagnostics
diag_tbl <- tibble::tibble(
  analysis = "ranger_sample_eval",
  timestamp = Sys.time(),
  n_train = nrow(df_train),
  n_test = nrow(df_test),
  rmse_test = rmse_test,
  mae_test = mae_test,
  r2_test = r2_test,
  num_trees = RANGER_NTREES,
  mtry = RANGER_MTRY
)
write_csv(diag_tbl, file.path(OUTPUT_DIR, "ranger_diagnostics.csv"))
message("Saved diagnostics to ", file.path(OUTPUT_DIR, "ranger_diagnostics.csv"))

# -----------------------
# 5. Variable importance & plot
# -----------------------
varimp <- data.frame(
  variable = names(rf_fit$variable.importance),
  importance = as.numeric(rf_fit$variable.importance),
  row.names = NULL,
  stringsAsFactors = FALSE
) %>% arrange(desc(importance))

write_csv(varimp, file.path(OUTPUT_DIR, "ranger_variable_importance.csv"))
message("Saved variable importance to ", file.path(OUTPUT_DIR, "ranger_variable_importance.csv"))

# Importance plot
p_imp <- ggplot(varimp, aes(x = reorder(variable, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  labs(title = "Ranger Variable Importance (permutation)", x = "", y = "Importance") +
  theme_minimal()
ggsave(file.path(OUTPUT_DIR, "ranger_variable_importance.png"), p_imp, width = 7, height = 4, dpi = 300)
message("Saved importance plot to ", file.path(OUTPUT_DIR, "ranger_variable_importance.png"))

# -----------------------
# 6. Predict on full dataset in chunks (memory-safe)
# -----------------------
message("Begin chunked predictions on full dataset...")

# Prepare output file: write header with id columns if available
pred_out_path <- file.path(OUTPUT_DIR, "unified_predictions.csv")
if (file.exists(pred_out_path)) file.remove(pred_out_path)  # overwrite

# We'll write columns: row_id, scene_id (if present), xco2, pred_xco2, pred_residual
n_rows <- nrow(df_work)
n_chunks <- ceiling(n_rows / PRED_CHUNK_ROWS)
message(glue("Total rows = {n_rows}; chunk rows = {PRED_CHUNK_ROWS}; chunks = {n_chunks}"))

start_row <- 1L
for (chunk_i in seq_len(n_chunks)) {
  end_row <- min(start_row + PRED_CHUNK_ROWS - 1L, n_rows)
  df_chunk <- df_work[start_row:end_row, , drop = FALSE]

  preds_chunk <- predict(rf_fit, data = df_chunk)$predictions

  out_chunk <- data.frame(
    row_id = start_row:end_row,
    scene_id = if ("scene_id" %in% names(df_chunk)) df_chunk$scene_id else NA,
    xco2 = df_chunk[[TARGET]],
    pred_xco2 = preds_chunk,
    residual = df_chunk[[TARGET]] - preds_chunk
  )

  # Append to CSV quickly using data.table::fwrite (fast and supports append)
  fwrite(out_chunk, pred_out_path, append = file.exists(pred_out_path), col.names = !file.exists(pred_out_path))

  message(glue("  Wrote chunk {chunk_i}/{n_chunks} rows {start_row}-{end_row}"))
  start_row <- end_row + 1L
}

message("All predictions written to: ", pred_out_path)

# -----------------------
# 7. Summary evaluation on full dataset (approx)
# -----------------------
# Compute overall RMSE/MAE on entire dataset by reading back the predictions file (fast)
preds_all <- fread(pred_out_path, select = c("xco2", "pred_xco2"))
rmse_full <- rmse(preds_all$xco2, preds_all$pred_xco2)
mae_full  <- mae(preds_all$xco2, preds_all$pred_xco2)
r2_full   <- 1 - sum((preds_all$xco2 - preds_all$pred_xco2)^2) / sum((preds_all$xco2 - mean(preds_all$xco2))^2)

message(glue("Full-data RMSE: {round(rmse_full,3)}; MAE: {round(mae_full,3)}; R2: {round(r2_full,3)}"))

# Save final summary
final_summary <- tibble::tibble(
  timestamp = Sys.time(),
  n_total = n_rows,
  model = "ranger",
  n_train = nrow(df_train),
  rmse_test = rmse_test,
  mae_test = mae_test,
  r2_test = r2_test,
  rmse_full = rmse_full,
  mae_full = mae_full,
  r2_full = r2_full
)
write_csv(final_summary, file.path(OUTPUT_DIR, "ranger_final_summary.csv"))
message("Saved final summary to ", file.path(OUTPUT_DIR, "ranger_final_summary.csv"))

# -----------------------
# 8. Optional: Save subset of predictions (Parquet/feather) for downstream use
# -----------------------
# If you want to keep a compact RDS with predictions for quick analysis:
preds_sample_for_debug <- fread(pred_out_path, nrows = min(200000, n_rows))
saveRDS(preds_sample_for_debug, file.path(OUTPUT_DIR, "predictions_sample_rds.rds"))
message("Saved a debug sample of predictions to predictions_sample_rds.rds")

message("\nDone. All outputs are in: ", OUTPUT_DIR)
# =============================================================================