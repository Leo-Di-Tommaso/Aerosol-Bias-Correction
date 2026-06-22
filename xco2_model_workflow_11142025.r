# ===============================
# MethaneSAT XCO2 — MLR (AICc Best-Subset) + Optional BIC Interaction Refinement + PLS
# With Selection Diagnostics + Pooled Model (scene-wise offset)
# + Model-Anomaly Adjustment of XCO2 and Checks
# Author: Eleanor Walker, Harvard University
# Date of last update: 2024-11-14
# ===============================

# ---- Libraries ----
library(ncdf4)
library(MASS)        # stepAIC
library(dplyr)
library(tidyr)
library(readr)
library(writexl)
library(splines)
library(car)         # vif
library(pls)         # plsr
library(broom)       # tidy/glance helpers
library(MuMIn)       # dredge, AICc
library(withr)
library(stringr)

# ---- Configuration ----
VIF_THRESHOLD <- 5
ALLOW_INTERACTION_REFINE <- TRUE     # set FALSE to skip BIC refinement
TOP_N_DREDGE <- 10                   # top-N models to log
MAX_SAVE_ROWS_ADJUSTED <- 200000L    # cap for writing adjusted time series

# ---- Inputs / Targets ----
cal_val_target_ids <- c(
  369, 363, 350, 377, 378, 23, 349, 356, 354, 26, 27, 24, 25,
  381, 379, 360, 368, 352, 362, 359, 348, 365, 355, 361, 375, 351, 371, 373
)

# ---- Helpers ----

tryGet <- function(nc, name, default = NA_real_) {
  out <- try(ncvar_get(nc, name), silent = TRUE)
  if (inherits(out, "try-error")) return(default)
  out
}

aggregate_matrix <- function(m, row_block_size = 16, col_block_size = 4) {
  if (!is.matrix(m) || all(is.na(m))) return(m)
  nr <- ceiling(nrow(m) / row_block_size)
  nc <- ceiling(ncol(m) / col_block_size)
  result <- matrix(NA_real_, nrow = nr, ncol = nc)
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      rows <- ((i - 1) * row_block_size + 1):min(i * row_block_size, nrow(m))
      cols <- ((j - 1) * col_block_size + 1):min(j * col_block_size, ncol(m))
      vals <- m[rows, cols, drop = FALSE]
      result[i, j] <- mean(vals, na.rm = TRUE)
    }
  }
  result
}

z <- function(x) as.numeric(scale(x))

# Compute VIFs among predictors only (no response). Returns named numeric vector.
vifs_named <- function(df, xvars) {
  xvars <- xvars[!is.na(xvars) & nzchar(xvars)]
  if (length(xvars) == 0) return(setNames(numeric(0), character(0)))
  if (length(xvars) == 1) return(setNames(1.0, xvars))
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

# Max VIF convenience
max_vif_num <- function(df, xvars) {
  v <- vifs_named(df, xvars)
  if (length(v) == 0) return(NA_real_)
  suppressWarnings(max(v, na.rm = TRUE))
}

# VIF-based pruning (uses same VIF definition as logging). Returns list(keep, trace)
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
      step = step,
      drop = drop_var,
      vif_max = vmax,
      current_predictors = paste(names(vifs), collapse = ", "),
      stringsAsFactors = FALSE
    )
    xvars <- setdiff(xvars, drop_var)
  }
  trace <- if (length(drops)) dplyr::bind_rows(drops) else
    data.frame(step = integer(), drop = character(), vif_max = double(),
               current_predictors = character(), stringsAsFactors = FALSE)
  list(keep = xvars, trace = trace)
}

# small helpers for checks
safe_var <- function(x) stats::var(x, na.rm = TRUE)
safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep]))
}

# ---- Results containers ----
mlr_results <- tibble(
  target_id  = integer(),
  scene_date = character(),
  file       = character(),
  n          = integer(),
  kept_predictors = character(),
  aicc       = double(),
  adj_r2     = double(),
  n_vars_input = integer(),
  n_vars_nzv   = integer(),
  n_vars_vif   = integer(),
  n_vars_best  = integer(),
  n_terms_selected = integer(),
  vif_max_before  = double(),
  vif_max_after   = double(),
  penalty_used    = character(),
  bic_delta       = double(),
  best_by         = character(),
  vif_drop_trace  = character()
)

dredge_results <- tibble(
  target_id   = integer(),
  scene_date  = character(),
  file        = character(),
  rank        = integer(),
  delta_aicc  = double(),
  aicc        = double(),
  n_terms     = integer(),
  adj_r2      = double(),
  terms       = character()
)

pls_results <- tibble(
  target_id  = integer(),
  scene_date = character(),
  file       = character(),
  n          = integer(),
  predictors = character(),
  ncomp      = integer(),
  rmsep_min  = double(),
  r2_val     = double()
)

rmsep_curves <- tibble(
  target_id  = integer(),
  scene_date = character(),
  file       = character(),
  comp       = integer(),
  rmsep      = double(),
  r2_cv      = double()
)

# NEW: per-scene checks after adjustment
adjust_checks <- tibble(
  target_id       = integer(),
  scene_date      = character(),
  file            = character(),
  n               = integer(),
  var_xco2        = double(),
  var_xco2_adj    = double(),
  var_reduction_pct = double(),
  max_abs_cor_before = double(),
  max_abs_cor_after  = double()
)

# pooled storage for pooled model
pooled_store <- list()

# ---- Main loop ----
scene_specific_id_old <- ""
for (scene_id in cal_val_target_ids) {
  message(sprintf("Processing target: t%d", scene_id))
  base_dir <- paste0("/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t", scene_id)
  l2_pp_files <- list.files(base_dir, recursive = TRUE, pattern = "\\.nc$", full.names = TRUE)
  if (length(l2_pp_files) == 0) next

  for (l2_pp in l2_pp_files) {
    # Skip duplicates by scene-specific folder element
    parts <- strsplit(l2_pp, "/")[[1]]
    scene_date <- parts[7]
    scene_specific_id <- parts[8]
    if (identical(scene_specific_id, scene_specific_id_old)) next
    scene_specific_id_old <- scene_specific_id

    # ---- Load NetCDF ----
    nc <- nc_open(l2_pp)

    qflag_o2 <- tryGet(nc, "product_o2proxy/main_quality_flag")
    qflag_c2 <- tryGet(nc, "product_co2proxy/main_quality_flag")

    xco2 <- tryGet(nc, "product_o2proxy/xco2_o2proxy_vcdcor")
    dp   <- tryGet(nc, "o2dp_fit_diagnostics/delta_pressure")

    alb1286 <- tryGet(nc, "o2dp_fit_diagnostics/retrieved_albedo_1286nm")
    alb1606 <- tryGet(nc, "co2proxy_fit_diagnostics/retrieved_albedo_1606nm")
    psfc    <- tryGet(nc, "co2proxy_fit_diagnostics/retrieved_surface_pressure")
    orog    <- tryGet(nc, "geolocation/terrain_height")

    sza <- tryGet(nc, "geolocation/sza")
    vza <- tryGet(nc, "geolocation/vza")
    saa <- tryGet(nc, "geolocation/saa")
    vaa <- tryGet(nc, "geolocation/vaa")

    nc_close(nc)

    # ---- QC ----
    xco2[xco2 > 1e36] <- NA
    if (is.matrix(qflag_o2)) {
      xco2[qflag_o2 == 1]    <- NA
      dp[qflag_o2 == 1]      <- NA
      alb1286[qflag_o2 == 1] <- NA
    }

    pct_good <- mean(is.finite(xco2)) * 100
    if (!is.finite(pct_good) || pct_good < 50) {
      message("  <50% valid XCO2; skipping file: ", l2_pp)
      next
    }

    # ---- Aggregate / bin ----
    xco2 <- aggregate_matrix(xco2)
    dp   <- aggregate_matrix(dp)
    alb1286 <- aggregate_matrix(alb1286)
    psfc <- aggregate_matrix(psfc)
    orog <- aggregate_matrix(orog)
    alb1606 <- aggregate_matrix(alb1606)
    if (!all(is.na(sza))) sza <- aggregate_matrix(sza)
    if (!all(is.na(vza))) vza <- aggregate_matrix(vza)
    if (!all(is.na(saa))) saa <- aggregate_matrix(saa)
    if (!all(is.na(vaa))) vaa <- aggregate_matrix(vaa)

    # ---- Build dataframe (vectorize) ----
    df <- data.frame(
      xco2 = as.vector(xco2),
      dp   = as.vector(dp),
      alb1286 = as.vector(alb1286),
      psfc    = as.vector(psfc),
      orog    = as.vector(orog),
      alb1606 = as.vector(alb1606),
      sza     = as.vector(if (all(is.na(sza))) matrix(NA, nrow(xco2), ncol(xco2)) else sza),
      vza     = as.vector(if (all(is.na(vza))) matrix(NA, nrow(xco2), ncol(xco2)) else vza)
    )

    # Primary angle (deg)
    if (!all(is.na(sza)) && !all(is.na(vza)) && !all(is.na(saa)) && !all(is.na(vaa))) {
      sza_v <- as.vector(sza); vza_v <- as.vector(vza)
      saa_v <- as.vector(saa); vaa_v <- as.vector(vaa)
      primary_angle <- acos(
        pmax(-1, pmin(1,
          cos(sza_v*pi/180) * cos(vza_v*pi/180) +
          sin(sza_v*pi/180) * sin(vza_v*pi/180) * cos((saa_v - vaa_v)*pi/180)
        ))
      ) * 180/pi
    } else {
      primary_angle <- rep(NA_real_, length(df$xco2))
    }
    df$primary_angle <- primary_angle

    # ---- Complete cases ----
    vars_req <- c("dp","alb1286","alb1606","psfc","orog","sza","vza","primary_angle")
    vars_req <- vars_req[vars_req %in% names(df)]
    mask_cc <- complete.cases(df[, c("xco2", vars_req)])
    if (sum(mask_cc) < 50) {
      message("  <50 complete cases across requested vars; skipping file: ", l2_pp)
      next
    }
    df_cc <- df[mask_cc, , drop = FALSE]
    n_cc  <- nrow(df_cc)

    # Hook: store rows for pooled model (raw scale; label scene)
    pooled_piece <- df_cc[, c("xco2", vars_req), drop = FALSE]
    pooled_piece$scene <- paste0("t", scene_id)
    pooled_store[[length(pooled_store) + 1L]] <- pooled_piece

    # Zero-variance pruning (original scale)
    n_vars_input <- length(vars_req)
    nzv <- vapply(vars_req, function(v) sd(df_cc[[v]], na.rm = TRUE) > 0, logical(1))
    vars_req <- vars_req[nzv]
    n_vars_nzv <- length(vars_req)
    if (n_vars_nzv == 0) { message("  All predictors zero-variance; skipping file: ", l2_pp); next }

    # Z-score survivors
    df_cc_z <- df_cc
    for (v in vars_req) df_cc_z[[v]] <- z(df_cc_z[[v]])

    # VIF pruning
    vif_before <- max_vif_num(df_cc_z, vars_req)
    vp <- vif_prune(df_cc_z, y = "xco2", xvars = vars_req, vif_thresh = VIF_THRESHOLD)
    kept_main  <- vp$keep
    vif_trace  <- vp$trace
    n_vars_vif <- length(kept_main)
    if (n_vars_vif == 0) { message("  All predictors removed by VIF; skipping file: ", l2_pp); next }
    vif_after  <- max_vif_num(df_cc_z, kept_main)

    # Best-subset over mains (exclude null)
    df_sub <- df_cc_z[, c("xco2", kept_main), drop = FALSE]
    global_main <- lm(xco2 ~ ., data = df_sub)
    dd <- with_options(list(na.action = "na.fail"), {
      MuMIn::dredge(global_main, rank = "AICc", m.lim = c(1, n_vars_vif), trace = FALSE)
    })
    if (!nrow(dd)) { message("  Dredge returned 0 models; skipping file: ", l2_pp); next }

    # pick best by AICc; ties -> highest adj-R2
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

    # Surviving mains
    terms_sel <- attr(terms(mlr_sel), "term.labels")
    main_in_model <- intersect(kept_main, terms_sel)

    # Stats
    aicc_val <- MuMIn::AICc(mlr_sel)
    adj_r2   <- summary(mlr_sel)$adj.r.squared

    # Log Top-N dredge
    N <- min(TOP_N_DREDGE, nrow(dd))
    models_N <- MuMIn::get.models(dd, subset = 1:N)
    for (i in seq_len(N)) {
      m_i <- models_N[[i]]
      trm <- attr(terms(m_i), "term.labels")
      dredge_results <- add_row(
        dredge_results,
        target_id  = scene_id,
        scene_date = scene_date,
        file       = l2_pp,
        rank       = i,
        delta_aicc = as.numeric(dd$delta[i]),
        aicc       = as.numeric(dd$AICc[i]),
        n_terms    = length(trm),
        adj_r2     = summary(m_i)$adj.r.squared,
        terms      = paste(trm, collapse = ", ")
      )
    }

    # Optional interaction refinement
    n_vars_best <- length(main_in_model)
    n_terms_selected <- length(terms_sel)
    penalty_used <- "AICc (mains)"
    bic_delta <- NA_real_
    if (ALLOW_INTERACTION_REFINE && length(main_in_model) >= 2) {
      rhs_mains <- paste(main_in_model, collapse = " + ")
      form_main <- as.formula(paste("xco2 ~", rhs_mains))
      form_up   <- as.formula(paste("xco2 ~ (", rhs_mains, ")^2"))
      start_mod <- lm(form_main, data = df_cc_z)
      k_pen <- log(n_cc)
      mlr_refined <- suppressWarnings(MASS::stepAIC(
        start_mod,
        scope = list(lower = ~ 1, upper = form_up),
        direction = "both",
        k = k_pen,
        trace = FALSE
      ))
      bic_start   <- AIC(start_mod,  k = k_pen)
      bic_refined <- AIC(mlr_refined, k = k_pen)
      adopt <- FALSE
      if (bic_refined + tol < bic_start) adopt <- TRUE
      else if (abs(bic_refined - bic_start) <= tol) {
        adj_start   <- summary(start_mod)$adj.r.squared
        adj_refined <- summary(mlr_refined)$adj.r.squared
        adopt <- is.finite(adj_refined) && (adj_refined > adj_start + tol)
      }
      if (adopt) {
        mlr_sel <- mlr_refined
        terms_sel <- attr(terms(mlr_sel), "term.labels")
        mains_present <- intersect(main_in_model, terms_sel)
        inter_terms <- grep(":", terms_sel, value = TRUE)
        if (length(inter_terms) > 0) {
          comps <- unique(unlist(strsplit(inter_terms, ":", fixed = TRUE)))
          mains_present <- union(mains_present, intersect(comps, main_in_model))
        }
        main_in_model <- unique(mains_present)
        aicc_val <- MuMIn::AICc(mlr_sel)
        adj_r2   <- summary(mlr_sel)$adj.r.squared
        n_terms_selected <- length(terms_sel)
        penalty_used <- "AICc (mains) + BIC refine"
        bic_delta <- as.numeric(bic_refined - bic_start)
        best_main_choice <- paste0(best_main_choice, " → BIC",
                                   if (abs(bic_refined - bic_start) <= tol) "+adjR2" else "")
      }
    }

    # ---- NEW: Adjust XCO2 by model anomaly & checks (per scene/file) ----
    # Predict on the SAME df_cc_z used to fit (has z-scored predictors)
    x_pred <- as.numeric(predict(mlr_sel, newdata = df_cc_z))
    model_anom <- x_pred - mean(x_pred, na.rm = TRUE)
    xco2_adj <- df_cc_z$xco2 - model_anom

    # simple diagnostics: variance reduction and correlation drops
    var_before <- safe_var(df_cc_z$xco2)
    var_after  <- safe_var(xco2_adj)
    red_pct    <- if (is.finite(var_before) && var_before > 0) 100*(var_before - var_after)/var_before else NA_real_

    # correlations w.r.t. selected mains (marginal, abs)
    if (length(main_in_model) > 0) {
      cors_before <- sapply(main_in_model, function(v) safe_cor(df_cc_z$xco2, df_cc_z[[v]]))
      cors_after  <- sapply(main_in_model, function(v) safe_cor(xco2_adj,     df_cc_z[[v]]))
      max_abs_before <- suppressWarnings(max(abs(cors_before), na.rm = TRUE))
      max_abs_after  <- suppressWarnings(max(abs(cors_after),  na.rm = TRUE))
    } else {
      max_abs_before <- NA_real_
      max_abs_after  <- NA_real_
    }

    adjust_checks <- add_row(
      adjust_checks,
      target_id = scene_id,
      scene_date = scene_date,
      file = l2_pp,
      n = n_cc,
      var_xco2 = var_before,
      var_xco2_adj = var_after,
      var_reduction_pct = red_pct,
      max_abs_cor_before = max_abs_before,
      max_abs_cor_after  = max_abs_after
    )

    # Optionally save the adjusted series (capped for safety)
    out_dir <- paste0("/home/ext_ewalker_g_harvard_edu/results/xco2_model_", format(Sys.Date(), "%Y%m%d"))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    if (length(xco2_adj) <= MAX_SAVE_ROWS_ADJUSTED) {
      out_adj <- file.path(out_dir, paste0("xco2_adjusted_", basename(l2_pp), ".csv"))
      readr::write_csv(
        tibble(xco2 = df_cc_z$xco2, xco2_model = x_pred, xco2_model_anom = model_anom, xco2_adjusted = xco2_adj),
        out_adj
      )
    }

    # Log MLR selection diagnostics
    mlr_results <- add_row(
      mlr_results,
      target_id = scene_id,
      scene_date = scene_date,
      file = l2_pp,
      n = n_cc,
      kept_predictors = paste(main_in_model, collapse = ", "),
      aicc = aicc_val,
      adj_r2 = adj_r2,
      n_vars_input = n_vars_input,
      n_vars_nzv   = n_vars_nzv,
      n_vars_vif   = n_vars_vif,
      n_vars_best  = n_vars_best,
      n_terms_selected = n_terms_selected,
      vif_max_before  = vif_before,
      vif_max_after   = vif_after,
      penalty_used    = penalty_used,
      bic_delta       = bic_delta,
      best_by         = best_main_choice,
      vif_drop_trace  = if (nrow(vif_trace)) {
        paste(sprintf("#%d drop=%s (VIF=%.2f) | pool={%s}",
                      vif_trace$step, vif_trace$drop, vif_trace$vif_max,
                      vif_trace$current_predictors), collapse = "  ||  ")
      } else ""
    )

    # ---- PLS on remaining main-effect variables ----
    if (length(main_in_model) >= 1) {
      pls_df <- df_cc_z[, c("xco2", main_in_model), drop = FALSE]
      max_comp <- max(1, min(10, length(main_in_model), nrow(pls_df) - 1))

      set.seed(42)
      pls_fit <- plsr(
        xco2 ~ .,
        data = pls_df,
        ncomp = max_comp,
        validation = "CV",
        segments = min(10, nrow(pls_df)),
        center = TRUE, scale = FALSE
      )

      pred_cv <- pls_fit$validation$pred  # [n, 1, K]
      ncomp_avail <- dim(pred_cv)[3]
      y <- pls_df$xco2
      y_pred <- sapply(seq_len(ncomp_avail), function(k) pred_cv[, 1, k])

      press_vals <- colSums((y - y_pred)^2)
      rmsep_vals <- sqrt(press_vals / length(y))
      tss <- sum((y - mean(y))^2)
      r2_curve <- pmax(0, 1 - press_vals / tss)

      best_k <- which.min(rmsep_vals)
      rm_min <- rmsep_vals[best_k]
      r2_cv  <- r2_curve[best_k]

      pls_results <- add_row(
        pls_results,
        target_id = scene_id,
        scene_date = scene_date,
        file = l2_pp,
        n = n_cc,
        predictors = paste(main_in_model, collapse = ", "),
        ncomp = best_k,
        rmsep_min = rm_min,
        r2_val = r2_cv
      )

      rmsep_curves <- bind_rows(
        rmsep_curves,
        tibble(
          target_id  = scene_id,
          scene_date = scene_date,
          file       = l2_pp,
          comp       = seq_len(ncomp_avail),
          rmsep      = rmsep_vals,
          r2_cv      = r2_curve
        )
      )
    }

    # Console breadcrumb
    message(sprintf(
      "  %s | n=%d | vars: input=%d, nzv=%d, vif=%d, best=%d | terms=%d | AICc=%.1f | VIFmax before/after=%.2f/%.2f | Var↓=%.1f%%",
      basename(l2_pp), n_cc, n_vars_input, n_vars_nzv, n_vars_vif, n_vars_best,
      n_terms_selected, aicc_val, vif_before, vif_after, red_pct
    ))
  }
}

# ===============================
# POOLED MODEL over all scenes (with scene-wise offset)
# ===============================
if (length(pooled_store) > 0) {
  pooled_df <- dplyr::bind_rows(pooled_store)
  pooled_df$scene <- factor(pooled_df$scene)  # per-scene intercepts

  vars_all <- c("dp","alb1286","alb1606","psfc","orog","sza","vza","primary_angle")
  vars_pool <- intersect(vars_all, names(pooled_df))

  mask_cc <- stats::complete.cases(pooled_df[, c("xco2", vars_pool)])
  pooled_cc <- pooled_df[mask_cc, , drop = FALSE]
  n_pool <- nrow(pooled_cc)
  n_scenes <- nlevels(pooled_cc$scene)

  if (n_pool >= 200 && length(vars_pool) > 0) {

    # zero-variance prune (original scale)
    nzv <- vapply(vars_pool, function(v) stats::sd(pooled_cc[[v]], na.rm = TRUE) > 0, logical(1))
    vars_pool <- vars_pool[nzv]
    if (length(vars_pool) > 0) {
      # z-score survivors (not xco2, not scene)
      pooled_cc_z <- pooled_cc
      for (v in vars_pool) pooled_cc_z[[v]] <- z(pooled_cc_z[[v]])

      # VIF pruning among predictors (scene excluded)
      vif_before_pool <- max_vif_num(pooled_cc_z, vars_pool)
      vp_pool <- vif_prune(pooled_cc_z, y = "xco2", xvars = vars_pool, vif_thresh = VIF_THRESHOLD)
      kept_main_pool <- vp_pool$keep
      vif_trace_pool <- vp_pool$trace

      if (length(kept_main_pool) > 0) {
        vif_after_pool <- max_vif_num(pooled_cc_z, kept_main_pool)

        # Best-subset mains with scene fixed; exclude null
        df_sub <- pooled_cc_z[, c("xco2", "scene", kept_main_pool), drop = FALSE]
        global_main <- lm(xco2 ~ scene + ., data = df_sub)
        dd_pool <- with_options(list(na.action = "na.fail"), {
          MuMIn::dredge(global_main, rank = "AICc", fixed = ~ scene,
                        m.lim = c(1, length(kept_main_pool)), trace = FALSE)
        })
        if (!nrow(dd_pool)) stop("Pooled dredge returned 0 models; check inputs.")

        aicc_min <- min(dd_pool$AICc); tol <- 1e-9
        cand_rows <- which(abs(dd_pool$AICc - aicc_min) <= tol)
        best_main_choice_pool <- "AICc"
        if (length(cand_rows) == 1) {
          mlr_sel_pool <- MuMIn::get.models(dd_pool, subset = cand_rows)[[1]]
        } else {
          cmods <- MuMIn::get.models(dd_pool, subset = cand_rows)
          cadjr <- sapply(cmods, function(m) summary(m)$adj.r.squared)
          best_i <- cand_rows[ which.max(cadjr) ]
          mlr_sel_pool <- MuMIn::get.models(dd_pool, subset = best_i)[[1]]
          best_main_choice_pool <- "AICc+adjR2"
        }

        terms_sel_pool <- attr(terms(mlr_sel_pool), "term.labels")
        main_in_model_pool <- intersect(kept_main_pool, terms_sel_pool)
        aicc_pool <- MuMIn::AICc(mlr_sel_pool)
        adjr2_pool <- summary(mlr_sel_pool)$adj.r.squared

        # Optional interaction refine (scene fixed)
        penalty_used_pool <- "AICc (mains + scene)"
        bic_delta_pool <- NA_real_
        n_terms_selected_pool <- length(terms_sel_pool)
        if (ALLOW_INTERACTION_REFINE && length(main_in_model_pool) >= 2) {
          rhs_mains <- paste(main_in_model_pool, collapse = " + ")
          form_main <- as.formula(paste("xco2 ~ scene +", rhs_mains))
          form_up   <- as.formula(paste("xco2 ~ scene + (", rhs_mains, ")^2"))
          start_mod <- lm(form_main, data = pooled_cc_z)
          k_pen <- log(n_pool)
          mlr_refined <- suppressWarnings(MASS::stepAIC(
            start_mod, scope = list(lower = ~ scene, upper = form_up),
            direction = "both", k = k_pen, trace = FALSE
          ))
          bic_start   <- AIC(start_mod,  k = k_pen)
          bic_refined <- AIC(mlr_refined, k = k_pen)
          adopt <- FALSE
          if (bic_refined + tol < bic_start) adopt <- TRUE
          else if (abs(bic_refined - bic_start) <= tol) {
            adj_start   <- summary(start_mod)$adj.r.squared
            adj_refined <- summary(mlr_refined)$adj.r.squared
            adopt <- is.finite(adj_refined) && (adj_refined > adj_start + tol)
          }
          if (adopt) {
            mlr_sel_pool <- mlr_refined
            terms_sel_pool <- attr(terms(mlr_sel_pool), "term.labels")
            mains_present <- intersect(main_in_model_pool, terms_sel_pool)
            inter_terms <- grep(":", terms_sel_pool, value = TRUE)
            if (length(inter_terms) > 0) {
              comps <- unique(unlist(strsplit(inter_terms, ":", fixed = TRUE)))
              mains_present <- union(mains_present, intersect(comps, main_in_model_pool))
            }
            main_in_model_pool <- unique(mains_present)
            aicc_pool <- MuMIn::AICc(mlr_sel_pool)
            adjr2_pool <- summary(mlr_sel_pool)$adj.r.squared
            n_terms_selected_pool <- length(terms_sel_pool)
            penalty_used_pool <- "AICc (mains + scene) + BIC refine"
            bic_delta_pool <- as.numeric(bic_refined - bic_start)
            best_main_choice_pool <- paste0(best_main_choice_pool, " → BIC",
                                            if (abs(bic_refined - bic_start) <= tol) "+adjR2" else "")
          }
        }

        # ---- NEW: Adjust XCO2 by model anomaly & checks (POOLED) ----
        # By default, subtract GLOBAL pooled mean of predictions:
        # If you prefer per-scene centering, replace mean(x_pred_pool) with ave(x_pred_pool, pooled_cc_z$scene).
        x_pred_pool <- as.numeric(predict(mlr_sel_pool, newdata = pooled_cc_z))
        model_anom_pool <- x_pred_pool - mean(x_pred_pool, na.rm = TRUE)
        xco2_adj_pool <- pooled_cc_z$xco2 - model_anom_pool

        var_before_pool <- safe_var(pooled_cc_z$xco2)
        var_after_pool  <- safe_var(xco2_adj_pool)
        red_pct_pool    <- if (is.finite(var_before_pool) && var_before_pool > 0)
          100*(var_before_pool - var_after_pool)/var_before_pool else NA_real_

        if (length(main_in_model_pool) > 0) {
          cors_b <- sapply(main_in_model_pool, function(v) safe_cor(pooled_cc_z$xco2, pooled_cc_z[[v]]))
          cors_a <- sapply(main_in_model_pool, function(v) safe_cor(xco2_adj_pool,  pooled_cc_z[[v]]))
          max_abs_b <- suppressWarnings(max(abs(cors_b), na.rm = TRUE))
          max_abs_a <- suppressWarnings(max(abs(cors_a), na.rm = TRUE))
        } else {
          max_abs_b <- NA_real_; max_abs_a <- NA_real_
        }

        # Save pooled adjusted (capped)
        out_dir <- paste0("/home/ext_ewalker_g_harvard_edu/results/xco2_model_", format(Sys.Date(), "%Y%m%d"))
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
        if (length(xco2_adj_pool) <= MAX_SAVE_ROWS_ADJUSTED) {
          readr::write_csv(
            tibble(scene = pooled_cc_z$scene,
                   xco2 = pooled_cc_z$xco2,
                   xco2_model = x_pred_pool,
                   xco2_model_anom = model_anom_pool,
                   xco2_adjusted = xco2_adj_pool),
            file.path(out_dir, "xco2_pooled_adjusted.csv")
          )
        }

        # Log pooled summaries
        mlr_pool_row <- dplyr::tibble(
          n_rows = n_pool,
          n_scenes = n_scenes,
          kept_predictors = paste(main_in_model_pool, collapse = ", "),
          aicc = aicc_pool,
          adj_r2 = adjr2_pool,
          n_terms_selected = n_terms_selected_pool,
          penalty_used = penalty_used_pool,
          best_by = best_main_choice_pool,
          vif_max_before = vif_before_pool,
          vif_max_after  = vif_after_pool,
          bic_delta = bic_delta_pool,
          var_xco2 = var_before_pool,
          var_xco2_adj = var_after_pool,
          var_reduction_pct = red_pct_pool,
          max_abs_cor_before = max_abs_b,
          max_abs_cor_after  = max_abs_a,
          vif_drop_trace = if (nrow(vif_trace_pool)) {
            paste(sprintf("#%d drop=%s (VIF=%.2f) | pool={%s}",
                          vif_trace_pool$step, vif_trace_pool$drop, vif_trace_pool$vif_max,
                          vif_trace_pool$current_predictors), collapse = "  ||  ")
          } else ""
        )

        # Pooled top-N dredge rows (predefine columns for add_row safety)
        Np <- min(TOP_N_DREDGE, nrow(dd_pool))
        models_N <- MuMIn::get.models(dd_pool, subset = 1:Np)
        dredge_pooled <- dplyr::tibble(
          rank      = integer(),
          delta_aicc= double(),
          aicc      = double(),
          n_terms   = integer(),
          adj_r2    = double(),
          terms     = character()
        )
        for (i in seq_len(Np)) {
          m_i <- models_N[[i]]
          trm <- attr(terms(m_i), "term.labels")
          dredge_pooled <- dplyr::add_row(
            dredge_pooled,
            rank       = i,
            delta_aicc = as.numeric(dd_pool$delta[i]),
            aicc       = as.numeric(dd_pool$AICc[i]),
            n_terms    = length(setdiff(trm, "scene")),
            adj_r2     = summary(m_i)$adj.r.squared,
            terms      = paste(setdiff(trm, "scene"), collapse = ", ")
          )
        }

        # Write pooled files
        readr::write_csv(mlr_pool_row,       file.path(out_dir, "xco2_pooled_mlr_aicc.csv"))
        readr::write_csv(dredge_pooled,      file.path(out_dir, "xco2_pooled_mlr_dredge_top_models.csv"))

        # PLS on mains only (no scene) + write curves
        if (length(main_in_model_pool) >= 1) {
          pls_df_pool <- pooled_cc_z[, c("xco2", main_in_model_pool), drop = FALSE]
          max_comp <- max(1, min(10, length(main_in_model_pool), nrow(pls_df_pool) - 1))
          set.seed(42)
          pls_fit_pool <- pls::plsr(
            xco2 ~ .,
            data = pls_df_pool,
            ncomp = max_comp,
            validation = "CV",
            segments = min(10, nrow(pls_df_pool)),
            center = TRUE, scale = FALSE
          )
          pred_cv <- pls_fit_pool$validation$pred
          ncomp_avail <- dim(pred_cv)[3]
          y  <- pls_df_pool$xco2
          yk <- sapply(seq_len(ncomp_avail), function(k) pred_cv[, 1, k])
          press_vals <- colSums((y - yk)^2)
          rmsep_vals <- sqrt(press_vals / length(y))
          tss <- sum((y - mean(y))^2)
          r2_curve <- pmax(0, 1 - press_vals / tss)
          best_k <- which.min(rmsep_vals)
          rm_min <- rmsep_vals[best_k]
          r2_cv  <- r2_curve[best_k]
          readr::write_csv(
            tibble(n_rows = n_pool, n_scenes = n_scenes,
                   predictors = paste(main_in_model_pool, collapse = ", "),
                   ncomp = best_k, rmsep_min = rm_min, r2_val = r2_cv),
            file.path(out_dir, "xco2_pooled_pls.csv")
          )
          readr::write_csv(
            tibble(comp = seq_len(ncomp_avail), rmsep = rmsep_vals, r2_cv = r2_curve),
            file.path(out_dir, "xco2_pooled_pls_rmsep_curve.csv")
          )
        }

        message(sprintf(
          "Pooled: n=%d rows across %d scenes | mains=%d | AICc=%.1f | adjR2=%.3f | VIFmax pre/post=%.2f/%.2f | Var↓=%.1f%%",
          n_pool, n_scenes, length(main_in_model_pool), aicc_pool, adjr2_pool,
          vif_before_pool, vif_after_pool, red_pct_pool
        ))
      }
    }
  } else {
    message("Pooled: not enough complete rows or no predictors after pruning; skipping pooled model.")
  }
} else {
  message("Pooled: no per-file data collected; skipping pooled model.")
}

# ---- Write per-file outputs ----
mlr_results <- mlr_results %>% arrange(target_id, scene_date, file)
pls_results <- pls_results %>% arrange(target_id, scene_date, file)
rmsep_curves <- rmsep_curves %>% arrange(target_id, scene_date, file, comp)
dredge_results <- dredge_results %>% arrange(target_id, scene_date, file, rank)
adjust_checks  <- adjust_checks  %>% arrange(target_id, scene_date, file)

out_dir <- paste0("/home/ext_ewalker_g_harvard_edu/results/xco2_model_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_xlsx <- file.path(out_dir, "xco2_mlr_aicc_and_pls_by_scene.xlsx")
out_csv1 <- file.path(out_dir, "xco2_mlr_aicc_by_scene.csv")
out_csv2 <- file.path(out_dir, "xco2_pls_by_scene.csv")
out_csv3 <- file.path(out_dir, "xco2_pls_rmsep_curves_by_scene.csv")
out_csv4 <- file.path(out_dir, "xco2_mlr_dredge_top_models_by_scene.csv")
out_csv5 <- file.path(out_dir, "xco2_adjustment_checks_by_scene.csv")   # NEW

writexl::write_xlsx(
  list(
    "mlr_aicc"         = mlr_results,
    "pls"              = pls_results,
    "pls_rmsep_curve"  = rmsep_curves,
    "mlr_dredge_topN"  = dredge_results,
    "xco2_adjust_checks" = adjust_checks
  ),
  path = out_xlsx
)

readr::write_csv(mlr_results, out_csv1)
readr::write_csv(pls_results, out_csv2)
readr::write_csv(rmsep_curves, out_csv3)
readr::write_csv(dredge_results, out_csv4)
readr::write_csv(adjust_checks, out_csv5)

message("Wrote Excel: ", out_xlsx)
message("Wrote CSVs: ", paste(out_csv1, out_csv2, out_csv3, out_csv4, out_csv5, sep = " | "))
