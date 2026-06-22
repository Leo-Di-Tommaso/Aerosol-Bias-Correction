library(ncdf4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mgcv)   # for spline smoothing via gam

output_dir <- "./results_bias_analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- ALL SCENES TO PROCESS ----
scene_files <- c(
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t349/20250516/c0BC915D0/p3838/MSAT_L2_c0BC915D0_p3838_v04001006_20250516T212530Z_212602Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t349/20250119/c069C15D0/p4210/MSAT_L2_c069C15D0_p4210_v04003001_20250119T213734Z_213806Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t348/20250220/c073C15C0/p3678/MSAT_L2_c073C15C0_p3678_v04001006_20250220T132419Z_132451Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t348/20250430/c0A9D15C0/p3665/MSAT_L2_c0A9D15C0_p3665_v04001006_20250430T132730Z_132802Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t361/20250320/c09041690/p1403/MSAT_L2_c09041690_p1403_v03004000_20250320T063412Z_063444Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t380/20250319/c08F617C0/p1524/MSAT_L2_c08F617C0_p1524_v03005001_20250319T080703Z_080735Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t381/20250426/c0A4217D0/p1618/MSAT_L2_c0A4217D0_p1618_v03005001_20250426T144316Z_144348Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t25/20241119/c04350190/p898/MSAT_L2_c04350190_p898_v03002000_20241119T122526Z_122558Z.nc",
  "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results/t26/20241124/c049501A0/p748/MSAT_L2_c049501A0_p748_v03002000_20241124T143538Z_143610Z.nc"
)

# Short ID from filename, e.g. "c0BC915D0_p3838_20250516" for titles/filenames
scene_id <- function(path) {
  base <- tools::file_path_sans_ext(basename(path))      # MSAT_L2_c..._p..._v..._DATE_TIME
  parts <- strsplit(base, "_")[[1]]
  coll <- parts[3]; pass <- parts[4]
  date <- sub("T.*", "", parts[6])                        # 20250516 from 20250516T212530Z
  paste(coll, pass, date, sep = "_")
}

# ---- PER-SCENE PROCESSING ----
process_scene <- function(scene_file) {
  sid <- scene_id(scene_file)
  message("Processing ", sid, " ...")

  nc <- nc_open(scene_file)
  on.exit(nc_close(nc), add = TRUE)   # always close, even on error

  get_flat <- function(name) {
    v <- ncvar_get(nc, name)
    fill <- ncatt_get(nc, name, "_FillValue")
    if (fill$hasatt) v[v == fill$value] <- NA
    as.vector(v)
  }

  ch4_vcd <- get_flat("co2proxy_fit_diagnostics/retrieved_ch4_vcd")
  co2_vcd <- get_flat("co2proxy_fit_diagnostics/retrieved_co2_vcd")
  albedo  <- get_flat("co2proxy_fit_diagnostics/retrieved_albedo_1606nm")
  qflag   <- get_flat("product_co2proxy/main_quality_flag")  # 0 = good

  # ---- CLOUD/QUALITY FILTER FIRST ----
  keep <- qflag == 0 & is.finite(albedo) & is.finite(ch4_vcd) & is.finite(co2_vcd)
  if (sum(keep) < 50) {
    warning("Scene ", sid, " has only ", sum(keep), " good pixels; skipping.")
    return(invisible(NULL))
  }
  albedo  <- albedo[keep]
  ch4_vcd <- ch4_vcd[keep]
  co2_vcd <- co2_vcd[keep]

  # ---- NORMALIZED PERCENT BIAS ----
  pct_bias <- function(x) (x - mean(x)) / mean(x) * 100
  df <- data.frame(
    albedo = albedo,
    CH4 = pct_bias(ch4_vcd),
    CO2 = pct_bias(co2_vcd)
  ) %>%
    pivot_longer(c(CH4, CO2), names_to = "gas", values_to = "bias_pct")

  # ---- PLOT ----
  p <- ggplot(df, aes(albedo, bias_pct, color = gas, fill = gas)) +
    geom_point(alpha = 0.15, size = 0.6) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE, linewidth = 1) +
    labs(title = paste0("Normalized % Bias vs Retrieved Albedo (1606 nm) — ", sid),
         subtitle = sprintf("N = %d good pixels | retrieved VCDs", nrow(df) / 2),
         x = "Retrieved albedo (CO2 band, 1606 nm)", y = "Normalized bias (%)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  out_path <- file.path(output_dir, paste0("bias_vs_albedo_", sid, ".png"))
  ggsave(out_path, p, width = 9, height = 6, dpi = 200)
  message("  saved ", out_path)
  invisible(p)
}

# ---- LOOP OVER ALL SCENES (one bad file won't stop the rest) ----
for (f in scene_files) {
  tryCatch(process_scene(f),
           error = function(e) message("FAILED ", basename(f), ": ", conditionMessage(e)))
}

message("Done. Plots are in ", normalizePath(output_dir))