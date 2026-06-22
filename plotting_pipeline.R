# =============================================================================
# MethaneSAT + VIIRS: Scene Visualization Pipeline (Fixed)
# =============================================================================
# This script:
# 1. Accepts a specific list of MSAT files.
# 2. Downloads/Finds matching VIIRS aerosol data.
# 3. Collocates data.
# 4. Generates diagnostic maps and scatterplots for presentation.
# TO RUN THIS SCRIPT: source ~/earthdata-env/bin/activate

# ---------------------------------------------------------------------------
# 0. Libraries and Constants
# ---------------------------------------------------------------------------
library(terra)
library(ncdf4)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(viridis)
library(MASS)   # for rlm()


raw_file_list <- c(
  # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250516/c0BC915D0/p3838/MSAT_L2_c0BC915D0_p3838_v04001006_20250516T212530Z_212602Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250119/c069C15D0/p4210/MSAT_L2_c069C15D0_p4210_v04003001_20250119T213734Z_213806Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t348/20250220/c073C15C0/p3678/MSAT_L2_c073C15C0_p3678_v04001006_20250220T132419Z_132451Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t108/20250120/c06B506C0/p3188/MSAT_L2_c06B506C0_p3188_v04001003_20250120T071415Z_071447Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t108/20250519/c0C1306C0/p3155/MSAT_L2_c0C1306C0_p3155_v04001003_20250519T071300Z_071332Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t51/20241104/c038C0330/p3866/MSAT_L2_c038C0330_p3866_v04001006_20241104T104958Z_105030Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t363/20250228/c07D816B0/p4002/MSAT_L2_c07D816B0_p4002_v04002001_20250228T204049Z_204114Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t60/20241228/c056403C0/p3860/MSAT_L2_c056403C0_p3860_v04001006_20241228T210922Z_210954Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t60/20240917/c017103C0/p3772/MSAT_L2_c017103C0_p3772_v04001006_20240917T211429Z_211501Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t56/20240928/c01CC0380/p3255/MSAT_L2_c01CC0380_p3255_v04001003_20240928T205946Z_210018Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t56/20240929/c01D50380/p3253/MSAT_L2_c01D50380_p3253_v04001003_20240929T210658Z_210730Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t56/20250320/c09060380/p3221/MSAT_L2_c09060380_p3221_v04001003_20250320T210203Z_210235Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t56/20241223/c05040380/p3231/MSAT_L2_c05040380_p3231_v04001003_20241223T204056Z_204128Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t25/20241007/c021E0190/p4566/MSAT_L2_c021E0190_p4566_v04003001_20241007T122638Z_122710Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t25/20240911/c01420190/p3215/MSAT_L2_c01420190_p3215_v04001003_20240911T122836Z_122908Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t25/20241009/c02300190/p3394/MSAT_L2_c02300190_p3394_v04001003_20241009T124038Z_124110Z.nc",
  # "gs://msat-prod-data-methanesat-level2/level2_results/t25/20241220/c04D90190/p3232/MSAT_L2_c04D90190_p3232_v04001003_20241220T122155Z_122227Z.nc"
    # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250516/c0BC915D0/p3838/MSAT_L2_c0BC915D0_p3838_v04001006_20250516T212530Z_212602Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20241230/c058A15D0/p4175/MSAT_L2_c058A15D0_p4175_v04003001_20241230T212120Z_212152Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250119/c069C15D0/p4210/MSAT_L2_c069C15D0_p4210_v04003001_20250119T213734Z_213806Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t349/20250430/c0A9C15D0/p3846/MSAT_L2_c0A9C15D0_p3846_v04001006_20250430T213320Z_213352Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t348/20250220/c073C15C0/p3678/MSAT_L2_c073C15C0_p3678_v04001006_20250220T132419Z_132451Z.nc",
    # "gs://msat-prod-data-methanesat-level2/level2_results/t348/20250430/c0A9D15C0/p3665/MSAT_L2_c0A9D15C0_p3665_v04001006_20250430T132730Z_132802Z.nc"
    
    
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



OUTPUT_DIR <- "presentation_plots"
VIIRS_DATA_DIR <- "viirs_data2"
PYTHON_SCRIPT_PATH <- "/home/ext_lditommaso_college_harvard_e/leo_data/analyses/xco2_aod/earthdata-search/download_viirs_deep_blue2.py"

# Convert gs:// paths to local mount paths
# Input File List (converted from gs:// to local mount)


MSAT_FILES <- gsub("gs://", "/mnt/gcs/", raw_file_list)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive=TRUE)
if (!dir.exists(VIIRS_DATA_DIR)) dir.create(VIIRS_DATA_DIR, recursive=TRUE)

# ---------------------------------------------------------------------------
# 1. Helper Functions
# ---------------------------------------------------------------------------

get_bbox_and_date <- function(msat_file) {
  nc <- try(nc_open(msat_file), silent=TRUE)
  if (inherits(nc, "try-error")) return(NULL)
  
  lon <- try(ncvar_get(nc, "geolocation/longitude"), silent=TRUE)
  lat <- try(ncvar_get(nc, "geolocation/latitude"), silent=TRUE)
  time_vals <- try(ncvar_get(nc, "geolocation/time"), silent=TRUE)
  nc_close(nc)
  
  if (inherits(lon, "try-error") || inherits(lat, "try-error") || inherits(time_vals, "try-error")) return(NULL)
  
  origin <- as.POSIXct("1970-01-01 00:00:00", tz="UTC")
  mid_time <- origin + median(time_vals, na.rm=TRUE)
  
  list(
    bbox = c(min(lon, na.rm=TRUE), min(lat, na.rm=TRUE), max(lon, na.rm=TRUE), max(lat, na.rm=TRUE)),
    date = as.character(as.Date(mid_time))
  )
}

find_or_download_viirs <- function(msat_file) {
  scene_info <- get_bbox_and_date(msat_file)
  if (is.null(scene_info)) return(NA_character_)
  
  date_obj <- as.Date(scene_info$date)
  year_doy <- format(date_obj, "%Y%j")
  search_pattern <- sprintf("AERDB_L2_VIIRS_(SNPP|NOAA20|NOAA21)\\.A%s.*\\.(nc|nc4)$", year_doy)
  
  local_files <- list.files(VIIRS_DATA_DIR, pattern=search_pattern, full.names=TRUE, recursive=TRUE)
  if (length(local_files) > 0) return(local_files[1])
  
  # Download via Python using explicit python3 call
  cat(sprintf(" -> Downloading VIIRS for date: %s...\n", scene_info$date))
  if (!file.exists(PYTHON_SCRIPT_PATH)) {
    warning("Downloader script not found.")
    return(NA_character_)
  }
  
  cmd_args <- c(
    PYTHON_SCRIPT_PATH,
    "--date", scene_info$date,
    "--outdir", VIIRS_DATA_DIR,
    "--bbox", as.character(scene_info$bbox)
  )
  
  res <- tryCatch({
    system2("python3", args = cmd_args, stdout=TRUE, stderr=TRUE)
  }, error=function(e) return(NULL))
  
  downloaded_file <- list.files(VIIRS_DATA_DIR, pattern=search_pattern, full.names=TRUE, recursive=TRUE)
  if(length(downloaded_file) > 0) return(downloaded_file[1])
  
  return(NA_character_)
}
calc_primary_angle <- function(sza, vza, saa, vaa) {
  dtor <- pi / 180
  rtod <- 180 / pi
  cos_phase <- cos(sza*dtor) * cos(vza*dtor) + 
               sin(sza*dtor) * sin(vza*dtor) * cos((saa - vaa)*dtor)
  acos(pmax(-1, pmin(1, cos_phase))) * rtod
}

extract_scene_data <- function(viirs_file, msat_file, scene_id) {
  # 1. Load VIIRS arrays and rasterize to VIIRS native grid
  if (is.na(viirs_file) || !file.exists(viirs_file)) stop("VIIRS file missing")
  v_nc <- nc_open(viirs_file)
  v_lon <- ncvar_get(v_nc, "Longitude")
  v_lat <- ncvar_get(v_nc, "Latitude")
  v_aod <- ncvar_get(v_nc, "Aerosol_Optical_Thickness_550_Land_Best_Estimate")
  v_ae  <- ncvar_get(v_nc, "Angstrom_Exponent_Land_Best_Estimate")
  nc_close(v_nc)

  v_df <- data.frame(lon = as.vector(v_lon), lat = as.vector(v_lat),
                     AOD550 = as.vector(v_aod), AE = as.vector(v_ae)) %>%
    filter(complete.cases(.))
  if (nrow(v_df) == 0) return(NULL)

  v_vect <- vect(v_df, geom=c("lon","lat"), crs="EPSG:4326")
  # Use VIIRS arrays to build a template raster (preserves scan ordering)
  r_template <- rast(nrows=dim(v_lon)[1], ncols=dim(v_lon)[2],
                     xmin=min(v_lon), xmax=max(v_lon),
                     ymin=min(v_lat), ymax=max(v_lat), crs="EPSG:4326")

  r_viirs <- rasterize(v_vect, r_template, field=c("AOD550", "AE"), fun="mean")

  # 2. Load MSAT and apply quality mask (same logic as in MSAT_VIIRS_model_pipeline.R)
  m_nc <- try(nc_open(msat_file), silent=TRUE)
  if (inherits(m_nc, "try-error")) return(NULL)
  get_var <- function(nc, v) { val <- try(ncvar_get(nc, v), silent=TRUE); if(inherits(val, "try-error")) return(NULL); val }

  m_lon <- get_var(m_nc, "geolocation/longitude")
  m_lat <- get_var(m_nc, "geolocation/latitude")
  xco2 <- get_var(m_nc, "product_o2proxy/xco2_o2proxy_vcdcor")
  dp   <- get_var(m_nc, "o2dp_fit_diagnostics/bias_corrected_delta_pressure")
  if (is.null(dp)) dp <- get_var(m_nc, "o2dp_fit_diagnostics/delta_pressure")
  alb1286 <- get_var(m_nc, "o2dp_fit_diagnostics/retrieved_albedo_1286nm")
  alb1606 <- get_var(m_nc, "co2proxy_fit_diagnostics/retrieved_albedo_1606nm")
  psfc    <- get_var(m_nc, "co2proxy_fit_diagnostics/retrieved_surface_pressure")
  orog    <- get_var(m_nc, "geolocation/terrain_height")
  sza     <- get_var(m_nc, "geolocation/sza")
  vza     <- get_var(m_nc, "geolocation/vza")
  saa     <- get_var(m_nc, "geolocation/saa")
  vaa     <- get_var(m_nc, "geolocation/vaa")
  q_flag  <- get_var(m_nc, "product_o2proxy/main_quality_flag")
  airglow_vcd <- get_var(m_nc, "o2dp_fit_diagnostics/airglow_vcd")
  xco2_prior  <- get_var(m_nc, "apriori_data/co2_profile")
  nc_close(m_nc)
  if (is.null(xco2) || is.null(m_lon)) return(NULL)

  # Apply quality mask (mask==1 indicates bad pixels)
  if (!is.null(q_flag)) {
    mask <- q_flag == 1
    if (!is.null(dp)) dp[mask] <- NA
    if (!is.null(xco2)) xco2[mask] <- NA
    if (!is.null(alb1286)) alb1286[mask] <- NA
    if (!is.null(alb1606)) alb1606[mask] <- NA
    if (!is.null(psfc)) psfc[mask] <- NA
    if (!is.null(orog)) orog[mask] <- NA
    if (!is.null(sza)) sza[mask] <- NA
    if (!is.null(vza)) vza[mask] <- NA
    if (!is.null(airglow_vcd)) airglow_vcd[mask] <- NA
    if (!is.null(xco2_prior)) xco2_prior[mask] <- NA
  }

  # Primary angle (optional) — compute if angles are present
  p_angle <- NULL
  if (!is.null(sza) && !is.null(vza) && !is.null(saa) && !is.null(vaa)) {
    p_angle <- calc_primary_angle(sza, vza, saa, vaa)
  }

  m_df <- data.frame(
    lon = as.vector(m_lon), lat = as.vector(m_lat),
    xco2 = as.vector(xco2), dp = as.vector(dp),
    alb1286 = as.vector(alb1286), alb1606 = as.vector(alb1606),
    psfc = as.vector(psfc), orog = as.vector(orog),
    sza = as.vector(sza), vza = as.vector(vza),
    airglow_vcd = as.vector(airglow_vcd), xco2_prior = as.vector(xco2_prior),
    primary_angle = if (is.null(p_angle)) NA_real_ else as.vector(p_angle)
  )
  # Keep only complete predictor rows
  m_df <- m_df %>% filter(complete.cases(.))
  if (nrow(m_df) == 0) return(NULL)

  # Rasterize MSAT to VIIRS template
  m_vect <- vect(m_df, geom=c("lon","lat"), crs="EPSG:4326")
  m_vars <- c("xco2", "dp", "alb1286", "alb1606", "psfc", "orog", "sza", "vza", "primary_angle")
  r_msat <- rasterize(m_vect, r_template, field=m_vars, fun="mean", na.rm=TRUE)

  # Collocate: crop to overlapping extent and stack
  valid_ext <- intersect(ext(r_viirs), ext(r_msat))
  if (is.null(valid_ext)) return(NULL)
  stack_all <- c(crop(r_viirs, valid_ext), crop(r_msat, valid_ext))

  # Ensure names align
  expected_names <- c("AOD550", "AE", "xco2", "dp", "alb1286", "alb1606", "psfc", "orog", "sza", "vza", "primary_angle")
  if (nlyr(stack_all) == length(expected_names)) names(stack_all) <- expected_names

  df_scene <- as.data.frame(stack_all, xy=TRUE, na.rm=TRUE)
  if (nrow(df_scene) < 50) return(NULL)

  # Compute AOD at 1286nm
  if (!all(c("AOD550", "AE") %in% names(df_scene))) return(NULL)
  df_scene$AOD_1286 <- df_scene$AOD550 * (1286/550)^(-df_scene$AE)
  df_scene$scene_id <- scene_id
  return(df_scene)
}

# ---------------------------------------------------------------------------
# 2. Plotting Helpers
# ---------------------------------------------------------------------------


create_map_pair <- function(df, var1, var2, title1, title2, main_title) {
  p1 <- ggplot(df, aes_string(x="x", y="y", fill=var1)) +
    geom_raster() + scale_fill_viridis_c(option="magma", name=title1) +
    coord_fixed() + labs(x=NULL, y=NULL, title=title1) + theme_minimal() +
    theme(legend.position="bottom", panel.grid=element_blank())
  
  p2 <- ggplot(df, aes_string(x="x", y="y", fill=var2)) +
    geom_raster() + scale_fill_viridis_c(option="viridis", name=title2) +
    coord_fixed() + labs(x=NULL, y=NULL, title=title2) + theme_minimal() +
    theme(legend.position="bottom", panel.grid=element_blank())
  
  grob <- arrangeGrob(p1, p2, ncol=2, top=textGrob(main_title, gp=gpar(fontsize=14)))
  return(grob)
}


create_scatter_plot <- function(df, scene_id) {
  if(!"AOD_1286" %in% names(df)) stop("AOD_1286 not found in dataframe")
  mod <- rlm(dp ~ AOD_1286, data=df, maxit=200)
  r2 <- round(summary(mod)$r.squared, 3)
  slope <- round(coef(mod)[2], 2)
  stats_label <- paste0("Regression: Slope = ", slope, ", R2 = ", r2)

  ggplot(df, aes(x=AOD_1286, y=dp)) +
    geom_point(alpha=0.3, color="steelblue") +
    geom_density2d(color="grey40", alpha=1) +
    geom_smooth(method = "rlm", method.args = list(maxit = 200),
            color = "red", se = FALSE)  +
    labs(title=paste(scene_id, ": dp vs AOD (1286nm)"),
         subtitle=stats_label,
         x="AOD (1286nm)", y="Delta Pressure (Pa)") +
    theme_bw()
}

create_scatter_plot2 <- function(df, scene_id) {
  mod2 <- lm(xco2 ~ dp, data = df, maxit=200)
  r22 <- round(summary(mod2)$r.squared, 3)
  slope2 <- round(coef(mod2)[2], 2)
  stats_label2 <- paste0("Regression: Slope = ", slope2, ", R2 = ", r22)

  ggplot(df, aes(x=dp, y=xco2)) +
    geom_point(alpha=0.3, color="steelblue") +
    geom_density2d(color="grey40", alpha=1) +
    geom_smooth(method = "lm", method.args = list(maxit = 200),
            color = "red", se = FALSE)  +
    labs(title=paste(scene_id, ": xco2 vs dp"),
         subtitle=stats_label2,
         x="Delta Pressure (Pa)", y="xco2") +
    theme_bw()
}

# ---------------------------------------------------------------------------
# 3. Main Execution
# ---------------------------------------------------------------------------
get_id <- function(path) {
  t_m <- str_extract(path, "t[0-9]+")
  p_m <- str_extract(path, "p[0-9]+")
  paste0(t_m, p_m)
}

files_df <- data.frame(MSAT_file=MSAT_FILES, stringsAsFactors=FALSE) %>%
  mutate(scene_id = sapply(MSAT_file, get_id)) %>%
  distinct(scene_id, .keep_all=TRUE) %>%
  filter(file.exists(MSAT_file))

message(sprintf("Found %d accessible files from list.", nrow(files_df)))

for(i in 1:nrow(files_df)) {
  scene_path <- files_df$MSAT_file[i]
  sid <- files_df$scene_id[i]
  
  cat(sprintf("\n[%d/%d] Processing %s...\n", i, nrow(files_df), sid))
  
  viirs_path <- find_or_download_viirs(scene_path)
  if(is.na(viirs_path)) {
    cat("  -> Skipping: VIIRS data not found.\n")
    next
  }
  
  df_scene <- tryCatch({
    extract_scene_data(viirs_path, scene_path, sid)
  }, error=function(e) {
    cat(sprintf("  -> Extraction Error: %s\n", e$message))
    return(NULL)
  })
  
  if(is.null(df_scene) || nrow(df_scene) < 50) {
    cat("  -> Skipping: Insufficient data.\n")
    next
  }
  
  scene_dir <- file.path(OUTPUT_DIR, sid)
  if(!dir.exists(scene_dir)) dir.create(scene_dir)
  
  cat("  -> Generating Plots...\n")
  
  p_aod_dp   <- create_map_pair(df_scene, "AOD_1286","dp","AOD (1286nm)","Delta Pressure",paste(sid,"- AOD (1286nm) vs DP"))
  ggsave(file.path(scene_dir,"map_AOD_vs_DP.png"), p_aod_dp, width=10, height=5, bg="white")
  p_aod_xco2 <- create_map_pair(df_scene, "AOD_1286","xco2","AOD (1286nm)","XCO2 (ppm)",paste(sid,"- AOD (1286nm) vs XCO2"))
  ggsave(file.path(scene_dir,"map_AOD_vs_XCO2.png"), p_aod_xco2, width=10, height=5, bg="white")
  
  p_dp_alb   <- create_map_pair(df_scene, "dp","alb1286","Delta Pressure","Albedo (1286nm)",paste(sid,"- DP vs Albedo"))
  ggsave(file.path(scene_dir,"map_DP_vs_Albedo.png"), p_dp_alb, width=10, height=5, bg="white")
  
  p_xco2_orog <- create_map_pair(df_scene, "xco2","orog","XCO2 (ppm)","Elevation (m)",paste(sid,"- XCO2 vs Orog"))
  ggsave(file.path(scene_dir,"map_XCO2_vs_Orog.png"), p_xco2_orog, width=10, height=5, bg="white")
  
  p_scatter <- create_scatter_plot(df_scene, sid)
  ggsave(file.path(scene_dir,"scatter_DP_vs_AOD.png"), p_scatter, width=6, height=5, bg="white")
    
  co2_scatter <- create_scatter_plot2(df_scene, sid)
  ggsave(file.path(scene_dir,"scatter_XCO2_vs_DP.png"), co2_scatter, width=6, height=5, bg="white")
  # Save scene collocated dataframe as compressed RDS so you can reload it later
  rds_path <- file.path(scene_dir, paste0(sid, "_scene_data.rds"))
  saveRDS(df_scene, file = rds_path, compress = "xz")
  cat(sprintf("  -> Saved scene data to %s\n", rds_path))

  cat("  -> Saved all plots.\n")
}

message("\nAll processing complete. Check folder: ", OUTPUT_DIR)

