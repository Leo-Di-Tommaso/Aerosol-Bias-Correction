library(ncdf4)
library(ggplot2)

MSAT_ROOT_DIR <- "/mnt/gcs/msat-prod-data-methanesat-level2/level2_results"

# ---- 1. Point at one scene ----
# A "scene" is typically one collection/target granule. Adjust the pattern to your naming.
scene_file <- list.files(
  MSAT_ROOT_DIR,
  pattern = "\\.nc$",
  recursive = TRUE,
  full.names = TRUE
)[1]                      # replace [1] with the specific granule you want

stopifnot(length(scene_file) == 1, file.exists(scene_file))

nc <- nc_open(scene_file)

# ---- 2. Inspect structure first (do not assume variable paths) ----
# Print every variable with its full group path so you can confirm names.
all_vars <- names(nc$var)
print(all_vars)
# Also check groups, since MethaneSAT L2 is group-structured:
print(names(nc$dim))

# ---- 3. Helper: case-insensitive find of a variable by keyword ----
find_var <- function(nc, keyword) {
  hits <- grep(keyword, names(nc$var), ignore.case = TRUE, value = TRUE)
  if (length(hits) == 0) return(NA_character_)
  hits
}

ch4_candidates    <- find_var(nc, "ch4")        # e.g. methane column / xch4
albedo_candidates <- find_var(nc, "albedo")     # surface albedo / reflectance
co2_candidates    <- find_var(nc, "co2")        # proxy CO2 column, if present

message("CH4 candidates:\n",    paste(ch4_candidates,    collapse = "\n"))
message("Albedo candidates:\n", paste(albedo_candidates, collapse = "\n"))
message("CO2 candidates:\n",    paste(co2_candidates,    collapse = "\n"))

# >>> STOP here on first run, read the printout, then set exact names below. <
ch4_var    <- ch4_candidates[1]      # replace with the confirmed name
albedo_var <- albedo_candidates[1]   # replace with the confirmed name
co2_var    <- if (length(co2_candidates) && !is.na(co2_candidates[1])) co2_candidates[1] else NA

# ---- 4. Read per-pixel arrays and flatten ----
read_flat <- function(nc, varname) {
  if (is.na(varname)) return(NULL)
  v <- ncvar_get(nc, varname)
  fill <- ncatt_get(nc, varname, "_FillValue")
  if (fill$hasatt) v[v == fill$value] <- NA
  as.vector(v)
}

ch4    <- read_flat(nc, ch4_var)
albedo <- read_flat(nc, albedo_var)
co2    <- read_flat(nc, co2_var)

nc_close(nc)

# ---- 5. Build a tidy per-pixel data frame ----
df <- data.frame(albedo = albedo, ch4 = ch4)
if (!is.null(co2)) df$co2 <- co2
df <- df[is.finite(df$albedo) & is.finite(df$ch4), ]

# ---- 6. Plot product vs albedo, same pixel ----
p_ch4 <- ggplot(df, aes(albedo, ch4)) +
  geom_point(alpha = 0.3, size = 0.6) +
  labs(x = "Retrieved albedo", y = "Retrieved CH4", title = "CH4 vs albedo (per pixel)") +
  theme_minimal()
print(p_ch4)

if (!is.null(df$co2)) {
  p_co2 <- ggplot(df, aes(albedo, co2)) +
    geom_point(alpha = 0.3, size = 0.6) +
    labs(x = "Retrieved albedo", y = "Proxy CO2 column", title = "CO2 vs albedo (per pixel)") +
    theme_minimal()
  print(p_co2)
}