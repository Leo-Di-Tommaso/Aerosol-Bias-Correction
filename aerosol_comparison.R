library(dplyr)
library(readr)

# Read the CSV
models <- read_csv("by_aerosol_data/model_comparison_all.csv")

# Get unique aerosol types
aerosol_types <- unique(models$aerosol_type)

for (aerosol in sort(aerosol_types)) {
  cat("\n", paste(rep("=", 80), collapse=""), "\n", sep="")
  cat("AEROSOL TYPE: ", aerosol, "\n")
  cat(paste(rep("=", 80), collapse=""), "\n", sep="")
  
  # Filter for this aerosol type
  aerosol_models <- models %>%
    filter(aerosol_type == aerosol) %>%
    arrange(aic)
  
  # Get unique AE categories for this aerosol
  ae_cats <- unique(aerosol_models$ae_category)
  
  # Overall stats
  cat(sprintf("\nTotal models: %d\n", nrow(aerosol_models)))
  cat(sprintf("AIC range: %.2f (min: %.2f, max: %.2f)\n", 
              max(aerosol_models$aic) - min(aerosol_models$aic),
              min(aerosol_models$aic), max(aerosol_models$aic)))
  cat(sprintf("AE Categories: %s\n\n", paste(ae_cats, collapse=", ")))
  
  # Show top 20 models with AE separation
  top_20 <- aerosol_models %>% slice(1:20)
  
  cat("TOP 20 MODELS BY AIC:\n")
  cat(sprintf("%-3s | %-12s | %-40s | %8s | %7s | %6s | %7s\n", 
              "Rank", "AE Category", "Model Spec", "AIC", "R2", "adjR2", "max_VIF"))
  cat(paste(rep("-", 115), collapse=""), "\n")
  
  for (i in 1:nrow(top_20)) {
    row <- top_20[i, ]
    cat(sprintf("%3d | %-12s | %-40s | %8.2f | %7.4f | %6.4f | %7.2f\n",
                i, row$ae_category, row$model_spec, row$aic, 
                row$r_squared, row$adj_r_squared, row$max_vif))
  }
  
  # Summary by AE category
  cat("\n\nSUMMARY BY AE CATEGORY:\n")
  cat(sprintf("%-15s | Count | Min AIC | Max AIC | Range | Mean R2 | Mean adjR2\n", "AE Category"))
  cat(paste(rep("-", 95), collapse=""), "\n")
  
  for (ae_cat in sort(ae_cats)) {
    subset_ae <- aerosol_models %>% filter(ae_category == ae_cat)
    cat(sprintf("%-15s | %5d | %7.2f | %7.2f | %5.2f | %7.4f | %7.4f\n",
                ae_cat, nrow(subset_ae),
                min(subset_ae$aic), max(subset_ae$aic),
                max(subset_ae$aic) - min(subset_ae$aic),
                mean(subset_ae$r_squared), mean(subset_ae$adj_r_squared)))
  }
}
