# ==============================================================================
# SCRIPT 03: EMPIRICAL GAMM CALIBRATION AND ADVANCED MODEL DIAGNOSTICS
# Project:  Accumulated Thermal Stress in Posidonia oceanica (Balearic Islands)
# Author:   Yeray González Marrero
# Purpose:  Calibrate the DHD-to-shoot-density-change relationship using field
#           data from the Xarxa de Monitoratge. Runs a model tournament
#           (Gaussian vs. Scaled-t, with and without depth), selects the best
#           model, extracts dynamic management thresholds from the model curve,
#           and exports all validation diagnostics and figures.
# Inputs:   outputs/rasters/dhd_posidonia_baleares.tif
#           data/raw/xarxa_densidad_haces_completo.xlsx
# Outputs:  models/modelo_scat_profundidad.rds  (winning model)
#           models/validation/                  (diagnostics and summaries)
#           outputs/tables/umbrales_gestion.xlsx
#           outputs/plots/thermal_stress_response.png/.svg
#           outputs/plots/temporal_impact_evolution.png/.svg
# ==============================================================================

# Delete workspace
rm(list = ls())

library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(mgcv)
library(writexl)
library(readxl)
library(stringr)
library(ggplot2)
library(svglite)

cat("Starting GAMM calibration and advanced model validation...\n")

# ==============================================================================
# SECTION 0: PUBLICATION-QUALITY GGPLOT2 THEME
# ==============================================================================
theme_publication <- function(base_size = 11, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title    = element_text(size = rel(1.3), face = "bold",  hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_text(size = rel(1.0), face = "plain", hjust = 0, margin = margin(b = 12)),
      plot.caption  = element_text(size = rel(0.8), face = "italic", hjust = 1, margin = margin(t = 10)),
      axis.title    = element_text(size = rel(1.0), face = "bold"),
      axis.title.x  = element_text(margin = margin(t = 10)),
      axis.title.y  = element_text(angle = 90, margin = margin(r = 10)),
      axis.text     = element_text(size = rel(0.9), color = "black"),
      legend.title  = element_text(size = rel(0.9), face = "bold"),
      legend.text   = element_text(size = rel(0.9)),
      legend.position = "bottom"
    )
}

# ==============================================================================
# SECTION 1: CREATE OUTPUT DIRECTORIES
# ==============================================================================
dirs <- c("models", "models/validation", "outputs/tables", "outputs/plots")
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# SECTION 2: DATA PREPARATION AND SPATIO-TEMPORAL MATCHING
# ==============================================================================
# Load the shoot density monitoring data (Xarxa de Monitoratge)
stations <- read_excel("data/raw/xarxa_densidad_haces_completo.xlsx")

# Convert station coordinates to an sf object for spatial extraction
unique_stations <- stations %>% select(id_estacion, lat, lon) %>% distinct()
stations_sf     <- st_as_sf(unique_stations, coords = c("lon", "lat"), crs = 4326)

# Load the DHD raster stack (one band per year, output of Script 02)
dhd_raster <- rast("outputs/rasters/dhd_posidonia_baleares.tif")
crs(dhd_raster) <- "EPSG:4326"

# Sanity check: ensure the DHD raster is not empty
max_vals <- global(dhd_raster, "max", na.rm = TRUE)
if (all(is.na(max_vals[, 1]) | is.nan(max_vals[, 1]) | is.infinite(max_vals[, 1]))) {
  stop("DHD raster 'dhd_posidonia_baleares.tif' is empty. Re-run Script 02.")
}

# Extract DHD values at each monitoring station location
cat("Extracting DHD values at monitoring station coordinates...\n")
dhd_extracted          <- terra::extract(dhd_raster, vect(stations_sf))
dhd_extracted$id_estacion <- unique_stations$id_estacion

# Reshape from wide to long format and parse the year from band names (DHD_YYYY)
dhd_long <- dhd_extracted %>%
  pivot_longer(cols = -c(ID, id_estacion),
               names_to  = "year_dhd",
               values_to = "dhd_accumulated") %>%
  mutate(
    year_dhd    = as.integer(sub("DHD_", "", year_dhd)),
    id_estacion = as.character(id_estacion)
  ) %>%
  filter(!is.na(dhd_accumulated) & !is.nan(dhd_accumulated)) %>%
  select(id_estacion, year_dhd, dhd_accumulated)

# Compute the annualised interannual rate of change in global shoot density
# (% yr-1), accounting for unequal sampling intervals between consecutive years.
delta_data <- stations %>%
  mutate(id_estacion = as.character(id_estacion), anyo = as.integer(anyo)) %>%
  filter(!is.na(densidad_global_feixos_m2)) %>%
  arrange(id_estacion, anyo) %>%
  group_by(id_estacion) %>%
  mutate(
    prev_year    = lag(anyo),
    prev_density = lag(densidad_global_feixos_m2),
    delta_pct_yr = ((densidad_global_feixos_m2 - prev_density) / prev_density * 100) /
                   (anyo - prev_year)
  ) %>%
  filter(!is.na(delta_pct_yr)) %>%
  ungroup() %>%
  mutate(id_pradera = as.factor(str_extract(id_estacion, "^[0-9]+")))

# Merge the biological response variable with the thermal predictor via
# a direct year-synchronised join (DHD of year Y matched to density change
# measured at the end of year Y).
model_data <- delta_data %>%
  inner_join(dhd_long, by = c("id_estacion" = "id_estacion", "anyo" = "year_dhd")) %>%
  filter(!is.na(dhd_accumulated) & !is.na(profundidad_m))

write.csv(model_data, "outputs/tables/model_calibration_data.csv", row.names = FALSE)
cat(sprintf("  [OK] Matched records (direct year synchronisation): %d rows.\n", nrow(model_data)))

# ==============================================================================
# SECTION 3: MODEL TOURNAMENT — GAUSSIAN VS. SCALED-T, WITH AND WITHOUT DEPTH
# ==============================================================================
# We fit four candidate GAMMs to select the most appropriate error distribution
# and covariate structure. The scaled t-distribution (scat) is included because
# the response variable (shoot density change) exhibits heavy-tailed behaviour
# with occasional extreme mortality and compensatory growth events.
# All models use REML estimation and include meadow identity as a random effect
# to account for the nested structure of the data (stations within meadows).
cat("\n[Phase 1-2] Fitting GAMM model tournament (Gaussian vs. Scaled-t)...\n")

mod_gaussian_depth <- gam(
  delta_pct_yr ~ s(dhd_accumulated, k = 3) + s(profundidad_m, k = 3) + s(id_pradera, bs = "re"),
  data = model_data, family = gaussian(link = "identity"), method = "REML", select = TRUE
)

mod_gaussian_simple <- gam(
  delta_pct_yr ~ s(dhd_accumulated, k = 3) + s(id_pradera, bs = "re"),
  data = model_data, family = gaussian(link = "identity"), method = "REML", select = TRUE
)

# Winning model: Scaled-t with depth covariate
mod_scat_depth <- gam(
  delta_pct_yr ~ s(dhd_accumulated, k = 3) + s(profundidad_m, k = 3) + s(id_pradera, bs = "re"),
  data = model_data, family = scat(link = "identity"), method = "REML", select = TRUE
)

mod_scat_simple <- gam(
  delta_pct_yr ~ s(dhd_accumulated, k = 3) + s(id_pradera, bs = "re"),
  data = model_data, family = scat(link = "identity"), method = "REML", select = TRUE
)

# Serialise all four models to disk for use in Script 04
saveRDS(mod_gaussian_depth,  "models/modelo_gaussiano_profundidad.rds")
saveRDS(mod_gaussian_simple, "models/modelo_gaussiano_simple.rds")
saveRDS(mod_scat_depth,      "models/modelo_scat_profundidad.rds")
saveRDS(mod_scat_simple,     "models/modelo_scat_simple.rds")

# ==============================================================================
# SECTION 4: MODEL SUMMARIES, PERFORMANCE METRICS, AND RESIDUAL DIAGNOSTICS
# ==============================================================================
cat("\n[Phase 3] Exporting model summaries, metrics, and residual diagnostics...\n")

# Export full model summaries as plain text
capture.output(summary(mod_gaussian_depth),  file = "models/validation/summary_gaussian_depth.txt")
capture.output(summary(mod_gaussian_simple), file = "models/validation/summary_gaussian_simple.txt")
capture.output(summary(mod_scat_depth),      file = "models/validation/summary_scat_depth.txt")
capture.output(summary(mod_scat_simple),     file = "models/validation/summary_scat_simple.txt")

sum_g_depth  <- summary(mod_gaussian_depth)
sum_g_simple <- summary(mod_gaussian_simple)
sum_s_depth  <- summary(mod_scat_depth)
sum_s_simple <- summary(mod_scat_simple)

# Compile smooth-term p-value tables across all models
df_g_depth  <- as.data.frame(sum_g_depth$s.table)  %>% mutate(Variable = rownames(.), Model = "Gaussian With Depth")
df_g_simple <- as.data.frame(sum_g_simple$s.table) %>% mutate(Variable = rownames(.), Model = "Gaussian Without Depth")
df_s_depth  <- as.data.frame(sum_s_depth$s.table)  %>% mutate(Variable = rownames(.), Model = "Scaled-t With Depth")
df_s_simple <- as.data.frame(sum_s_simple$s.table) %>% mutate(Variable = rownames(.), Model = "Scaled-t Without Depth")

standardise_names <- function(df) {
  colnames(df)[1:4] <- c("edf", "Ref.df", "F_Statistic", "p_value")
  return(df)
}

smooth_terms_table <- bind_rows(
  standardise_names(df_g_depth), standardise_names(df_g_simple),
  standardise_names(df_s_depth), standardise_names(df_s_simple)
) %>% select(Model, Variable, edf, Ref.df, F_Statistic, p_value, everything())

write_xlsx(smooth_terms_table, "models/validation/smooth_terms_p_values.xlsx")

# Helper functions to safely extract convergence and GCV scores
extract_convergence <- function(model) {
  if (!is.null(model$converged)) return(model$converged)
  return(TRUE)
}
extract_gcv <- function(model) {
  if (!is.null(model$gcv.ubre) && length(model$gcv.ubre) == 1) return(as.numeric(model$gcv.ubre))
  return(NA)
}
extract_r_squared <- function(model_summary) {
  if (!is.null(model_summary$r.sq)) return(as.numeric(model_summary$r.sq))
  return(NA)
}

# Global performance metrics table
global_metrics <- data.frame(
  Model = c("Gaussian With Depth", "Gaussian Without Depth",
            "Scaled-t With Depth", "Scaled-t Without Depth"),
  Adjusted_R_Squared    = c(extract_r_squared(sum_g_depth), extract_r_squared(sum_g_simple),
                             extract_r_squared(sum_s_depth), extract_r_squared(sum_s_simple)),
  Deviance_Explained_Pct = c(sum_g_depth$dev.expl * 100, sum_g_simple$dev.expl * 100,
                              sum_s_depth$dev.expl * 100, sum_s_simple$dev.expl * 100),
  AIC                   = c(AIC(mod_gaussian_depth), AIC(mod_gaussian_simple),
                             AIC(mod_scat_depth),    AIC(mod_scat_simple)),
  GCV_Score             = c(extract_gcv(mod_gaussian_depth), extract_gcv(mod_gaussian_simple),
                             extract_gcv(mod_scat_depth),    extract_gcv(mod_scat_simple)),
  Optimiser_Converged   = c(extract_convergence(mod_gaussian_depth), extract_convergence(mod_gaussian_simple),
                             extract_convergence(mod_scat_depth),    extract_convergence(mod_scat_simple))
)
write_xlsx(global_metrics, "models/validation/global_validation_metrics.xlsx")

# Residual diagnostic plots (Q-Q, residuals vs. fitted, histogram, response)
export_diagnostics <- function(model, label) {
  png(paste0("models/validation/residual_diagnostics_", label, ".png"),
      width = 1200, height = 900, res = 150)
  par(mfrow = c(2, 2))
  capture.output(gam.check(model),
                 file = paste0("models/validation/gam_check_console_", label, ".txt"))
  dev.off()
}

export_diagnostics(mod_gaussian_depth,  "gaussian_depth")
export_diagnostics(mod_gaussian_simple, "gaussian_simple")
export_diagnostics(mod_scat_depth,      "scat_depth")
export_diagnostics(mod_scat_simple,     "scat_simple")

# ==============================================================================
# SECTION 5: DYNAMIC EXTRACTION OF OPERATIONAL MANAGEMENT THRESHOLDS
# ==============================================================================
# Two thresholds are derived from the winning model (Scaled-t with depth):
#
#   Critical threshold: the DHD value at which the model prediction crosses
#   zero (i.e., the transition from net growth to net regression). Identified
#   as the minimum absolute value of the prediction curve over a fine grid.
#
#   Alert threshold: the 75th percentile of DHD values observed at stations
#   that had not yet reached the critical threshold. This represents the upper
#   boundary of the sub-lethal stress range documented in the field.
cat("\n[Phase 4] Extracting dynamic management thresholds from the model curve...\n")

median_depth <- median(model_data$profundidad_m, na.rm = TRUE)

# Build a fine prediction grid over the full DHD range, holding depth constant
# at the dataset median and excluding the meadow random effect.
dhd_grid <- data.frame(
  dhd_accumulated = seq(0, max(model_data$dhd_accumulated, na.rm = TRUE) + 5,
                        length.out = 10000),
  profundidad_m   = median_depth,
  id_pradera      = model_data$id_pradera[1]
)
dhd_grid$pred <- predict(mod_scat_depth, newdata = dhd_grid,
                         type = "response", exclude = "s(id_pradera)")

# Critical threshold: zero-crossing of the prediction curve
zero_crossing_idx <- which.min(abs(dhd_grid$pred))
critical_threshold <- dhd_grid$dhd_accumulated[zero_crossing_idx]
cat(sprintf("  -> Critical threshold (zero-crossing): %.2f DHD\n", critical_threshold))

# Alert threshold: 75th percentile of sub-lethal observations
alert_threshold <- quantile(
  model_data$dhd_accumulated[model_data$dhd_accumulated < critical_threshold],
  0.75, na.rm = TRUE
)

thresholds_table <- data.frame(
  Category  = c("Alert Threshold (Green-Yellow)", "Critical Threshold (Yellow-Red)"),
  DHD_Value = c(alert_threshold, critical_threshold),
  Justification = c(
    "75th percentile of DHD values at stations below the critical threshold (upper bound of sub-lethal stress range)",
    "Zero-crossing of the Scaled-t GAMM prediction curve (onset of net shoot regression)"
  )
)
write_xlsx(thresholds_table, "outputs/tables/management_thresholds.xlsx")

# ==============================================================================
# SECTION 6: VALIDATION FIGURES
# ==============================================================================
cat("\n[Phase 5] Generating validation and temporal evolution figures...\n")

# --- 6.1 Spatial Cross-Validation (Leave-One-Island-Out) ---
# Assesses the transferability of the thermal signal to unseen geographic regions.
# Each island is held out in turn as the test set; the model is refitted on the
# remaining islands and predictions are made for the held-out island.
islands    <- unique(as.character(model_data$isla))
cv_results <- data.frame()

for (test_island in islands) {
  train_set <- model_data %>% filter(isla != test_island)
  test_set  <- model_data %>% filter(isla == test_island)

  if (nrow(train_set) < 10 || nrow(test_set) == 0) next

  cv_model <- tryCatch(
    gam(delta_pct_yr ~ s(dhd_accumulated, k = 3) + s(profundidad_m, k = 3),
        data = train_set, family = scat(link = "identity"),
        method = "REML", select = TRUE),
    error = function(e) NULL
  )
  if (is.null(cv_model)) next

  cv_predictions <- predict(cv_model, newdata = test_set, type = "response")
  cv_results <- rbind(cv_results,
                      data.frame(island_held_out = test_island,
                                 observed        = test_set$delta_pct_yr,
                                 predicted       = as.numeric(cv_predictions)))
}

if (nrow(cv_results) > 0) {
  rmse_cv <- sqrt(mean((cv_results$observed - cv_results$predicted)^2))
  r2_cv   <- cor(cv_results$observed, cv_results$predicted)^2

  p_cv <- ggplot(cv_results, aes(x = observed, y = predicted, colour = island_held_out)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.8) +
    scale_colour_brewer(palette = "Set1", name = "Island held out") +
    annotate("text", x = min(cv_results$observed), y = max(cv_results$predicted),
             label = sprintf("R²_CV = %.2f\nRMSE = %.1f %%/yr", r2_cv, rmse_cv),
             hjust = 0, vjust = 1, size = 4) +
    labs(title    = "Spatial Cross-Validation (Leave-One-Island-Out)",
         subtitle = "Predictive capacity of the thermal component for unseen geographic regions",
         x        = "Observed density change (% yr\u207b\u00b9)",
         y        = "Predicted density change (% yr\u207b\u00b9)") +
    theme_publication()

  ggsave("models/validation/spatial_cv_observed_vs_predicted.png",
         plot = p_cv, width = 6, height = 6, dpi = 300, bg = "white")
}

# --- 6.2 External Validation: Comparison with Marbà & Duarte (2010) ---
# The empirical linear model from Marbà & Duarte (2010) is applied to the
# current dataset to assess whether the historical relationship still holds
# under contemporary Mediterranean thermal conditions.
marba_slope  <- -0.1
model_data   <- model_data %>% mutate(pred_marba = marba_slope * dhd_accumulated)

rmse_marba <- sqrt(mean((model_data$delta_pct_yr - model_data$pred_marba)^2))
r2_marba   <- cor(model_data$delta_pct_yr, model_data$pred_marba)^2

p_marba <- ggplot(model_data, aes(x = pred_marba, y = delta_pct_yr)) +
  geom_point(aes(colour = isla), size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_smooth(method = "lm", se = TRUE, colour = "steelblue", linewidth = 0.8) +
  scale_colour_brewer(palette = "Set1", name = "Island") +
  annotate("text", x = min(model_data$pred_marba), y = max(model_data$delta_pct_yr),
           label = sprintf("R² = %.2f\nRMSE = %.1f %%/yr", r2_marba, rmse_marba),
           hjust = 0, vjust = 1, size = 4) +
  labs(title    = "External Validation: Marbà & Duarte (2010) applied to the Balearic Islands",
       subtitle = "Historical empirical predictions versus current field observations",
       x        = "Predicted density change — Marbà & Duarte (% yr\u207b\u00b9)",
       y        = "Observed density change in the Balearic Islands (% yr\u207b\u00b9)") +
  theme_publication()

ggsave("models/validation/external_validation_marba_duarte.png",
       plot = p_marba, width = 6, height = 6, dpi = 300, bg = "white")

# --- 6.3 GAMM Calibration Curve vs. Marbà & Duarte (2010) ---
dev_explained_scat <- sum_s_depth$dev.expl * 100

p_calib <- ggplot(model_data, aes(x = dhd_accumulated, y = delta_pct_yr)) +
  geom_point(aes(colour = isla, shape = isla), size = 3, alpha = 0.8) +
  geom_line(data = dhd_grid,
            aes(y = pred, linetype = "GAMM (this study)"),
            colour = "darkred", linewidth = 1) +
  geom_line(data = dhd_grid,
            aes(y = marba_slope * dhd_accumulated, linetype = "Marbà & Duarte (2010)"),
            colour = "steelblue", linewidth = 1) +
  geom_vline(xintercept = critical_threshold, linetype = "dashed", colour = "black", alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dotted") +
  annotate("text", x = critical_threshold + 0.5, y = max(model_data$delta_pct_yr),
           label = sprintf("Critical threshold\n(%.1f DHD)", critical_threshold),
           hjust = 0, size = 3.5) +
  scale_linetype_manual(name   = "Model",
                        values = c("GAMM (this study)" = "solid",
                                   "Marbà & Duarte (2010)" = "dashed")) +
  scale_colour_brewer(palette = "Set1", name = "Island") +
  labs(
    title    = "Calibration: Posidonia oceanica Response to Thermal Stress",
    subtitle = sprintf("Scaled-t GAMM with depth (deviance explained = %.1f%%) vs. Marbà & Duarte (2010)",
                       dev_explained_scat),
    x        = "Accumulated DHD (summer season, > 28 °C)",
    y        = "Global shoot density change (% yr\u207b\u00b9)"
  ) +
  theme_publication()

ggsave("outputs/plots/thermal_stress_response.png",
       plot = p_calib, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("outputs/plots/thermal_stress_response.svg",
       plot = p_calib, width = 8, height = 6, bg = "white")

# --- 6.4 Regional Temporal Evolution of Thermal Impact ---
# Computes the mean predicted density change across all stations for each year,
# providing a regional-scale index of thermal impact over the study period.
model_data$regional_pred <- predict(mod_scat_depth, newdata = model_data,
                                    type = "response", exclude = "s(id_pradera)")

temporal_data <- model_data %>%
  group_by(anyo) %>%
  summarise(mean_pred = mean(regional_pred, na.rm = TRUE)) %>%
  mutate(status = ifelse(mean_pred >= 0, "positive", "negative"))

p_temporal <- ggplot(temporal_data, aes(x = anyo, y = mean_pred)) +
  # Shaded region for years without field sampling (model-estimated period)
  geom_rect(aes(xmin = 2020, xmax = 2021, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.5, inherit.aes = FALSE) +
  annotate("text", x = 2020.5, y = 4.5,
           label = "Estimated period\n(no field sampling)",
           colour = "grey40", size = 2.8, fontface = "italic",
           angle = 0, vjust = 1, hjust = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_smooth(method = "loess", colour = "black", fill = "grey80", alpha = 0.5) +
  geom_point(aes(colour = status), size = 4) +
  scale_colour_manual(values = c("positive" = "#387635", "negative" = "#B7372F")) +
  labs(
    title    = "Regional Temporal Evolution of Thermal Impact",
    subtitle = "Mean expected shoot density change across the Balearic ecoregion",
    x        = "Year",
    y        = "Mean regional expected density change (% yr\u207b\u00b9)"
  ) +
  theme_publication() +
  theme(legend.position = "none")

ggsave("outputs/plots/temporal_impact_evolution.png",
       plot = p_temporal, width = 8, height = 5, dpi = 300, bg = "white")
ggsave("outputs/plots/temporal_impact_evolution.svg",
       plot = p_temporal, width = 8, height = 5, bg = "white")

cat("\n[COMPLETE] Script 03 finished: GAMM calibration and validation successful.\n")
