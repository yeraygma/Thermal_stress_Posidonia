# ==============================================================================
# SCRIPT 02: THERMAL PROCESSING AND DHD INDEX CALCULATION
# Project:  Accumulated Thermal Stress in Posidonia oceanica (Balearic Islands)
# Author:   Yeray Gonzalez-Marrero
# Purpose:  Process daily SST NetCDF files from Copernicus to compute the
#           Degree-Heating-Days (DHD) index for each year in the study period.
#           DHD quantifies cumulative thermal stress above the acute stress
#           threshold for Posidonia oceanica (28.0 °C; Marbà & Duarte, 2010).
# Input:    data/raw/sst_baleares_2015_2024.nc
# Output:   outputs/rasters/dhd_posidonia_baleares.tif  (multiband GeoTIFF,
#           one band per year, named DHD_YYYY)
# ==============================================================================


# Delete workspace
rm(list = ls())

library(terra)
library(lubridate)

cat("--- STARTING COPERNICUS NETCDF THERMAL PROCESSING ---\n")

# --- Locate the input NetCDF file ---
nc_file <- "data/raw/sst_baleares_2015_2024.nc"
if (!file.exists(nc_file)) {
  stop("Input file not found: data/raw/sst_baleares_2015_2024.nc\n",
       "Run Script 01 and download SST data from Copernicus before proceeding.")
}

# ==============================================================================
# STEP 1: Load the SST variable from the NetCDF file
# ==============================================================================
# We explicitly request the 'analysed_sst' subdataset to avoid loading
# ancillary variables (e.g., quality flags, error estimates) by accident.
cat("Loading variable 'analysed_sst'...\n")
sst_raw <- rast(nc_file, subds = "analysed_sst")

# ==============================================================================
# STEP 2: Automatic unit detection and Kelvin-to-Celsius conversion
# ==============================================================================
# Copernicus L4 SST products are distributed in Kelvin. We detect the unit
# automatically by inspecting the maximum value of the first layer: values
# above 100 are unambiguously Kelvin (ocean SST never exceeds ~35 °C ≈ 308 K).
val_max <- global(sst_raw[[1]], "max", na.rm = TRUE)[1, 1]
if (is.na(val_max)) {
  stop("The first layer of the NetCDF is entirely NA. ",
       "Check the Copernicus download for corruption or an empty spatial subset.")
}

if (val_max > 100) {
  cat(sprintf("Kelvin units detected (max value: %.1f K). Converting to Celsius...\n", val_max))
  sst_celsius <- sst_raw - 273.15
} else {
  cat(sprintf("Celsius units detected (max value: %.1f °C). No conversion needed.\n", val_max))
  sst_celsius <- sst_raw
}

# ==============================================================================
# STEP 3: Extract temporal metadata and filter to the summer growing season
# ==============================================================================
# The DHD index is computed only for the summer growing season (May–October),
# when Posidonia oceanica is physiologically active and susceptible to
# heat-induced mortality (Marbà & Duarte, 2010).
time_vector <- time(sst_celsius)
if (all(is.na(time_vector))) {
  stop("terra could not parse the time dimension from the NetCDF calendar. ",
       "Verify that the file contains a valid CF-compliant time axis.")
}

months_all <- month(time_vector)
years_all  <- year(time_vector)

# Retain only layers falling within the summer window (May = 5, October = 10)
summer_idx      <- which(months_all >= 5 & months_all <= 10)
sst_summer      <- sst_celsius[[summer_idx]]
years_summer    <- years_all[summer_idx]

cat(sprintf("Processing %d summer-season daily layers...\n", length(summer_idx)))

# ==============================================================================
# STEP 4: Compute daily thermal excess above the acute stress threshold
# ==============================================================================
# For each day, the thermal excess is defined as max(SST - 28.0, 0).
# We use terra::clamp() rather than ifelse() to avoid introducing spurious NAs
# in cells where SST is below the threshold.
cat("Computing daily thermal exceedances above 28.0 °C...\n")
thermal_excess <- sst_summer - 28.0
thermal_excess <- clamp(thermal_excess, lower = 0, values = TRUE)

# ==============================================================================
# STEP 5: Accumulate daily exceedances to annual DHD totals (map algebra)
# ==============================================================================
# For each year, sum all daily thermal exceedances to produce the annual DHD.
# A land mask derived from the first SST layer is applied after summation to
# restore NA values over land: terra::sum(..., na.rm = TRUE) converts land
# cells (which are NA in the SST data) to zero, which is ecologically incorrect.
unique_years   <- unique(years_summer)
annual_dhd_list <- list()

land_mask <- sst_summer[[1]] * 0 + 1   # 1 over sea, NA over land

for (yr in unique_years) {
  cat(sprintf("  -> Accumulating DHD for year %d...\n", yr))
  year_idx  <- which(years_summer == yr)
  dhd_year  <- sum(thermal_excess[[year_idx]], na.rm = TRUE)

  # Re-apply the land mask to restore NA values over land cells
  dhd_year  <- dhd_year * land_mask

  annual_dhd_list[[as.character(yr)]] <- dhd_year
}

# ==============================================================================
# STEP 6: Assemble the multiband output raster and assign metadata
# ==============================================================================
dhd_stack <- rast(annual_dhd_list)
names(dhd_stack) <- paste0("DHD_", unique_years)
crs(dhd_stack)   <- "EPSG:4326"   # Explicitly set CRS to WGS84

# ==============================================================================
# STEP 7: Write output and report summary statistics
# ==============================================================================
output_tif <- "outputs/rasters/dhd_posidonia_baleares.tif"
writeRaster(dhd_stack, output_tif, overwrite = TRUE)

max_dhd_by_year <- global(dhd_stack, "max", na.rm = TRUE)
cat("\n[OK] DHD raster successfully generated and saved to:\n")
cat(sprintf("     %s\n", output_tif))
cat(sprintf("     Peak historical DHD: %.2f degree-days\n",
            max(max_dhd_by_year[, 1], na.rm = TRUE)))
