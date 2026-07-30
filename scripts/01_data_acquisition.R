# ==============================================================================
# SCRIPT 01: SPATIAL AND TEMPORAL DATA ACQUISITION
# Project:  Accumulated Thermal Stress in Posidonia oceanica (Balearic Islands)
# Author:   Yeray Gonzalez-Marrero
# Purpose:  Download the Posidonia oceanica distribution layer from IDEBaleares
#           and provide instructions for downloading SST data from Copernicus.
# Outputs:
#   data/raw/posidonia_total_baleares.geojson       <- raw downloaded layer
#   data/raw/posidonia_total_baleares_valida.geojson <- geometry-repaired layer
# ==============================================================================

# Delete workspace
rm(list = ls())

library(sf)
library(dplyr)
library(httr)
library(jsonlite)

# Create the project directory structure if it does not already exist
dirs <- c("data/raw", "data/processed", "outputs/maps", "outputs/rasters",
          "outputs/tables", "outputs/plots", "scripts", "cache")
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ==============================================================================
# SECTION 1: DOWNLOAD OF POSIDONIA OCEANICA DISTRIBUTION POLYGONS (IDEIB)
# ==============================================================================
# Source:   IDEBaleares ArcGIS REST service — GOIB_Posidonia_IB, Layer 2
# URL:      https://ideib.caib.es/geoserveis/rest/services/public/GOIB_Posidonia_IB/MapServer/2
#
# TECHNICAL NOTE: The IDEBaleares ArcGIS REST server enforces a hard limit of
# 1,000 features per request. The layer contains 32,839 polygons, so full
# coverage requires paginated requests using the `resultOffset` parameter.
# The `arcpullr` package does not support automatic pagination for this service;
# therefore, requests are handled directly via `httr` and `jsonlite`.
# ==============================================================================

cat("Downloading Posidonia oceanica distribution layer from IDEBaleares (full pagination)...\n")

url_query      <- "https://ideib.caib.es/geoserveis/rest/services/public/GOIB_Posidonia_IB/MapServer/2/query"
page_size      <- 500   # Conservative page size to avoid server-side timeouts
max_retries    <- 3     # Maximum number of retry attempts per page on network error

# --- Step 1: Retrieve the total feature count from the server ---
resp_count <- GET(url_query,
                  query = list(where = "1=1", returnCountOnly = "true", f = "json"),
                  timeout(30))
total_features <- content(resp_count, as = "parsed")$count
cat(sprintf("  Total polygons in layer: %d\n", total_features))

# --- Step 2: Paginated download loop ---
# Iterates over the full feature set in pages of `page_size`, with automatic
# retry logic to handle transient network failures.
all_features <- list()
offset       <- 0

while (offset < total_features) {
  success <- FALSE

  for (attempt in seq_len(max_retries)) {
    resp <- tryCatch({
      GET(url_query,
          query = list(
            where             = "1=1",
            outFields         = "*",
            f                 = "geojson",
            outSR             = "4326",       # Request coordinates in WGS84
            resultOffset      = offset,
            resultRecordCount = page_size
          ),
          timeout(120))
    }, error = function(e) NULL)

    if (!is.null(resp) && status_code(resp) == 200) {
      data_page <- tryCatch(
        fromJSON(rawToChar(resp$content), simplifyVector = FALSE),
        error = function(e) NULL
      )
      if (!is.null(data_page) && length(data_page$features) > 0) {
        all_features <- c(all_features, data_page$features)
        success      <- TRUE
        break
      }
    }
    cat(sprintf("    Attempt %d failed at offset %d. Retrying in 5 s...\n", attempt, offset))
    Sys.sleep(5)
  }

  if (!success) {
    warning(sprintf("Could not download page at offset %d after %d attempts. Aborting.",
                    offset, max_retries))
    break
  }

  offset <- offset + page_size
  cat(sprintf("  Downloaded %d / %d polygons...\n",
              min(offset, total_features), total_features))
}

# --- Step 3: Reconstruct the full GeoJSON and convert to an sf object ---
geojson_full <- list(type = "FeatureCollection", features = all_features)
geojson_txt  <- toJSON(geojson_full, auto_unbox = TRUE)
posidonia_sf <- st_read(geojson_txt, quiet = TRUE)

# Ensure the CRS is set to WGS84 (EPSG:4326)
if (!is.na(st_crs(posidonia_sf)) && st_crs(posidonia_sf)$epsg != 4326) {
  posidonia_sf <- st_transform(posidonia_sf, crs = 4326)
}

# Save the raw downloaded layer to disk
st_write(posidonia_sf,
         "data/raw/posidonia_total_baleares.geojson",
         delete_dsn = TRUE,
         quiet      = TRUE)

cat(sprintf("\n[OK] Download complete: %d polygons saved to data/raw/posidonia_total_baleares.geojson\n",
            nrow(posidonia_sf)))

# ==============================================================================
# SECTION 1b: GEOMETRY REPAIR AND EXPORT OF VALIDATED LAYER
# ==============================================================================
# The layer downloaded from IDEBaleares frequently contains invalid geometries
# (duplicate vertices, self-intersections) arising from manual digitisation.
# These are repaired here once so that Script 04 can load the validated layer
# directly without recomputing the repair at runtime.
# ==============================================================================

cat("\nRepairing invalid geometries...\n")

n_invalid <- sum(!st_is_valid(posidonia_sf), na.rm = TRUE)
cat(sprintf("  Invalid geometries detected: %d\n", n_invalid))

# st_make_valid() uses the GEOS library to repair duplicate vertices and
# incorrectly wound rings, which are the most common issues in this dataset.
posidonia_valid <- st_make_valid(posidonia_sf)

# Remove any empty geometries that may result from the repair process
posidonia_valid <- posidonia_valid[!st_is_empty(posidonia_valid), ]

n_invalid_post <- sum(!st_is_valid(posidonia_valid), na.rm = TRUE)
cat(sprintf("  Invalid geometries after repair: %d\n", n_invalid_post))
cat(sprintf("  Final feature count: %d\n", nrow(posidonia_valid)))

# Save the repaired layer. This file is the spatial mask used in Script 04.
st_write(posidonia_valid,
         "data/raw/posidonia_total_baleares_valida.geojson",
         delete_dsn = TRUE,
         quiet      = TRUE)

cat("\n[OK] Validated layer saved to: data/raw/posidonia_total_baleares_valida.geojson\n")
cat("     This file is the spatial mask input for Script 04.\n")

# ==============================================================================
# SECTION 2: XARXA DE MONITORATGE SHOOT DENSITY DATA
# ==============================================================================
# The biological calibration data (global shoot density per station and year)
# are compiled manually from the annual monitoring reports of the Xarxa de
# Monitoratge de les Praderies de Posidònia (Govern de les Illes Balears).
#
# Reference:
#   Govern de les Illes Balears (2017, 2019, 2022, 2023, 2024).
#   Memòria de la Xarxa de Monitoratge de les Praderies de Posidònia.
#   Direcció General de Medi Natural i Gestió Forestal.
#   Available at: https://www.caib.es/sites/posidonia/ca/definicia/
#
# The compiled dataset must be placed at:
#   data/raw/xarxa_densidad_haces_completo.xlsx
# ==============================================================================

if (file.exists("xarxa_densidad_haces_completo.xlsx")) {
  file.copy("xarxa_densidad_haces_completo.xlsx",
            "data/raw/xarxa_densidad_haces_completo.xlsx", overwrite = TRUE)
  cat("[OK] Xarxa shoot density data copied to data/raw/\n")
} else {
  warning("File 'xarxa_densidad_haces_completo.xlsx' not found in the working directory. ",
          "Please place the compiled dataset at data/raw/ before running Script 03.")
}

# ==============================================================================
# SECTION 3: SST DATA DOWNLOAD FROM COPERNICUS MARINE SERVICE
# ==============================================================================
#
# VERIFIED DATASET INFORMATION:
# ------------------------------------------------------------------------------
# Dataset ID : cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021
# Variable   : analysed_sst  (in Kelvin; subtract 273.15 to convert to °C)
# Coverage   : 1982-01-01 to present (daily updates)
# Resolution : 0.05° × 0.05° (~5 km)
# Bounding   : lat 30.1–46.0°N, lon -18.1–36.3°E
# Service    : arco-geo-series (supports spatiotemporal subsetting)
#
# IMPORTANT: Use the DATASET ID (prefixed with "cmems_"), not the PRODUCT ID
# shown in the Copernicus web catalogue. The two identifiers differ.
#
# OPTION A — Single download (full period 2015-2024, ~500 MB):
# ------------------------------------------------------------------------------
# Run the following command in the RStudio Terminal:
#
# copernicusmarine subset \
#   --dataset-id cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021 \
#   --variable analysed_sst \
#   --minimum-longitude 1.0 \
#   --maximum-longitude 4.5 \
#   --minimum-latitude 38.0 \
#   --maximum-latitude 40.5 \
#   --start-datetime 2015-05-01T00:00:00 \
#   --end-datetime 2024-09-30T00:00:00 \
#   --output-filename sst_baleares_2015_2024.nc \
#   --output-directory data/raw/
#
# OPTION B — Year-by-year download (recommended; ~50 MB per file):
# ------------------------------------------------------------------------------
# for YEAR in 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024; do
#   copernicusmarine subset \
#     --dataset-id cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021 \
#     --variable analysed_sst \
#     --minimum-longitude 1.0 \
#     --maximum-longitude 4.5 \
#     --minimum-latitude 38.0 \
#     --maximum-latitude 40.5 \
#     --start-datetime ${YEAR}-05-01T00:00:00 \
#     --end-datetime ${YEAR}-09-30T00:00:00 \
#     --output-filename sst_baleares_${YEAR}.nc \
#     --output-directory data/raw/
# done
#
# BOUNDING BOX NOTE:
# The rectangle [lon: 1.0–4.5, lat: 38.0–40.5] covers the entire Balearic
# archipelago, including the Cabrera Archipelago (39.15°N, 2.95°E).
# ==============================================================================

cat("Script 01 complete.\n")
cat("Run the Copernicus command above in the Terminal to download the SST data.\n")
