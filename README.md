# Accumulated Thermal Stress Traffic Light (DHD) for *Posidonia oceanica* Meadows in the Balearic Islands

[![Read Final Report DHD Posidonia](https://img.shields.io/badge/Leer_Informe-PDF-red?style=for-the-badge&logo=adobeacrobatreader&logoColor=white)](documents/report_DHD_Posidonia.pdf)


## Project Description

This repository contains the complete R analysis pipeline for the thermal stress indicator on *Posidonia oceanica* in the Balearic Islands.

The indicator integrates satellite remote sensing data, official benthic cartography, and in situ empirical calibration to quantify the accumulated physiological stress exerted by thermal anomalies on *Posidonia oceanica* meadows. Using a Generalised Additive Mixed Model (GAMM) with a scaled t-distribution, two operational management thresholds are identified: the alert threshold (4.9 DHD) and the critical threshold for net regression (5.6 DHD). These thresholds enable the spatial classification of the meadows into three states of thermal vulnerability.



## Project Structure

```
posidonia-dhd-baleares/
├── data/
│   └── raw/
│       ├── posidonia_total_baleares_valida.geojson  # Official Posidonia distribution (IDEIB, Layer 2)
│       ├── xarxa_densidad_haces_completo.xlsx       # Global shoot density data (Xarxa de Monitoratge)
│       └── sst_baleares_YYYY.nc                     # Annual SST data (Copernicus CMEMS)
├── documents/
    └── report_DHD_Posidonia.pdf                     # Final report
├── models/
│   ├── modelo_scat_profundidad.rds                  # Winning GAMM model (Scaled-t with depth)
│   └── validation/                                  # Diagnostics and validation metrics
├── outputs/
│   ├── maps/                                        # Thermal traffic light maps (PNG/SVG)
│   ├── rasters/                                     # DHD rasters and classified traffic lights (GeoTIFF)
│   └── tables/                                      # Management thresholds and calibration data
├── scripts/
│   ├── 01_data_acquisition.R                        # Download and repair of Posidonia layer (IDEIB)
│   ├── 02_calculate_dhd.R                           # DHD index calculation from SST data (NetCDF)
│   ├── 03_gamm_calibration.R                        # GAMM calibration, model selection, and thresholds
│   └── 04_visualization_and_mapping.R               # Generation of the 5 cartographic concepts
├── .gitignore
├── README.md                                        # This file
└── posidonia-dhd-baleares.Rproj                     # RStudio Project
```

## Data Sources

| Data | Source | Access |
|---|---|---|
| Sea Surface Temperature (SST) | Copernicus Marine Service (CMEMS) — `cmems_SST_MED_SST_L4_REP_OBSERVATIONS_010_021` | Public, requires free registration at [marine.copernicus.eu](https://marine.copernicus.eu) |
| *Posidonia oceanica* distribution | IDEBaleares — GOIB_Posidonia_IB, Layer 2 | Public, automated download in `01_data_acquisition.R` |
| Global shoot density | Xarxa de Monitoratge de les Praderies de Posidònia (Govern de les Illes Balears, 2017, 2019, 2022, 2023, 2024) | Public, available in PDF format at [ibanat.caib.es](https://www.ibanat.caib.es) |

## Analysis Structure

### 1. Data Acquisition (`01_data_acquisition.R`)
Downloads the official *Posidonia oceanica* distribution layer from the IDEBaleares ArcGIS REST server using full pagination (32,839 polygons), repairs invalid geometries with `st_make_valid()`, and saves the validated layer in `data/raw/`.

### 2. DHD Index Calculation (`02_calculate_dhd.R`)
Processes daily SST NetCDF files from Copernicus to calculate the Degree-Heating-Days (DHD) index, summing daily thermal anomalies above 28.0 °C during the summer growing season (May–September) for each available year. The result is a multiband raster masked to the actual distribution of the meadow.

### 3. GAMM Calibration (`03_gamm_calibration.R`)
Calibrates the relationship between DHD and changes in shoot density through a tournament of four GAMMs (Gaussian and Scaled-t, with and without depth as a covariate). The winning model (Scaled-t with depth) is validated using spatial cross-validation (Leave-One-Island-Out) and residual diagnostics. Management thresholds are extracted dynamically from the model curve.

### 4. Visualisation and Cartography (`04_visualization_and_mapping.R`)
Generates five cartographic concepts for each analysed year, all masked to the actual distribution of *Posidonia oceanica*. Each concept serves a specific analytical or communication purpose:

- **Concept 1: Categorical Risk Map.** Classifies the marine space based solely on the critical threshold (5.6 DHD). It provides a binary view (safe vs. risk of collapse), useful for identifying areas requiring immediate intervention.
- **Concept 2: Distance to Critical Threshold Map.** Displays the remaining thermal margin (in DHD) before reaching the biological deficit point. It highlights areas that, whilst not yet in regression, are dangerously close to the limit.
- **Concept 3: Predictive Impact Map (GAMM).** Projects the expected interannual change in shoot density across the archipelago, based on the calibrated GAMM. It translates thermal stress directly into expected biological loss (or gain).
- **Concept 4: Discrete Management Traffic Light.** Divides the space into three distinct zones (Green, Yellow, Red) based on the alert and critical thresholds. It offers a simplified, highly communicable overview for policymakers and the general public.
- **Concept 5: Continuous Management Traffic Light (Pure DHD).** Maintains the three-colour traffic light logic but uses a continuous gradient within each category to reflect true spatial variability. This is the most scientifically robust and visually informative option, recommended for publication.

## Initial Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/[username]/posidonia-dhd-baleares.git
   ```
2. **Open the project:** Open `posidonia-dhd-baleares.Rproj` in RStudio.
3. **Install dependencies:**
   ```r
   install.packages(c("terra", "sf", "ggplot2", "dplyr", "tidyr", "mgcv",
                      "tidyterra", "writexl", "readxl", "scales", "httr",
                      "jsonlite", "gridExtra", "svglite"))
   ```
4. **Download SST data:** Execute the `copernicusmarine subset` command provided in the comments of `01_data_acquisition.R` to download temperature data from the Copernicus Marine Service.

## System Requirements

- R (version 4.1 or higher)
- RStudio (recommended)
- Git
- Free account on [Copernicus Marine Service](https://marine.copernicus.eu) for SST data download

### Main R Packages
`terra`, `sf`, `mgcv`, `ggplot2`, `tidyterra`, `dplyr`, `httr`, `jsonlite`, `scales`, `gridExtra`

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Contact

Yeray González Marrero · yeraygma@gmail.com · ORCID: [0000-0003-0306-5208](https://orcid.org/0000-0003-0306-5208)
