# Testing a general theory for flowering time shift as a function of growing season length

## John S. Park, John Jackson, Anna Bergsten, Jon Ågren


#### 2024-11-05
#### Repository created by John Jackson and John S. Park

---

This repository is a reproducible analysis pipeline for the manuscript *Testing a general theory for flowering time shift as a function of growing season length*. It contains all data used in analysis, 
as well as scripts for the full analysis pipeline used in the study. The associated archived Zenodo repository associated with can be found at the doi
above. Below, please find full file descriptions:

### 1. `scripts/`

The `scripts/` folder contains the full analysis pipeline for the current study. Files include:

* `A1_lythrum_data_aquisition.R` - raw data processing for *Lythrum salicaria*
* `A2_solidago_data_aquisition.R` - raw data processing for *Solidago virgaurea*
* `B_manuscript_schematics.R` - plotting script to generate schematics used in Figure 1 and Figure 2.
* `C_optimal_flowering_models.R` - Bayesian non-linear regression framework to explore the relationship between growing season length and flowering time, present models, and produce Figure 3.
* `D_insitu_flowering_estimates.R` - Using information on the start of spring to predict changes in *in situ* optimal flowering time between different populations, and forecast with increases to the growing season, Figure 4.
* `S1_spatial_data.R` - Using geo-referenced information from source populations to show spatial distribution of the raw data
* `S2_historical_estimates_SoS_GSL.R` - Using site-specific climate data to estimate the growing season length and start of spring across all populations through time.

### 2. `data/`

This folder contains all data used in the current study. Files include:

* `lythrum_phenology_99.csv` - raw phenological data of flowering time from populations of *Lythrum salicaria*
* `PopCoordinatesGatheredMay2023.csv` - geo-location information for populations of *Lythrum salicaria*
* `SMHI_pthbv_t_1961_2022_daily_4326-lythrumpops.csv` - Full local climate data for Umea, the site for *Lythrum salicaria*
* `umea_commongarden_1999_2001_daily_4326.csv` - Site/year specific climate information for Umea, the site for *Lythrum salicaria*
* `lythrum_99.RData` - clean phenological data for *Lythrum salicaria*
* `lythrum_population_coordinates.RData` - spatially explicit data object for *Lythrum salicaria*
* `lythrum_seasoncaps_predict.list.RData` - output from growing season length calculations for populations of *Lythrum salicaria*
* `umea_predict.list.RData` - output from growing season length calculations for Umea, for *in situ* calculations of *Lythrum salicaria*
* `Solidago_FloweringStartUppsalaGarden2005.csv` - raw phenological data of flowering time from populations of *Solidago virgaurea*
* `Solidago_pops_GPS-coordinates.csv` - geo-location information for populations of *Solidago virgaurea*
* `SMHI_pthbv_t_1961_2022_daily_4326-Solidagopops.csv` - Full local climate data for Uppsala, the site for *Solidago virgaurea*
* `uppsala_commongarden_1961_2005_daily_4326.csv` - Site/year specific climate information for Uppsala, the site for *Solidago virgaurea*
* `solidago_05.RData` - clean phenological data for *Solidago virgaurea*
* `solidago_population_coordinates.RData` - spatially explicit data object for *Solidago virgaurea*
* `solidago_seasoncaps_predict.list.RData` - output from growing season length calculations for populations of *Solidago virgaurea*
* `uppsala_predict.list.RData` - output from growing season length calculations for Uppsala, for *in situ* calculations of *Lythrum salicaria*
* `Natural_Earth_Land_Data/` - folder of open source spatial data from https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-land/

### 3. `output/`

This folder contains all output figures used in the manuscript.

