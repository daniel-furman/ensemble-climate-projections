## Ensemble Climate Projections 

---

All code and data required to reproduce analyses presented at the SICB 2021 conference and SCCUR 2019 conference, as well as Supporting Information.

* How can models of climate change effects on geographic distributions incorporate symbiotic species relationships? 
* What is the best quantification of uncertainty for climate forecasts in Southwest deserts?

## Workflow

---

The programming workflow follows the path `CMIP6_preprocessing.py` -> `shapefile_preprocessing.R` -> `ML_sdms_train.py` -> `ML_sdms_predict.py ` -> `sdms_2020.R`. The two ML_ files perform model selection with PyCaret and sk-learn, analyzed in `Comparing_MLs.ipynb`. Geospatial predictions with the resulting models are then performed in `sdms_2020.R`. The two pre-processing files efficiently piped climate data into the GitHub repository linked below - which should be cloned to re-run the workflow (<1 gb).


## Figures from SICB 2021

---

Figures: To be filled

---

## Data

---
Data is housed in its own [GitHub repository](https://github.com/daniel-furman/xantusia-data). Climate data obtained from [Worldclim version 2](https://www.worldclim.org/) and presence data from GBIF. Future climate forecasts obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html), while the train_tifs represent near current conditions, averaged from 1970-2000. The raster data are on a 2.5x2.5 minutes grid and have the following extent: (-125.0208, -92.00083, 20, 46.9975).
