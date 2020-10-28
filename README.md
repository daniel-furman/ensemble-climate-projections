## Species Distribution Modeling: Ensemble Climate Projections for Southwestern Deserts

*How can models of climate change effects on geographic distributions incorporate symbiotic species relationships? * What is the best quantification of uncertainty for climate forecasts in Southwest deserts?*

---

All code and data required to reproduce analyses presented at the SICB 2021 conference and SCCUR 2019 conference, as well as Supporting Information.


### Workflow

---

The programming workflow chronologically followed `CMIP6_preprocessing.py` -> `shapefile_preprocessing.R` -> `ML_sdms_train.py` -> `ML_sdms_predict.py ` -> `sdms_2020.R`. The two `ML_sdms_.py` files perform model selection with PyCaret and sk-learn, with results printed and further analyzed in `Comparing_MLs.ipynb`. Geospatial predictions with the best resulting models was then performed in `sdms_2020.R`. The two pre-processing files efficiently piped climate data into the GitHub repository linked below - which should be cloned to re-run the workflow (<1 GB).


### Figures from SICB 2021

---

Figures: To be filled

---

### Data

---
Data is housed in its own [GitHub repository](https://github.com/daniel-furman/xantusia-data). Climate data obtained from [Worldclim version 2](https://www.worldclim.org/) and presence data from GBIF. Future climate forecasts obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html), while the train_tifs represent near current conditions, averaged from 1970-2000. The raster data are on a 2.5x2.5 minutes grid and have the following extent: (-125.0208, -92.00083, 20, 46.9975).

### Abstract

---

The desert night lizard (*Xantusia vigilis*) is a habitat specialist abundantly spread across arid regions of the North American southwest, often reliant on Joshua tree branches (*Yucca brevifolia*) for shelter. Future climate change impacts on *Y. brevifolia* are thus of particular concern for the *X.vigilis'* ecological conservation. 

Here, we explored the impacts of climate change on their geographical distributions by classifying the species' climatic niches. We trained the Species Distribtion Models (SDMs) with a set of ten uncorrelated WorldClim Bioclimatic variables (1970-2000 averages) and presence locations of each species (>1000 unique locations). 

A random forest classifier performed best from a set of over fifteen candidates (including Maxent), emerging as the most predictive model of the current geographic distribution (e.g., OOB misclassification error ~ 2-3%). We then projected the SDM to future climate conditions, simulated with eight climate models from CMIP6 over four Shared Socioeconomic Pathways, for the years 2040-2100. Under these scenarios, the range of *X. vigilis* was predicted to decline to between 11% to 55% of its current distribution, assuming little or no dispersal, overlapped with projections of *Y. brevifolias'* distribution. In addition, a single climate model, CanESM5, consistently predicted the direst scenario of future habitat suitability. 

Our results highlight the importance of including symbiotic and other ecologically important species into models of climate change effects on geographic distributions, with conservation risks possibly heightened for localities which face extreme climate events, such as wildfires.  
