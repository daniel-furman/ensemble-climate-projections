## Ensemble Climate Projections 

---

All code and data required to reproduce analyses presented at the SICB 2021 conference and SCCUR 2019 conference, as well as the Supporting Information.

* How can models of climate change effects on geographic distributions incorporate symbiotic species relationships? 
* What is the best quantification of uncertainty for climate forecasts in Southwest deserts?

## Workflow

---

The programming workflow follows the path `CMIP6_preprocessing.py` -> `shapefile_preprocessing.R` -> `ML_sdms_train.py` -> `ML_sdms_predict.py ` -> `sdms_2020.R`. The two ML_ files perform model selection with PyCaret and sk-learn, analyzed in `Comparing_MLs.ipynb`. Geospatial predictions with the resulting models are then performed in `sdms_2020.R`, using the [dismo](https://cran.r-project.org/web/packages/dismo/dismo.pdf) package. The first two pre-processing files were employed to create the data repository below, yet are not intended to be re-ran as is. 


## Figures from SICB 2021

---

Figures: To be filled

---

## Data

---

All data 
