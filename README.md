## Assessing Climate Change Impacts on a Pair of Symbiotic Species with Ensemble Species Distribution Models


*How can models of climate change effects on geographic distributions incorporate symbiotic species relationships? What underlying uncertainties are contained in geospatial climate change forecasts?*

---

See the online notebook first: [`Comparing_MLs.ipynb`](https://nbviewer.jupyter.org/github/daniel-furman/ensemble-climate-projections/blob/main/Comparing_MLs.ipynb). All code and data required to reproduce research presented at the SICB 2021 and [SCCUR 2019](https://drive.google.com/file/d/114wmqQgjkc5DHLQmVI19AvlTw4K_daYQ/view?usp=sharing) conferences. Spatial classification was powered by the PyImpute library, see my open-source [contributions](https://github.com/perrygeo/pyimpute/pull/21) therein. Also see my corresponding <a target="_blank" rel="noopener noreferrer" href="https://daniel-furman.github.io/py-sdms-tutorial/"> Python tutorial</a>.


### Data

---
The input data is located in the `data/` subfolder. Climate information from 19 bioclimatic features (2.5 arc-minute resolution; baseline 1970-2000; with an extent from 109.3°W to 122.8°W and 31.9°N to 38.2°N) was downloaded from the publicly available [WorldClim database](https://www.worldclim.org) (v. 2.1, Fick & Hijmans, 2017). Presence data were downloaded from the publicly available Global Biodiversity Information Facility database ([GBIF](https://www.gbif.org), downloaded November 1, 2020) and Leavitt et al., 2007. In addition, the most recently updated climate forecasts were obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html), these are contained in a seperate [data repository](https://github.com/daniel-furman/xantusia-data).

### Workflow

---

The three pre-processing scripts (`CMIP6_preprocessing.py`, `shapefile_preprocessing.R`, and `SDMs_data_piping.R`) ETL'ed the geospatial data (1000's of GB). The `ML_sdms_.py` scripts train and validate ML classifiers using the PyCaret library. Modelling also relied on `recursive-ranker.py`, function recursively selects de-correlated features for modeling below a Spearman's metric threshold, using the rank of feature importance scores. These outputs are available in `Comparing_MLs.ipynb` with geospatial predictions and future suitability ensemble analyses. Lastly, `pca.R` is a rought draft of the script that will find Jaccard similarity among 3d hulls for PCA features of climatic niche.
