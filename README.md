## Species Distribution Modeling: Ensemble Climate Projections for Southwestern Deserts

*How can models of climate change effects on geographic distributions incorporate symbiotic species relationships? What is the best quantification of uncertainty for climate forecasts in Southwestern deserts?*

---

All code and data required to reproduce research presented at the SICB 2021 and [SCCUR 2019](https://drive.google.com/file/d/1Z9HJnW3p1tecLkUbK6zODJH1dMgzz06j/view?usp=sharing) conferences. See the online notebook first: [`Comparing_MLs.ipynb`](https://nbviewer.jupyter.org/github/daniel-furman/ensemble-climate-projections/blob/main/Comparing_MLs.ipynb). Spatial predictions were powered by the pyimpute library, see my open-source [contributions](https://github.com/perrygeo/pyimpute/pull/21) therein. Also see my corresponding <a target="_blank" rel="noopener noreferrer" href="https://daniel-furman.github.io/py-sdms-tutorial/"> tutorial in Python</a>.


### Data

---
Link to [data repository](https://github.com/daniel-furman/xantusia-data). Climate information from 19 bioclimatic features (2.5 arc-minute resolution; 1970-2000; 92°W to 125°W and 20°N, 47°N) was downloaded from the publicly available [WorldClim database](https://www.worldclim.org) (v. 2.1, Fick & Hijmans, 2017). Presence data were downloaded from the publicly available Global Biodiversity Information Facility database ([GBIF](https://www.gbif.org), downloaded November 1, 2020). In addition, the most recently updated climate forecasts were obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html).  

### Workflow

---

The two pre-processing scripts (`CMIP6_preprocessing.py` and `shapefile_preprocessing.R`) efficiently piped the geospatial data (1000's of GB). The `ML_sdms_.py` scripts train and validate ML classifiers using the PyCaret library. The outputs are available in `Comparing_MLs.ipynb` alongside geospatial predictions and future suitability ensemble analyses.
