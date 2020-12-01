## Assessing Climate Change Impacts on a Pair of Symbiotic Species with Ensemble Species Distribution Models


*What underlying uncertainties are contained in geospatial climate change forecasts? How can models of climate change effects on geographic distributions incorporate symbiotic species relationships?*

---

See the online notebook first: [`Comparing_MLs.ipynb`](https://nbviewer.jupyter.org/github/daniel-furman/ensemble-climate-projections/blob/main/Comparing_MLs.ipynb). All code and data required to reproduce research presented at the SICB 2021 and [SCCUR 2019](https://drive.google.com/file/d/114wmqQgjkc5DHLQmVI19AvlTw4K_daYQ/view?usp=sharing) conferences. Spatial classification was powered by the pyimpute library, see my open-source [contributions](https://github.com/perrygeo/pyimpute/pull/21) therein. Also see my corresponding <a target="_blank" rel="noopener noreferrer" href="https://daniel-furman.github.io/py-sdms-tutorial/"> tutorial</a> on predicting geospatial distributions with machine learning models in Python.

### Introduction 
---

* `Fundamental niche:` The entire envelope of abiotic and biotic factors suitable to the survival of a species
* `Species Distribution Models (SDMs):` Associates presence locations to environmental variables, an estimate of fundamental niche
* `Symbiotic species:` A close and long-term biological interaction between two species, e.g. the Desert Night Lizard (*X. vigilis*) is often reliant on Joshua Trees (*Y. brevifolia*) for shelter
* `Workflow:` Data pre-processing -> Model fitting -> Assessment -> Baseline interpolation (1970-2000) -> Extrapolation across 21st century

**Question 1: *What underlying uncertainties are contained in geospatial climate change forecasts?*** An ensemble of Species Distribution Models were extrapolated to unseen data, across eight Global Climate Models, four shared socioeconomic pathways, and three bi-decade time periods. Across these future climate forecasts, we predicted the same magnitude of habitat degredation for the two study species (49% to 90.7% Night Lizard decline from baseline; 52% to 91.4% Joshua Tree decline from baseline). For the figure below, we consider the mean where at least five/eight GCMs were in agreement, per scenario, assuming negligible species dispersal from the baseline distribution (area of intersection / area of baseline).

<p align="center"><img src="data/ensemble_extrapolation.png" width = 630/>

**Question 2: *How can models of climate change effects on geographic distributions best incorporate symbiotic species relationships?*** We minimized modelling error by using a soft voting ensemble of well-fit classifiers, as well as by benchmarking climatic change between interpolation and extrapolation data, with Jaccard Similarity among principal components. While the magnitude of future habitat decline was roughly equivalent for the two species (as per above), the distributions were predicted to diverge, with declining overlap between the two species across the 21st century (~56% decrease from current conditions, on average). By 2090, spatial overlap between the two distributions may decrease by as much as ~87% of current conditions. 

**Conclusion:** Our results reveal the importance of symbiotic species relationships in Species Distribution Models to more confidently pinpoint suitable areas for conservation, helping to protect the study area's deserts in the face of climate change. We hypothesize that habitat degradation will be heightened for areas with severe change in climate paired with environmental catastrophe, such as strong wildfire.

### Programming Workflow

---

The `ML_sdms_.py` train and then validate ML classifiers. The `recursive-ranker.py` function recursively selected which features to use for the modeling, such that they were below a Spearman's threshold. We used the rank of feature importance scores to decide which variables to drop at each recursive call. These outputs are available in `Comparing_MLs.ipynb`, along with the geospatial predictions for the baseline and future climates. Lastly, `pca_benchmark.R` calculates the similarity between the model interpolation and extrapolation data using a Jacard similarity metric among principal components. 


### Data

---

The input data is located in the `data/` folder. Climate information was stored in 19 bioclimatic features (2.5 arc-minute resolution; baseline 1970-2000; with an extent from 109.3°W to 122.8°W and 31.9°N to 38.2°N), downloaded from the publicly available [WorldClim database](https://www.worldclim.org) (v. 2.1, Fick & Hijmans, 2017). Presence data were downloaded from the publicly available Global Biodiversity Information Facility database ([GBIF](https://www.gbif.org), downloaded November 1, 2020) and from Leavitt et al., 2007. In addition, the most recently updated climate forecasts were obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html), also from the WorldClim database. 

