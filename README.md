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

**Question 1: *What underlying uncertainties are contained in geospatial climate change forecasts?***  An ensemble of SDMs were extrapolated to unseen data, across eight Global Climate Models (GCMs), four shared socioeconomic pathways, and three bi-decade time periods. Across these forecasts, we predicted roughly equal magnitudes of habitat degradation for the study species (49% to 90.7% Night Lizard decline from baseline; 52% to 91.4% Joshua Tree decline from baseline). We explored uncertainty in these projections through spatial overlap among the eight GCMs, assuming negligible species dispersal.


**Question 2: *How can models of climate change effects on geographic distributions best incorporate symbiotic species relationships?*** We first benchmarked climatic change between interpolation and extrapolation data, with Jaccard Similarity among principal components. With the models thus benchmarked, we help confirm our measurements of habitat decline (see above). Furthermore, we predicted the species distributions as spatially diverging across the 21st century, with substantial decline in the pair’s overlap (~56% decrease from baseline climate, on average). By 2090, geographic overlap between the two species may decrease by as much as ~87% from current conditions. 



### Programming Workflow

---

The `ML_sdms_.py` train and then validate ML classifiers. The `recursive-ranker.py` function recursively selected which features to use for the modeling, such that they were below a Spearman's threshold. We used the rank of feature importance scores to decide which variables to drop at each recursive call. These outputs are available in `Comparing_MLs.ipynb`, along with the geospatial predictions for the baseline and future climates. Lastly, `pca_benchmark.R` calculates the similarity between the model interpolation and extrapolation data using a Jacard similarity metric among principal components. 


### Data

---

The input data is located in the `data/` folder. Climate information was stored in 19 bioclimatic features (2.5 arc-minute resolution; baseline 1970-2000; with an extent from 109.3°W to 122.8°W and 31.9°N to 38.2°N), downloaded from the publicly available [WorldClim database](https://www.worldclim.org) (v. 2.1, Fick & Hijmans, 2017). Presence data were downloaded from the publicly available Global Biodiversity Information Facility database ([GBIF](https://www.gbif.org), downloaded November 1, 2020) and from Leavitt et al., 2007. In addition, the most recently updated climate forecasts were obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html), also from the WorldClim database. 

### Conclusion:

---

Our results reveal the importance of symbiotic species relationships for SDMs, so to more confidently select areas within predictions of truly suitable habitat. We hypothesize that habitat degradation for our two study species will be heightened where severely changing climate is paired with environmental catastrophe, such as strong wildfire. We then pinpointed areas across the region to target conservation for the two species.
