## Assessing Climate Change Impacts on a Pair of Symbiotic Species with Ensemble Species Distribution Models


*What underlying uncertainties are contained in geospatial climate change forecasts? How can models of climate change effects on geographic distributions incorporate symbiotic species relationships?*

---

See the online notebook first: [`Comparing_MLs.ipynb`](https://nbviewer.jupyter.org/github/daniel-furman/ensemble-climate-projections/blob/main/Comparing_MLs.ipynb). All code and data required to reproduce research presented at the SICB 2021 and [SCCUR 2019](https://drive.google.com/file/d/114wmqQgjkc5DHLQmVI19AvlTw4K_daYQ/view?usp=sharing) conferences. Spatial classification was powered by the pyimpute library, see my open-source [contributions](https://github.com/perrygeo/pyimpute/pull/21) therein. Also see my corresponding <a target="_blank" rel="noopener noreferrer" href="https://daniel-furman.github.io/py-sdms-tutorial/"> tutorial</a> on predicting geospatial distributions with machine learning models in Python.

### Introduction 
---

* Niche: An envelope of abiotic and biotic factors suitable to the survival of a species
* SDMs: Associates presence locations to environmental variables, ideally an estimate of fundamental niche
* Workflow: `Data pre-processing` -> `Model fitting` -> `Assessment` -> `Baseline interpolation (1970-2000)` -> `Extrapolation across time`

**Question 1: *What underlying uncertainties are contained in geospatial climate change forecasts?*** An ensemble of Species Distribution Models were extrapolated across eight Global Climate Models, four shared socioeconomic pathways, and three bi-decade time periods. Across these conditions, we predicted similar decline in suitable habitat for *X. vigilis* (51% to 9.3% of baseline) and *Y. brevifolia* (48% to 8.6% of baseline), considering areas where at least five GCMs overlapped, assuming negligible species dispersal (area intersection / baseline).

<img src="data/ensemble_extrapolation.png" width = 630/>

**Question 2: *How can models of climate change effects on geographic distributions best incorporate symbiotic species relationships?*** We minimized modelling error by using a soft voting ensemble of well-fit classifiers, as well as by benchmarking climatic change between interpolation and extrapolation data, with Jaccard Similarity among principal components. Over the above forecasts, we predicted constricting spatial overlap between the species (~56% decrease from baseline climate, on average), which worsened across time. 

**Conclusion:** Our results reveal the importance of symbiotic species relationships in Species Distribution Models of climate change effects. We hypothesize that habitat degradation will be heightened for areas with both severe change in climate and environmental catastrophe, such as strong wildfire. We next strive to identify the areas most overlapped across time, primarily so to target ecological conservation for local species populations.  

### Programming Workflow

---

The `ML_sdms_.py` train and then validate ML classifiers. The `recursive-ranker.py` function recursively selected which features to use for the modeling, such that they were below a Spearman's threshold. We used the rank of feature importance scores to decide which variables to drop at each recursive call. These outputs are available in `Comparing_MLs.ipynb`, along with the geospatial predictions for the baseline and future climates. Lastly, `pca_benchmark.R` calculates the similarity between the model interpolation and extrapolation data using a Jacard similarity metric among principal components. 


### Data

---

The input data is located in the `data/` folder. Climate information was stored in 19 bioclimatic features (2.5 arc-minute resolution; baseline 1970-2000; with an extent from 109.3°W to 122.8°W and 31.9°N to 38.2°N), downloaded from the publicly available [WorldClim database](https://www.worldclim.org) (v. 2.1, Fick & Hijmans, 2017). Presence data were downloaded from the publicly available Global Biodiversity Information Facility database ([GBIF](https://www.gbif.org), downloaded November 1, 2020) and from Leavitt et al., 2007. In addition, the most recently updated climate forecasts were obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html), also from the WorldClim database. 

