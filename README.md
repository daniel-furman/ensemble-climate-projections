## Assessing Climate Change Impacts on Xantusia vigilis lizards and their Joshua tree habitats with Ensemble Species Distribution Models

*How well do a variety of species distribution models (SDMs) perform for current geographic ranges of the two species? How much will climate change constrict the two species’ distributions, as well as their overlap. What can we kearn about jointly projecting SDMs of climate change impacts on these and other symbiotic species.*

---

See the online notebook first: [`Comparing_MLs.ipynb`](https://nbviewer.jupyter.org/github/daniel-furman/ensemble-climate-projections/blob/main/Comparing_MLs.ipynb). All code and data required to reproduce research presented at the SICB 2021 and [SCCUR 2019](https://drive.google.com/file/d/114wmqQgjkc5DHLQmVI19AvlTw4K_daYQ/view?usp=sharing) conferences. Spatial classification was powered by the pyimpute library, see my open-source [contributions](https://github.com/perrygeo/pyimpute/pull/21) therein. Also see my corresponding <a target="_blank" rel="noopener noreferrer" href="https://daniel-furman.github.io/py-sdms-tutorial/"> tutorial</a> on predicting geospatial distributions with machine learning models in Python.

### Discussion 
---

A blend of tree-based, machine learning classifiers performed best at predicting the current geographic ranges for Xant. and Ybrev. The overlap in baseline geographic range was roughly fifty percent (intersection/union). Substantial bioclimatic change is predicted to constrict the geographic ranges for both species (Fig. 6), as well as reduce species overlap by ~56% from baseline conditions, on average. Areas with a high likelihood of future species overlap were then pinpointed, to target ecological conservation for Mojave Desert populations of Xant. reliant on fallen Joshua trees for shelter.



### Programming Workflow

---

* Data pre-processing -> Model fitting -> Assessment -> Baseline interpolation (1970-2000) -> Extrapolation across 21st century

The `ML_sdms_.py` train and then validate ML classifiers. The `recursive-ranker.py` function recursively selected which features to use for the modeling, such that they were below a Spearman's threshold. We used the rank of feature importance scores to decide which variables to drop at each recursive call. These outputs are available in `Comparing_MLs.ipynb`, along with the geospatial predictions for the baseline and future climates. Lastly, `pca_benchmark.R` calculates the similarity between the model interpolation and extrapolation data using a Jacard similarity metric among principal components. 


### Data

---

The input data is located in the `data/` folder. Climate information was stored in 19 bioclimatic features (2.5 arc-minute resolution; baseline 1970-2000; with an extent from 109.3°W to 122.8°W and 31.9°N to 38.2°N), downloaded from the publicly available [WorldClim database](https://www.worldclim.org) (v. 2.1, Fick & Hijmans, 2017). Presence data were downloaded from the publicly available Global Biodiversity Information Facility database ([GBIF](https://www.gbif.org), downloaded November 1, 2020) and from Leavitt et al., 2007. In addition, the most recently updated climate forecasts were obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html), also from the WorldClim database. 

### Conclusion:

---

Our results reveal the importance of symbiotic species relationships for SDMs, so to more confidently select areas within predictions of truly suitable habitat. We hypothesize that habitat degradation for our two study species will be heightened where severely changing climate is paired with environmental catastrophe, such as strong wildfire. We then pinpointed areas across the region to target conservation for the two species.
