## Assessing Climate Change Impacts on *Xantusia vigilis* lizards and their Joshua tree habitats with Ensemble Species Distribution Models

*How well do a variety of species distribution models (SDMs) perform for current geographic ranges of the two species? How much will climate change shift the spatial overlap between the two species’ distributions? What can we kearn about jointly projecting SDMs of climate change impacts on these and other symbiotic species?*

---

See the online notebook first: [`Comparing_MLs.ipynb`](https://nbviewer.jupyter.org/github/daniel-furman/ensemble-climate-projections/blob/main/Comparing_MLs.ipynb). All code and data required to reproduce research presented at the SICB 2021 and [SCCUR 2019](https://drive.google.com/file/d/114wmqQgjkc5DHLQmVI19AvlTw4K_daYQ/view?usp=sharing) conferences. Spatial classification was powered by the pyimpute library, see my open-source [contributions](https://github.com/perrygeo/pyimpute/pull/21) therein. Also see my corresponding <a target="_blank" rel="noopener noreferrer" href="https://daniel-furman.github.io/py-sdms-tutorial/"> tutorial</a> on predicting geospatial distributions with machine learning models in Python.

### Discussion 
---

Blends of tree-based machine learning classifiers performed best at predicting the geographic range for *Xant.* and *Ybrev.*, a pair of iconic Mojave desert species. The spatial overlap between the two species’ near-current (1970-2000) distributions was ~fifty percent for overall intersection/union and ~x percent solely within *Xant.’s* geographic range. Across a multitude of future bioclimatic scenarios, we predicted severe constriction for both species' distributions (Fig. 6) and a significant reduction in their geographic overlap (down by ~56% from baseline conditions, on average). Equipped with these predictions across the 21st century, we pinpoint areas with a high likelihood of retaining *Xant.* and *Ybrev.* overlap and suggest ways to target their ecological conservation.

Our research can serve as a case study for projecting joint distributions of symbiotic species to explore ecological climate change affects. Habitat degradation for Mojave Desert species will likely result from changing climate and from discrete events of environmental catastrophe, such as wildfire and drought. Future inquiry will include an assessment of lag between *Xant.* and *Ybrev.* range shifts, movement in topography among predictions, and field work at Joshua tree strands recently burned by wildfire (see [NYT article](https://www.nytimes.com/interactive/2020/12/09/climate/redwood-sequoia-tree-fire.html?) near the southern California/Arizona border.



### Programming Workflow

---

* Data pre-processing -> Model fitting -> Assessment -> Baseline interpolation (1970-2000) -> Extrapolation across 21st century

The `ML_sdms_.py` train and then validate ML classifiers. The `recursive-ranker.py` function recursively selected which features to use for the modeling, such that they were below a Spearman's threshold. We used the rank of feature importance scores to decide which variables to drop at each recursive call. These outputs are available in `Comparing_MLs.ipynb`, along with the geospatial predictions for the baseline and future climates. Lastly, `pca_benchmark.R` calculates the similarity between the model interpolation and extrapolation data using a Jacard similarity metric among principal components. 


### Data

---

The input data is located in the `data/` folder. Climate information was stored in 19 bioclimatic features (2.5 arc-minute resolution; baseline 1970-2000; with an extent from 109.3°W to 122.8°W and 31.9°N to 38.2°N), downloaded from the publicly available [WorldClim database](https://www.worldclim.org) (v. 2.1, Fick & Hijmans, 2017). Presence data were downloaded from the publicly available Global Biodiversity Information Facility database ([GBIF](https://www.gbif.org), downloaded November 1, 2020) and from Leavitt et al., 2007. In addition, the most recently updated climate forecasts were obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html), also from the WorldClim database. 

