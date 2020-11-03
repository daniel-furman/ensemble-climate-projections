## Species Distribution Modeling: Ensemble Climate Projections for Southwestern Deserts

*How can models of climate change effects on geographic distributions incorporate symbiotic species relationships? What is the best quantification of uncertainty for climate forecasts in Southwestern deserts?*

---

All code and data required to reproduce research presented at the SICB 2021 and [SCCUR 2019](https://drive.google.com/drive/u/0/folders/15nZUMuGLiINuhSuP6DJ6hg27YKZxeC9A) conferences. See [`Comparing_MLs.ipynb`](https://nbviewer.jupyter.org/github/daniel-furman/ensemble-climate-projections/blob/main/Comparing_MLs.ipynb) first. Spatial predictions were powered by the pyimpute library, see my [contributions](https://github.com/perrygeo/pyimpute/pull/21) merged at version 0.2 (_main.py). 

### Workflow

---

The `ML_sdms_.py` scripts train and validate ML classifiers using the PyCaret library. The outputs are available in `Comparing_MLs.ipynb` alongside geospatial predictions and future suitability ensemble analyses. Predictions were checked by also performing the predcitions in RStudio, with the traditionally used dismo package, in `sdms_2020.R`. The two pre-processing scripts (`CMIP6_preprocessing.py` and `shapefile_preprocessing.R`) efficiently piped geospatial data to its own [repository](https://github.com/daniel-furman/xantusia-data) which can be cloned before re-running the workflow (<1 GB). Lastly, `_main.py` fixed deprecated code within the pyimpute library, critical to performing geospatial classification in Python, which was [merged](https://github.com/perrygeo/pyimpute/pull/21) to the master branch for version 0.2.

### Data

---
All data required for the analyses is contained in a separate [repository](https://github.com/daniel-furman/xantusia-data). Climate data obtained from [Worldclim version 2](https://www.worldclim.org/) and presence data from GBIF. Future climate forecasts obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html), while the train_tifs represent near current conditions, averaged from 1970-2000. The raster data are on a 2.5 arc-minute grid and have the following extent: (92°W to 125°W and 20°N, 47°N).


### Abstract

---

The desert night lizard (*Xantusia vigilis*) is a habitat specialist abundantly spread across arid regions of the North American southwest, often reliant on Joshua tree branches (*Yucca brevifolia*) for shelter. Future climate change impacts on *Y. brevifolia* are thus of particular concern for the *X.vigilis'* ecological conservation. 

Here, we explored the impacts of climate change on their geographical distributions by classifying the species' climatic niches. We trained the Species Distribution Models (SDMs) with a set of ten uncorrelated WorldClim Bioclimatic variables (1970-2000 averages) and presence locations of each species (>1000 unique locations). 

A random forest classifier performed best from a set of over fifteen candidates (including Maxent), emerging as the most predictive model of the current geographic distribution (e.g., OOB misclassification error ~ 3%). We then projected the SDM to future climate conditions, simulated with eight climate models from CMIP6 over four Shared Socioeconomic Pathways, for the years 2040-2100. Under these scenarios, the range of *X. vigilis* was predicted to decline to between 11% to 55% of its current distribution, assuming little or no dispersal, overlapped with projections of *Y. brevifolias'* distribution. In addition, a single climate model, CanESM5, consistently predicted the direst scenario of future habitat suitability. 

Our results highlight the importance of including symbiotic and other ecologically important species into models of climate change effects on geographic distributions, with conservation risks possibly heightened for localities which face extreme climate events, such as wildfires.  

### Requirements

---

Python dependencies are listed in a `requirements-py.txt` file, including the library version numbers. You can replicate the environment your codebase needs by using virtualenv:

```
# This creates the virtual environment
cd $PROJECT-PATH
virtualenv ensemble-climate-projections

# Then install the dependencies by referring to the requirements-py.txt:

# This installs the modules
pip install -r requirements-py.txt

# This activates the virtual environment
source ensemble-climate-projections/bin/activate
```
R dependencies are listed in a `requirements-R.txt` file. You can replicate the environment your codebase needs by using install.packages():

```
reqs <- read.table('PROJECT-PATH/requirements-R.txt', col.names = c('package', 'version'))
for (i in 1:length(reqs$package)){
  install.packages(reqs$package[i]) # devtools::install_version if versions desired
  }
```

### Poster Figures (SICB 2021)

---

Figures: To be filled closer to the conference date. 

---
