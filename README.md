## Species Distribution Modeling: Ensemble Climate Projections for Southwestern Deserts

*How can models of climate change effects on geographic distributions incorporate symbiotic species relationships? What is the best quantification of uncertainty for climate forecasts in Southwestern deserts?*

---

All code and data required to reproduce research presented at the SICB 2021 and [SCCUR 2019](https://drive.google.com/drive/u/0/folders/15nZUMuGLiINuhSuP6DJ6hg27YKZxeC9A) conferences. See [`Comparing_MLs.ipynb`](https://nbviewer.jupyter.org/github/daniel-furman/ensemble-climate-projections/blob/main/Comparing_MLs.ipynb) first. Spatial predictions were powered by the pyimpute library, see my [contributions](https://github.com/perrygeo/pyimpute/pull/21) merged at v.0.2 (_main.py). 

### Workflow

---

The `ML_sdms_.py` scripts train and validate ML classifiers using the PyCaret library. The outputs are available in `Comparing_MLs.ipynb` alongside geospatial predictions and future suitability ensemble analyses. Predictions were checked by also performing the predcitions in RStudio, with the traditionally used dismo package, in `sdms_2020.R`. The two pre-processing scripts (`CMIP6_preprocessing.py` and `shapefile_preprocessing.R`) efficiently piped geospatial data to its own [repository](https://github.com/daniel-furman/xantusia-data) which can be cloned before re-running the workflow (<1 GB). Lastly, `_main.py` fixed deprecated code within the pyimpute library, critical to performing geospatial classification in Python, which was [merged](https://github.com/perrygeo/pyimpute/pull/21) to the master branch for version 0.2.

### Data

---
All data required for the analyses is contained in a separate [repository](https://github.com/daniel-furman/xantusia-data). Climate data obtained from [Worldclim v. 2.1](https://www.worldclim.org/) and presence data from GBIF. Future climate forecasts obtained from [CMIP6](https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html), while the train_tifs represent near current conditions, averaged from 1970-2000. The raster data are on a 2.5 arc-minute grid and have the following extent: (92°W to 125°W and 20°N, 47°N).

### Select Writings

---

#### 1 | Abstract
The desert night lizard (*Xantusia vigilis*) is a habitat specialist abundantly spread across arid regions of the North American southwest, often reliant on Joshua tree branches (*Yucca brevifolia*) for shelter. Predictions of range restriction for *Y. brevifolia* are thus of particular concern for the night lizard’s ecological conservation. 

Here, we explored climate change impacts on *X. vigilis’* geographical distribution by associating species presences to climatic niche, aiming to relate the night lizard’s responses to those for *Y. brevifolia*. We employed ten decorrelated Bioclimatic variables (1970-2000 averages) with sets of presence locations for each species (>1000 unique locations). Five classifiers were selected from a set of over twenty candidates (including Maxent) (validation-set misclassification error ~ 3%), and spatially blended to obtain a regionally accurate current distribution. Future climate conditions were simulated with eight climate models (CMIP6) across twelve temporal and SSP scenarios. 

Assuming little to no dispersal, we found the climatic niche of *X. vigilis* declining from ~55% to ~11% of its current range and *Y. brevifolia* similarly declining. We also measured spatial similarity between projections for the two species, as evidenced by overlapping PCA transformations across the predictions of future climatic niche. Our results highlight the importance of including symbiotic and other ecologically important species into models of climate change effects on geographic distributions, with conservation risks possibly heightened for local populations which face extreme climate events, such as wildfires.

#### 2 | Methods

#### 2.1 | Spatial Data 
Spatial analyses were performed in the R statistical environment (v. 4.0.3; R Core Team, 2020) with the Raster (v. 3.3-13; Hijmans et al., 2012), Geospatial Data Abstraction Library (“rgdal”; v. 1.5-18; Bivand, Keitt, & Rowlingson, 2019), and Spatial Points packages (“sp”; Pebesma & Bivand., 2005). 
				
#### 2.1.1 | Raster Data
We downloaded 19 bioclimatic rasters (2.5 arc-minute resolution; 1970-2000) from the publicly available WorldClim database (www.worldclim.org; v. 2.1, Fick & Hijmans, 2017), which contain ecologically meaningful axes of species’ climatic niche commonly used for SDMs. The spatial extent of the analyses was 92°W to 125°W and 20°N, 47°N, cropped so to include the entire latitudinal and longitudinal distribution of X. vigilis and Y. brevifolia and a background sampling space suited for projection. To avoid high levels of collinearity among model features, we selected the ten most relevant bioclimatic variables from Maxent’s jackknife resampled importance measures, filtered below a Spearman’s correlation threshold (r < 0.5) across the entire modelled environmental space. 

#### 2.1.2 | Presence Data
Occurrence records were obtained from the publicly available Global Biodiversity Information Facility database (GBIF; www.gbif.org, downloaded November 1, 2020). We processed occurrence records to remove duplicates (within 1e-4 decimal degrees), erroneous (belonging to a different subspecies), and/or spatially biased (grid sampled) coordinates. The final sets of presences contained at most 25 presences per raster pixel. Background points, also known as pseudo-absences, were sampled at random to twice the number of species presences across the default extent factor (5%).

#### 2.2 | Species Distribution Modeling
Modeling was conducted using the Python programming language (v. 3.7.6) with the SciKit-Learn (v. 0.23.2;  Pedregosa et al., 2011), PyCaret (v. 2.2; Moez, 2019) and PyImpute libraries (v. 0.2 ; Perry, 2015), using R with the dismo package (v. 1.1-4 ; site), and using the Maxent software (v. 3.4.1; Phillips et al., 2006; Phillips & Dudik, 2008). The covariate data were first split into training (80%) and validation (20%) sets to assess model stability and sensitivity; however, final models were built using 100% of the available data, definitively preserving small subgroups of presence data (e.g. the isolated X. vigilis Pinnacles clade). Five classifiers were selected from a set of over twenty candidates: Random Forest (rf), Extremely Randomized Trees (et), Light Gradient Boosting Machine (lightgbm), Extreme Gradient Boosting (xgboost), and Catboost Gradient Boosting (catboost), trained over ten partitions of random selections with replacement (10-fold cross validation), with a validation set performance of {F score > 0.939; AUC > .991, false negatives ~3%} . Spatial outputs of binary presence/absence were then blended so as to only include areas where all five classifiers predicted suitability, resulting in a prediction that best matched previous estimates of the species’ geospatial distribution.

#### 2.3 | Future Distribution Projections
Models for X. vigilis and Y. Brevifolia were projected onto future climate , including four greenhouse gas scenarios, or shared socio-economic pathways (SSPs) for each of eight global climate models (GCMs) (CMIP6, Eyring et al., 2016), between two decade periods across the 21st century (2040-2100). The GCMs included BCC-CSM2-MR, CanESM5, CNRM-CM6-1, CNRM-ESM2-1, IPSL-CM6A-LR, MIROC-ES2L, MIROC6, and MRI-ESM2-0. The final consensus prediction used 480 models (spatial blend of five algorithms, eight GCMs, three time periods, and four SSPs). We then used an ensemble approach to identify levels of agreement among the projections, estimating the proportion of geographic range that remains climatically suitable with negligible species dispersal to novel areas. 
