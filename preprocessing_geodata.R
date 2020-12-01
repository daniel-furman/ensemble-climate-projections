# ---
# title: "sdms_data_piping"
# author: "Daniel Furman"
# date: "2020"
# ---

#' 
#' This R script pipes presence and absence data for xantusia
#' vigilis and yucca brevifolia from the raw Worldclim and 
#' GBIF data.
#' 
# X. vigilis citation:
# GBIF.org (01 November 2020) GBIF Occurrence
# Download https://doi.org/10.15468/dl.ceutzd 

# Y. brevifolia citation:
# GBIF.org (01 November 2020) GBIF Occurrence
# Download https://doi.org/10.15468/dl.g6swrm 


## ------------------------------------------------------------------------

# Libraries needed for this file, download all dependencies from 
# requirements-R.txt

library(raster)
library(rgdal)
library(rgeos)
library(maptools)
library(USAboundaries)
library(sf)
library(dismo)

set.seed(100) # seed using random number generator

## ------------------------------------------------------------------------

# Raster data, for the current distribution 

# set study space extent:
ext <- extent(-122.767, -109.3192, 31.91825, 38.1922)
# download raw bioclimatic tifs first from WorldClim 2.1
# grab those at the 2.5 arcminutes 
train.files <- list.files(
  path='data_2.0/raw/wc2.1_2.5m_bio', full.names=TRUE)
train.layers <- stack(train.files)
train.layers <- crop(train.layers, ext)
bcvar.index <- c(1,10:19, 2:9)
for (i in c(1:19)){
  names(train.layers[[i]])<-paste('bclim', bcvar.index[i], sep = "")
}
names(train.layers)
print(train.files)
#dir.create(path = "data/train-rasters-2.5m/")
#writeRaster(train.layers, filename=paste('data_2.0/train-rasters-2.5m/',
#names(train.layers), sep = ""), bylayer=TRUE,format="ascii")


## ------------------------------------------------------------------------

# Xantusia vigilis:
# download GBIF database:
xv_raw <- read.csv( 'data_2.0/raw/xv_GBIF_presences.csv', header = TRUE)
xv <- data.frame(matrix(ncol = 2, nrow = length(xv_raw$decimalLongitude)))
xv[,1] <- xv_raw$decimalLongitude
xv[,2] <- xv_raw$decimalLatitude
colnames(xv) <- c('lon','lat')
xv.duplicates <- duplicated(xv)
xv <- xv[!xv.duplicates,]
# next load the Leavit et al., 2007 database
xv_raw_Leavitt <- read.csv( 'data_2.0/raw/Leavitt_2007.csv', header = FALSE)
xv_raw_Leavitt <- xv_raw_Leavitt[,1:2]
colnames(xv_raw_Leavitt) <- c('lon','lat')

# merge the two sets of Xantusia presences, filter one remaining duplicate
xv <- rbind(xv_raw_Leavitt, xv)
xv.duplicates <- duplicated(data.frame(xv$lon, xv$lat))
xv <- xv[!xv.duplicates,]
xv <- xv[complete.cases(xv),]

# remove subspecies of other xantusiid subspecies
# first, take only those within study extent: 
xv <- xv[which(
  xv$lon>=ext[1] & xv$lon<=ext[2] & xv$lat>=ext[3] & xv$lat<=ext[4]),] 
# remove Yucca Valley and San Jacinto clades:
xv <- xv[-which(
  xv$lon>=-117 & xv$lon<=-115.5 & xv$lat>=33.25 & xv$lat<=34.25),]
xv <- xv[-which(
  xv$lon>=-116.8 & xv$lon<=-116 & xv$lat>=34 & xv$lat<=34.5),]
# remove X. wigginsi:
xv <- xv[-which(xv$lon>=-118 & xv$lon<=-114.2 & xv$lat<=32.8),]
# remove X. arizonae and bezyi:
xv <- xv[-which(xv$lon>=-113.2 & xv$lat<=34.6),]

# print the final extent (5 % larger than finalized xv presences)
print(extent(min(xv$lon), max(xv$lon), min(
  xv$lat),max(xv$lat)) * 1.1)
ext = extent(min(xv$lon), max(xv$lon), min(xv$lat),max(xv$lat)) * 1.1

# gridded sampling
r <- raster(ext, res = res(train.layers), crs = crs(train.layers))
acsel <- gridSample(data.frame(xv$lon, xv$lat), r, n=2)
head(acsel)
length(acsel$xv.lon) # 907 after down from 1187 before
colnames(acsel) <- c('lon', 'lat')

# sample background points from a slightly wider extent
bg <- randomPoints(train.layers[[1]], length(acsel$lon),
                   ext=ext, extf = 1) 
colnames(bg) <- c('lon','lat')
bg <- data.frame(bg)
bg <- gridSample(data.frame(bg$lon, bg$lat), r, n=1)
head(bg)
colnames(bg) <- c('lon','lat')
length(bg$lon) #907

# plot our points to check
map_states <- us_boundaries("2000-12-31", states = c("California",
                                                     "Arizona",
                                                     "Nevada", "Utah"))
col2rgb("grey29")
mycol <- rgb(74,74,74, max = 255, alpha = 200, names = "grey29")
plot(train.layers[[1]], main='Bioclim 1')
points(bg, col= mycol, pch = 14,cex=.25)
points(acsel, col='black', pch = 16,cex=.35)
plot(st_geometry(map_states), add = TRUE)

train <- rbind(acsel, bg) 
pa_train <- c(rep(1, nrow(acsel)), rep(0, nrow(bg)))
train <- data.frame(cbind(CLASS=pa_train, train))
# make sure no duplicates between background and presences
xv.train.duplicates <- duplicated(data.frame(train$lon,train$lat))
train  <- train[!xv.train.duplicates,]
head(train)

# create spatial points object with 100% of presences
crs <- crs(train.layers[[1]])
train <- train[sample(nrow(train)),]
class.pa <- data.frame(train[,1])
colnames(class.pa) <- c('CLASS')
dataMap.xv  <- SpatialPointsDataFrame(train[,c(2,3)], class.pa,
                                      proj4string =crs)
# write as shp for final modeling (100% of data)
#dir.create(path = "data/geofile-xv-presences/")
writeOGR(dataMap.xv, 'data_2.0/geofile-xv-presences/xv.shp','xv',
         driver='ESRI Shapefile', overwrite = TRUE)

## ------------------------------------------------------------------------

# split into train and test dataframes with extracted values
xv_final <- data.frame(train$lon[train$CLASS==1],train$lat[train$CLASS==1])
colnames(xv_final) <- c('lon', 'lat')
bg_final <- data.frame(train$lon[train$CLASS==0],train$lat[train$CLASS==0])
colnames(bg_final) <- c('lon', 'lat')
train_layers <- list.files(path = 'data_2.0/train-rasters-2.5m/',
                           pattern='asc', full.names=TRUE)
predictors <- stack(train_layers) #, layers = feature_names)

vars <- c('bclim11', 'bclim12', 'bclim14', 'bclim15', 'bclim18', 'bclim2',
          'bclim3', 'bclim4', 'bclim6', 'bclim7', 'bclim8', 'bclim9')
predictors_xv <- predictors[[vars]]
names(predictors_xv) 

group <- kfold(xv_final, 5)
pres_train <- xv_final[group != 1, ]
pres_test <- xv_final[group == 1, ]
write.csv(pres_train, 'data_2.0/xv_train_coords.csv')
group <- kfold(bg_final, 5)
backg_train <- bg_final[group != 1, ] 
backg_test <- bg_final[group == 1, ] 
envtrain_init <- rbind(pres_train, backg_train) 
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train))) 
envtrain <- extract(predictors_xv, envtrain_init)
envtrain <- data.frame(cbind(pa=pb_train, envtrain)) 
envtrain <- envtrain[sample(nrow(envtrain)),]
head(envtrain)
envtest_init <- rbind(pres_test, backg_test)
pb_test <- c(rep(1, nrow(pres_test)), rep(0, nrow(backg_test))) 
envtest <- extract(predictors_xv, envtest_init)
envtest <- data.frame(cbind(pa=pb_test, envtest))
envtest <- envtest[sample(nrow(envtest)),]
head(envtest)
length(envtrain$pa)
length(envtest$pa)

write.csv(envtrain, file ='data_2.0/envtrain_xv.csv', row.names=FALSE)
write.csv(envtest, file ='data_2.0/envtest_xv.csv', row.names=FALSE)

## ------------------------------------------------------------------------

# plot pycaret/sci-kit learn modeling results with false negatives

py_blend <- raster('
    outputsexp_id=101, xantusia_after/blender-baseline/responses.tif')
py_blend <- py_blend*train.layers[[1]]/train.layers[[1]] 
py_blend  <- calc(py_blend , fun=function(x){ x[x < 0.5] <- NA; return(x)})

wrong <- extract(py_blend, data.frame(xv_final$lon, xv_final$lat))
wrong <- xv_final[which(is.na(wrong)),]
wrong

leav <- extract(py_blend, data.frame(xv_raw_Leavitt$lon, xv_raw_Leavitt$lat))
wrong_leav <- xv_raw_Leavitt[which(is.na(leav)),]
wrong_leav

py_blend <- raster('
   outputsexp_id=101, xantusia_after/blender-baseline/probability_1.0.tif')
py_blend <- py_blend*train.layers[[1]]/train.layers[[1]] 

#png(filename = 'outputsexp_id=101, xantusia_after/cont_baseline_xv.png',
#=5, width=2800, height=2000, res=800)
plot(py_blend, main='')#, col = 'forest green')
points(xv_final, col=mycol, pch = 16, cex=.2)
#points(bg_final, col=mycol, pch = 16, cex=.6)
plot(st_geometry(map_states), add = TRUE)
#dev.off()

# compute collinearity for recursive ranker function

mat_1 <-round(cor(extract(predictors, xv_final), method = 'spearman'),7)
mat_1 <-abs(mat_1)
write.csv(data.frame(mat_1), 'data_2.0/collinearity/xv-baseline.csv',
          row.names=FALSE)
write.csv(extract(predictors, xv_final), 'data_2.0/collinearity/raw_data.csv',
          row.names=FALSE)

## ------------------------------------------------------------------------
# yucca brevifolia
yb_raw <- read.csv('data_2.0/raw/yb_GBIF_presences.csv', header = TRUE)
yb_raw <- yb_raw[which(yb_raw$countryCode=='US'),] 
yb <- data.frame(matrix(ncol = 2, nrow = length(yb_raw$decimalLongitude)))
yb[,1] <- yb_raw$decimalLongitude
yb[,2] <- yb_raw$decimalLatitude
yb <- unique(yb)
yb <- yb[complete.cases(yb),]
colnames(yb) <- c('lon','lat')
length(yb$lon)
yb <- yb[which(yb$lon>=ext[1] & yb$lon<=ext[2]),] 
yb <- yb[which(yb$lat>=ext[3] & yb$lat<=ext[4]),] 
length(yb$lon)

# gridded sampling
r_yb <- raster(ext, res = res(train.layers), crs = crs(train.layers))
acsel_yb <- gridSample(data.frame(yb$lon, yb$lat), r_yb, n=2)
head(acsel_yb)
colnames(acsel_yb) <- c('lon', 'lat')
length(acsel_yb$lon) # 1325 after down from 3883 before

# sample background points from a slightly wider extent
bg_yb <- randomPoints(train.layers[[1]], length(acsel_yb$lon),
                      ext=ext, extf = 1) 
colnames(bg_yb) <- c('lon','lat')
bg_yb <- data.frame(bg_yb)
bg_yb <- gridSample(data.frame(bg_yb$lon, bg_yb$lat), r_yb, n=2)
head(bg_yb)
length(bg_yb$bg_yb.lon) # 1325
colnames(bg_yb) <- c('lon', 'lat')

train_yb <- rbind(acsel_yb, bg_yb)
pa_train_yb <- c(rep(1, nrow(acsel_yb)), rep(0, nrow(bg_yb)))
train_yb <- data.frame(cbind(CLASS=pa_train_yb, train_yb))
# make sure no duplicates between background and presences
yb.train.duplicates <- duplicated(data.frame(train_yb$lon,train_yb$lat))
train_yb  <- train_yb[!yb.train.duplicates,]
head(train_yb)

# create spatial points object with 100% of presences
crs <- crs(train.layers[[1]])
train_yb <- train_yb[sample(nrow(train_yb)),]
class.pa <- data.frame(train_yb[,1])
colnames(class.pa) <- c('CLASS')
dataMap.yb  <- SpatialPointsDataFrame(train_yb[,c(2,3)], class.pa,
                                      proj4string = crs)
# write as shp for final modeling (100% of data)
#dir.create(path = "data/geofile-yb-presences/")
writeOGR(dataMap.yb, 'data_2.0/geofile-yb-presences/yb.shp','yb',
         driver='ESRI Shapefile', overwrite = TRUE)

plot(train.layers[[1]], main='Bioclim 1')
points(bg_yb, col= mycol, pch = 16,cex=.25)
points(acsel_yb, col='black', pch = 16,cex=.35)
plot(st_geometry(map_states), add = TRUE)

## ------------------------------------------------------------------------
# split into train and test dataframes with extracted values
yb_final <- data.frame(train_yb$lon[train_yb$CLASS==1],
                       train_yb$lat[train_yb$CLASS==1])
colnames(yb_final) <- c('lon', 'lat')
bg_final_yb <- data.frame(train_yb$lon[train_yb$CLASS==0],
                          train_yb$lat[train_yb$CLASS==0])
colnames(bg_final_yb) <- c('lon', 'lat')
train_layers <- list.files(path = 'data_2.0/train-rasters-2.5m/',
                           pattern='asc', full.names=TRUE)
predictors <- stack(train_layers) #, layers = feature_names)

vars <- c('bclim10', 'bclim12', 'bclim14', 'bclim15', 'bclim17',
          'bclim18', 'bclim2', 'bclim3', 'bclim4', 'bclim6', 'bclim7',
          'bclim8', 'bclim9')
predictors_yb <- predictors[[vars]]
names(predictors_yb) 

group <- kfold(yb_final, 5) 
pres_train_yb <- yb_final[group != 1, ]
pres_test_yb <- yb_final[group == 1, ] 
write.csv(pres_train, 'data_2.0/yb_train_coords.csv')
group <- kfold(bg_final_yb, 5) 
backg_train_yb <- bg_final_yb[group != 1, ] 
backg_test_yb <- bg_final_yb[group == 1, ] 
envtrain_init_yb <- rbind(pres_train_yb, backg_train_yb) 
pb_train_yb <- c(rep(1, nrow(pres_train_yb)), rep(0, nrow(backg_train_yb))) 
envtrain_yb <- extract(predictors_yb, envtrain_init_yb)
envtrain_yb <- data.frame(cbind(pa=pb_train_yb, envtrain_yb)) 
envtrain_yb <- envtrain_yb[sample(nrow(envtrain_yb)),]
head(envtrain_yb)
envtest_init_yb <- rbind(pres_test_yb, backg_test_yb)
pb_test_yb <- c(rep(1, nrow(pres_test_yb)), rep(0, nrow(backg_test_yb))) 
envtest_yb <- extract(predictors_yb, envtest_init_yb)
envtest_yb <- data.frame(cbind(pa=pb_test_yb, envtest_yb))
envtest_yb <- envtest_yb[sample(nrow(envtest_yb)),]
head(envtest_yb)
length(envtrain_yb$pa)
length(envtest_yb$pa)

write.csv(envtrain, file ='data_2.0/envtrain_yb.csv', row.names=FALSE)
write.csv(envtest, file ='data_2.0/envtest_yb.csv', row.names=FALSE)

## ------------------------------------------------------------------------

# plot pycaret/sci-kit learn modeling results

py_blend <- raster(
  'outputsexp_id=101, yucca_b/blender-baseline/responses.tif')
py_blend <- py_blend*train.layers[[1]]/train.layers[[1]] 
py_blend  <- calc(py_blend , fun=function(x){ x[x < 0.5] <- NA; return(x)})

wrong <- extract(py_blend, data.frame(yb_final$lon, yb_final$lat))
wrong <- yb_final[which(is.na(wrong)),]
wrong

py_blend <- raster('
     outputsexp_id=101, yucca_b/blender-baseline/probability_1.0.tif')
py_blend <- py_blend*train.layers[[1]]/train.layers[[1]] 

#png(filename = 'outputsexp_id=101, yucca_b/cont_baseline_yb.png',
#pointsize=5, width=2800, height=2000, res=800)
plot(py_blend, main='') #col = 'forest green')
points(yb_final, col= mycol, pch = 16,cex=.6)
#points(wrong, col = 'red', pch = 18,cex=1)
plot(st_geometry(map_states), add = TRUE)
#dev.off()

# compute collinearity for recursive ranker function
mat_1 <-round(cor(extract(predictors, yb_final), method = 'spearman'),7)
mat_1 <-abs(mat_1)
write.csv(data.frame(mat_1), 'data_2.0/collinearity/yb-baseline.csv',
          row.names=FALSE)
write.csv(extract(predictors, yb_final), 'data_2.0/collinearity/raw_data_yb.csv',
          row.names=FALSE)

# unable to be re-run following code without additional data,
# hence we comment it out, contact for zipped data:

# ## ------------------------------------------------------------------------
# # Future analyses: xantusia
# 
# scenarios = c('ssp126_2041-2060', 'ssp126_2061-2080', 'ssp126_2081-2100',
#               'ssp245_2041-2060', 'ssp245_2061-2080', 'ssp245_2081-2100',
#               'ssp370_2041-2060', 'ssp370_2061-2080', 'ssp370_2081-2100',
#               'ssp585_2041-2060', 'ssp585_2061-2080', 'ssp585_2081-2100')
# 
# overlapped_stats = data.frame(ssp126_2041_2060=c(0,0,0),
#                               ssp126_2061_2080=c(0,0,0),
#                               ssp126_2081_2100=c(0,0,0),
#                               ssp245_2041_2060=c(0,0,0),
#                               ssp245_2061_2080=c(0,0,0),
#                               ssp245_2081_2100=c(0,0,0),
#                               ssp370_2041_2060=c(0,0,0),
#                               ssp370_2061_2080=c(0,0,0),
#                               ssp370_2081_2100=c(0,0,0),
#                               ssp585_2041_2060=c(0,0,0),
#                               ssp585_2061_2080=c(0,0,0),
#                               ssp585_2081_2100=c(0,0,0),
#                               row.names = c('at least two', ' at least five'
#                                             , 'at least eight'))
# 

py_blend <- raster(
   paste('outputsexp_id=101, xantusia',
         '_after/blender-baseline/responses.tif', sep = ""))
py_blend <- py_blend*train.layers[[1]]/train.layers[[1]] #remove nas
py_blend  <- calc(py_blend , fun=function(x){ x[x < 0.5] <- NA; return(x)})
prf100 <- calc(py_blend, fun=function(x){ x[x==1] <- 100; return(x)})
curr_xv <- cellStats(py_blend, function(i, ...) sum(!is.na(i)))
curr_xv

# 
# l = 1
# for (i in scenarios) {
#   BCC <- raster(paste('data/CMIP6/',
#                     i,'/BCC/responses.tif', sep = ""))
#   BCC <- BCC*train.layers[[1]]/train.layers[[1]]
#   Can <- raster(paste('data/CMIP6/',
#                     i,'/Can/responses.tif', sep = ""))
#   Can <- Can*train.layers[[1]]/train.layers[[1]]
#   CNRMcm6 <- raster(paste('data/CMIP6/',
#                     i,'/CNRMcm6/responses.tif', sep = ""))
#   CNRMcm6 <- CNRMcm6*train.layers[[1]]/train.layers[[1]]
#   CNRMesm2 <- raster(paste('data/CMIP6/',
#                     i,'/CNRMesm2/responses.tif', sep = ""))
#   CNRMesm2 <- CNRMesm2*train.layers[[1]]/train.layers[[1]]
#   IPSL <- raster(paste('data/CMIP6/',
#                     i,'/IPSL/responses.tif', sep = ""))
#   IPSL <- IPSL*train.layers[[1]]/train.layers[[1]]
#   Miroc6 <- raster(paste('data/CMIP6/',
#                     i,'/Miroc6/responses.tif', sep = ""))
#   Miroc6 <- Miroc6*train.layers[[1]]/train.layers[[1]]
#   MirocES2L <- raster(paste('data/CMIP6/',
#                     i,'/MirocES2L/responses.tif', sep = ""))
#   MirocES2L <- MirocES2L*train.layers[[1]]/train.layers[[1]]
#   MRI <- raster(paste('data/CMIP6/',
#                     i,'/MRI/responses.tif', sep = ""))
#   MRI <- MRI*train.layers[[1]]/train.layers[[1]]
#   future_overlap <- BCC+Can+CNRMcm6+CNRMesm2+IPSL+Miroc6+MirocES2L+MRI
#   # calculate ratio of overlap
#   proj_and_current <- (prf100 + future_overlap)
#   two_plus_suitable_pixels <- calc(proj_and_current, fun=function(x){ x[x < 102]<-
#     NA; return(x)} )
#   two_plus_suitable_pixels <- cellStats(two_plus_suitable_pixels, 
#                                       function(i, ...) sum(!is.na(i)))
#   overlap_two_pixels <- two_plus_suitable_pixels / curr_xv
#   overlap_two_pixels # the upper error bar
#   five_plus_suitable_pixels <- calc(proj_and_current, fun=function(x){ x[x < 105] <-
#     NA; return(x)} )
#   writeRaster(five_plus_suitable_pixels, filename = paste(
#     'data/CMIP6/', i, '/five_plus_xv.tif', sep = ""), format="GTiff", overwrite = TRUE)
#   five_plus_suitable_pixels <- cellStats(five_plus_suitable_pixels,
#                                        function(i, ...) sum(!is.na(i)))
#   overlap_five_pixels <- five_plus_suitable_pixels / curr_xv
#   overlap_five_pixels # the quasi median
# 
#   
#   eight_plus_suitable_pixels <- calc(proj_and_current, fun=function(x){ x[x < 108] <-
#     NA; return(x)} )
#   eight_plus_suitable_pixels <- cellStats(eight_plus_suitable_pixels,
#                                         function(i, ...) sum(!is.na(i)))
#   overlap_eight_pixels <- eight_plus_suitable_pixels / curr_xv
#   overlap_eight_pixels # the lower error bar
#   overlapped_stats[,l] <- c(overlap_two_pixels, overlap_five_pixels,
#                            overlap_eight_pixels)
#   l = l+1
# }
# 
# overlapped_stats
# 
# #write.csv(overlapped_stats,"/Volumes/HardDrive/overlapped_stats.csv",
#           #row.names = TRUE)
# #png(filename = 'data/CMIP6/ssp370_2081-2100/future_overlap.png',
#     #pointsize=5, width=2800, height=2000, res=800)
# plot(future_overlap)
# plot(proj_and_current) #, col = 'forest green')
# plot(st_geometry(map_states), add = TRUE)
# #dev.off()
# 
# ## ------------------------------------------------------------------------
# # Future analyses: yucca b
# # unable to be re-run without additional data, contact for zipped data
# 
# overlapped_stats_yb = data.frame(ssp126_2041_2060=c(0,0,0),
#                               ssp126_2061_2080=c(0,0,0),
#                               ssp126_2081_2100=c(0,0,0),
#                               ssp245_2041_2060=c(0,0,0),
#                               ssp245_2061_2080=c(0,0,0),
#                               ssp245_2081_2100=c(0,0,0),
#                               ssp370_2041_2060=c(0,0,0),
#                               ssp370_2061_2080=c(0,0,0),
#                               ssp370_2081_2100=c(0,0,0),
#                               ssp585_2041_2060=c(0,0,0),
#                               ssp585_2061_2080=c(0,0,0),
#                               ssp585_2081_2100=c(0,0,0),
#                               row.names = c('at least two', ' at least five'
#                                             , 'at least eight'))
# 
# py_blend <- raster(
#   paste('/Volumes/HardDrive/xantusia-codebase/outputsexp_id=101, yucca',
#         '_b/blender-baseline/responses.tif', sep = ""))
# py_blend <- py_blend*train.layers[[1]]/train.layers[[1]] #remove nas
# py_blend  <- calc(py_blend , fun=function(x){ x[x < 0.5] <- NA; return(x)})
# prf100 <- calc(py_blend, fun=function(x){ x[x==1] <- 100; return(x)})
# curr_yb <- cellStats(py_blend, function(i, ...) sum(!is.na(i)))
# curr_yb
# 
# l = 1
# for (i in scenarios) {
#   BCC <- raster(paste('/Volumes/HardDrive/xantusia-codebase/data/CMIP6/',
#                       i,'/BCC/yb/responses.tif', sep = ""))
#   BCC <- BCC*train.layers[[1]]/train.layers[[1]]
#   Can <- raster(paste('/Volumes/HardDrive/xantusia-codebase/data/CMIP6/',
#                       i,'/Can/yb/responses.tif', sep = ""))
#   Can <- Can*train.layers[[1]]/train.layers[[1]]
#   CNRMcm6 <- raster(paste('/Volumes/HardDrive/xantusia-codebase/data/CMIP6/',
#                       i,'/CNRMcm6/yb/responses.tif', sep = ""))
#   CNRMcm6 <- CNRMcm6*train.layers[[1]]/train.layers[[1]]
#   CNRMesm2 <- raster(paste('/Volumes/HardDrive/xantusia-codebase/data/CMIP6/',
#                       i,'/CNRMesm2/yb/responses.tif', sep = ""))
#   CNRMesm2 <- CNRMesm2*train.layers[[1]]/train.layers[[1]]
#   IPSL <- raster(paste('/Volumes/HardDrive/xantusia-codebase/data/CMIP6/',
#                       i,'/IPSL/yb/responses.tif', sep = ""))
#   IPSL <- IPSL*train.layers[[1]]/train.layers[[1]]
#   Miroc6 <- raster(paste('/Volumes/HardDrive/xantusia-codebase/data/CMIP6/',
#                       i,'/Miroc6/yb/responses.tif', sep = ""))
#   Miroc6 <- Miroc6*train.layers[[1]]/train.layers[[1]]
#   MirocES2L <- raster(paste('/Volumes/HardDrive/xantusia-codebase/data/CMIP6/',
#                       i,'/MirocES2L/yb/responses.tif', sep = ""))
#   MirocES2L <- MirocES2L*train.layers[[1]]/train.layers[[1]]
#   MRI <- raster(paste('/Volumes/HardDrive/xantusia-codebase/data/CMIP6/',
#                       i,'/MRI/yb/responses.tif', sep = ""))
#   MRI <- MRI*train.layers[[1]]/train.layers[[1]]
#   future_overlap <- BCC+Can+CNRMcm6+CNRMesm2+IPSL+Miroc6+MirocES2L+MRI
#   # calculate ratio of overlap
#   proj_and_current <- (prf100 + future_overlap)
#   two_plus_suitable_pixels <- calc(proj_and_current, fun=function(x){ x[x < 102]<-
#     NA; return(x)} )
#   two_plus_suitable_pixels <- cellStats(two_plus_suitable_pixels, 
#                                         function(i, ...) sum(!is.na(i)))
#   overlap_two_pixels <- two_plus_suitable_pixels / curr_xv
#   overlap_two_pixels # the upper error bar
#   five_plus_suitable_pixels <- calc(proj_and_current, fun=function(x){ x[x < 105] <-
#     NA; return(x)} )
#   writeRaster(five_plus_suitable_pixels, filename = paste(
#     'data/CMIP6/', i, '/five_plus_yb.tif', sep = ""), format="GTiff")
#   five_plus_suitable_pixels <- cellStats(five_plus_suitable_pixels,
#                                          function(i, ...) sum(!is.na(i)))
#   overlap_five_pixels <- five_plus_suitable_pixels / curr_xv
#   overlap_five_pixels # the quasi median
#   eight_plus_suitable_pixels <- calc(proj_and_current, fun=function(x){ x[x < 108] <-
#     NA; return(x)} )
#   eight_plus_suitable_pixels <- cellStats(eight_plus_suitable_pixels,
#                                           function(i, ...) sum(!is.na(i)))
#   overlap_eight_pixels <- eight_plus_suitable_pixels / curr_xv
#   overlap_eight_pixels # the lower error bar
#   overlapped_stats_yb[,l] <- c(overlap_two_pixels, overlap_five_pixels,
#                             overlap_eight_pixels)
#   l = l+1
# }
# 
# overlapped_stats_yb
# 
# # calculate species overlap
# 
# py_blend <- raster(paste('outputsexp_id=101, yucca',
#         '_b/blender-baseline/responses.tif', sep = ""))
# py_blend1 <- raster(paste('outputsexp_id=101, xantusia',
#         '_after/blender-baseline/responses.tif', sep = ""))
# both_species = py_blend1 + py_blend
# both_species  <- calc(both_species  , fun=
#                 function(x){ x[x < 0.5] <- NA; return(x)})
# plot(both_species)
# union <- cellStats(both_species, function(i, ...) sum(!is.na(i)))
# both_species  <- calc(both_species  , fun=
#                         function(x){ x[x < 1.5] <- NA; return(x)})
# intersection <- cellStats(both_species, function(i, ...) sum(!is.na(i)))
# intersection/union
# 
# plot(both_species) #, color = 'greens')
# plot(st_geometry(map_states), add = TRUE)
# 
# overlapped_sp = c(0,0,0,0,0,0,0,0,0,0,0,0)
# 
# l = 1
# for (i in scenarios) {
#   py_blend <- raster(paste('data/CMIP6/', i, '/five_plus_xv.tif', sep = ""))
#   py_blend  <- calc(py_blend  , fun=
#         function(x){ x[is.na(x)] <- 0; return(x)})
#   py_blend1 <- raster(paste('data/CMIP6/', i, '/five_plus_yb.tif', sep = ""))
#   py_blend1  <- calc(py_blend1  , fun=
#                       function(x){ x[is.na(x)] <- 0; return(x)})
#   both_species = py_blend1 + py_blend
#   both_species  <- calc(both_species  , fun=
#                           function(x){ x[x < 105] <- NA; return(x)})
#   union <- cellStats(both_species, function(i, ...) sum(!is.na(i)))
#   both_species  <- calc(both_species  , fun=
#                           function(x){ x[x < 210] <- NA; return(x)})
#   intersection <- cellStats(both_species, function(i, ...) sum(!is.na(i)))
#   overlapped_sp[l] <- intersection/union
#   l <- l+1
# }
# #species overlap:
# mean(as.numeric(overlapped_sp))
