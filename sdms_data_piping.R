#' ---
#' title: "sdms_data_piping"
#' author: "Daniel Furman"
#' date: "2020"
#' #' ---
#' 
#' This R script pipes presence and absence data for xantusia
#' vigilis and yucca brevifolia from the raw Worldclim and 
#' GBIF data.
#' 
#' It is part of a new codebase, not meant to be run with 
#' the other "ensemble-climate-projections" repo directly
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
library(maptools)

set.seed(100) # seed using random number generator

## ------------------------------------------------------------------------

# Raster data, for the current distribution 

# set study space extent:
ext <- extent(-122.767 , -109.3192 , 31.91825  , 38.1922)
# download raw bioclimatic tifs first from WorldClim 2.1
# grab those at the 2.5 arcminutes 
train.files <- list.files(path='data/raw/wc2.1_2.5m_bio/',
                           pattern='tif', full.names=TRUE)
train.layers <- stack(train.files)
train.layers <- crop(train.layers, ext)
bcvar.index <- c(1,10:19, 2:9)
for (i in c(1:19)){
  names(train.layers[[i]])<-paste('bclim', bcvar.index[i], sep = "")
}
names(train.layers)
print(train.files)
dir.create(path = "data/train-rasters-2.5m/")
writeRaster(train.layers, filename=paste('data/train-rasters-2.5m/',
  names(train.layers), sep = ""), bylayer=TRUE,format="ascii")
# also load a finer scale, at 30 arcseconds 
train.files <- list.files(path='data/raw/wc2.1_30s_bio/',
                           pattern='tif', full.names=TRUE)
train.fine <- stack(train.files)
train.fine <- crop(train.fine, ext)
for (i in c(1:19)){
  names(train.fine[[i]])<-paste('bclim', bcvar.index[i], sep = "")
}
dir.create(path = "data/train-rasters-30s/")
writeRaster(train.fine, filename=paste('data/train-rasters-30s/',
  names(train.fine), sep = ""), bylayer=TRUE,format="ascii")

## ------------------------------------------------------------------------

# Xantusia vigilis:
# download GBIF database:
xv_raw <- read.csv( 'data/raw/xv_GBIF_presences.csv', header = TRUE)
xv <- data.frame(matrix(ncol = 2, nrow = length(xv_raw$decimalLongitude)))
xv[,1] <- xv_raw$decimalLongitude
xv[,2] <- xv_raw$decimalLatitude
colnames(xv) <- c('lon','lat')
xv.duplicates <- duplicated(xv)
for (i in c(1:length(xv_raw$gbifID))){
  xv$Id[i] <- 'GBIF - Id ' + as.String(xv_raw$gbifID[i])
}
xv <- xv[!xv.duplicates,]
# next load the Leavit et al., 2007 database
xv_raw_Leavitt <- read.csv( 'data/raw/Leavitt_2007.csv', header = FALSE)
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
bg <- gridSample(data.frame(bg$lon, bg$lat), r, n=2) # grid sampled
head(bg)
length(bg$bg.lon) #907
colnames(bg) <- c('lon', 'lat')

train <- rbind(acsel, bg)  # combine with presences 
pa_train <- c(rep(1, nrow(acsel)), rep(0, nrow(bg)))
train <- data.frame(cbind(CLASS=pa_train, train)) # final dataframe
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
dir.create(path = "data/geofile-xv-presences/")
writeOGR(dataMap.xv, 'data/geofile-xv-presences/xv.shp','xv',
         driver='ESRI Shapefile')
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

## ------------------------------------------------------------------------

# split into train and test dataframes with extracted values
xv_final <- data.frame(train$lon[train$CLASS==1],train$lat[train$CLASS==1])
colnames(xv_final) <- c('lon', 'lat')
bg_final <- data.frame(train$lon[train$CLASS==0],train$lat[train$CLASS==0])
colnames(bg_final) <- c('lon', 'lat')
train_layers <- list.files(path = 
  'data/train-rasters-2.5m/', pattern='asc', full.names=TRUE)
predictors <- stack(train_layers) #, layers = feature_names)

#vars <- c('bclim11', 'bclim12', 'bclim14', 'bclim15', 'bclim18', 'bclim2',
          #'bclim3', 'bclim4', 'bclim6', 'bclim7', 'bclim8', 'bclim9')
#predictors <- predictors[[vars]]
names(predictors) # confirm the set of variables you want

group <- kfold(xv_final, 5) # divide into 5 groups
pres_train <- xv_final[group != 1, ] # make a train set of 80 percent
pres_test <- xv_final[group == 1, ] # make a test set of 20 percent
write.csv(pres_train, 'data/xv_train_coords.csv')
group <- kfold(bg_final, 5) # make background set of 
backg_train <- bg_final[group != 1, ] # make back train
backg_test <- bg_final[group == 1, ] # make back test
envtrain_init <- rbind(pres_train, backg_train) 
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train))) 
envtrain <- extract(predictors, envtrain_init)
envtrain <- data.frame(cbind(pa=pb_train, envtrain)) 
envtrain <- envtrain[sample(nrow(envtrain)),]
head(envtrain)
envtest_init <- rbind(pres_test, backg_test)
pb_test <- c(rep(1, nrow(pres_test)), rep(0, nrow(backg_test))) 
envtest <- extract(predictors, envtest_init)
envtest <- data.frame(cbind(pa=pb_test, envtest))
envtest <- envtest[sample(nrow(envtest)),]
head(envtest)
length(envtrain$pa)
length(envtest$pa)

write.csv(envtrain, file ='data/envtrain_xv.csv', row.names=FALSE)
write.csv(envtest, file ='data/envtest_xv.csv', row.names=FALSE)

## ------------------------------------------------------------------------

# plot pycaret/sci-kit learn modeling results with false negatives

py_blend <- raster('
        outputsexp_id=101, xantusia_after/blender-baseline/responses.tif')
py_blend <- py_blend*train.layers[[1]]/train.layers[[1]] #remove nas
py_blend  <- calc(py_blend , fun=function(x){ x[x < 0.5] <- NA; return(x)} )

wrong <- extract(py_blend, data.frame(xv_final$lon, xv_final$lat))
wrong <- xv_final[which(is.na(wrong)),]
wrong

leav <- extract(py_blend, data.frame(xv_raw_Leavitt$lon, xv_raw_Leavitt$lat))
wrong_leav <- xv_raw_Leavitt[which(is.na(leav)),]
wrong_leav

py_blend <- raster('
    outputsexp_id=101, xantusia_after/blender-baseline/probability_1.0.tif')
py_blend <- py_blend*train.layers[[1]]/train.layers[[1]] #remove nas

#png(filename = 'outputsexp_id=101, xantusia_after/cont_baseline_xv.png',
    #pointsize=5, width=2800, height=2000, res=800)
plot(py_blend, main='')#, col = 'forest green')
points(xv_final, col=mycol, pch = 16, cex=.3)
#points(wrong_leav, col = 'red', pch = 16,cex=1)
#points(wrong, col = 'red', pch = 18,cex=1)
plot(st_geometry(map_states), add = TRUE)
#dev.off()

# compute collinearity for recursive ranker function

#mat_1 <-round(cor(extract(predictors, xv_final), method = 'spearman'),7)
#mat_1 <-abs(mat_1)
#write.csv(data.frame(mat_1), 'data/collinearity/xv-baseline.csv',
          #row.names=FALSE)
#write.csv(extract(predictors, xv_final), 'data/collinearity/raw_data.csv',
          #row.names=FALSE)

## ------------------------------------------------------------------------
# yucca brevifolia
yb_raw <- read.csv('data/raw/yb_GBIF_presences.csv', header = TRUE)
yb_raw <- yb_raw[which(yb_raw$countryCode=='US'),] # restrict to US
yb <- data.frame(matrix(ncol = 2, nrow = length(yb_raw$decimalLongitude)))
yb[,1] <- yb_raw$decimalLongitude
yb[,2] <- yb_raw$decimalLatitude
yb <- unique(yb) # xantusia without duplicates
yb <- yb[complete.cases(yb),]
colnames(yb) <- c('lon','lat')
length(yb$lon)
yb <- yb[which(yb$lon>=ext[1] & yb$lon<=ext[2]),] 
yb <- yb[which(yb$lat>=ext[3] & yb$lat<=ext[4]),] 
length(yb$lon)

# gridded sampling
r <- raster(ext, res = res(train.layers), crs = crs(train.layers))
acsel <- gridSample(data.frame(yb$lon, yb$lat), r, n=2)
head(acsel)
colnames(acsel) <- c('lon', 'lat')
length(acsel$lon) # 1325 after down from 3883 before

# sample background points from a slightly wider extent
bg <- randomPoints(train.layers[[1]], length(acsel$lon),
                   ext=ext, extf = 1) 
colnames(bg) <- c('lon','lat')
bg <- data.frame(bg)
bg <- gridSample(data.frame(bg$lon, bg$lat), r, n=2) # grid sampled
head(bg)
length(bg$bg.lon) # 1325
colnames(bg) <- c('lon', 'lat')

train <- rbind(acsel, bg)  # combine with presences 
pa_train <- c(rep(1, nrow(acsel)), rep(0, nrow(bg)))
train <- data.frame(cbind(CLASS=pa_train, train)) # final dataframe
# make sure no duplicates between background and presences
yb.train.duplicates <- duplicated(data.frame(train$lon,train$lat))
train  <- train[!yb.train.duplicates,]
head(train)

# create spatial points object with 100% of presences
crs <- crs(train.layers[[1]])
train <- train[sample(nrow(train)),]
class.pa <- data.frame(train[,1])
colnames(class.pa) <- c('CLASS')
dataMap.yb  <- SpatialPointsDataFrame(train[,c(2,3)], class.pa,
                                      proj4string = crs)
# write as shp for final modeling (100% of data)
dir.create(path = "data/geofile-yb-presences/")
writeOGR(dataMap.yb, 'data/geofile-yb-presences/yb.shp','yb',
         driver='ESRI Shapefile')
  
plot(train.layers[[1]], main='Bioclim 1')
points(bg, col= mycol, pch = 16,cex=.25)
points(acsel, col='black', pch = 16,cex=.35)
plot(st_geometry(map_states), add = TRUE)

## ------------------------------------------------------------------------
# split into train and test dataframes with extracted values
yb_final <- data.frame(train$lon[train$CLASS==1],train$lat[train$CLASS==1])
colnames(yb_final) <- c('lon', 'lat')
bg_final <- data.frame(train$lon[train$CLASS==0],train$lat[train$CLASS==0])
colnames(bg_final) <- c('lon', 'lat')
train_layers <- list.files(path = 
  'data/train-rasters-2.5m/', pattern='asc', full.names=TRUE)
predictors <- stack(train_layers) #, layers = feature_names)

vars <- c('bclim10', 'bclim12', 'bclim14', 'bclim15', 'bclim17',
          'bclim18', 'bclim2', 'bclim3', 'bclim4', 'bclim6', 'bclim7',
          'bclim8', 'bclim9')
predictors <- predictors[[vars]]
names(predictors) # confirm the set of variables you want

group <- kfold(yb_final, 5) # divide into 5 groups
pres_train <- yb_final[group != 1, ] # make a train set of 80 percent
pres_test <- yb_final[group == 1, ] # make a test set of 20 percent
write.csv(pres_train, 'data/yb_train_coords.csv')
group <- kfold(yb_final, 5) # make background set of 
backg_train <- bg_final[group != 1, ] # make back train
backg_test <- bg_final[group == 1, ] # make back test
envtrain_init <- rbind(pres_train, backg_train) 
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train))) 
envtrain <- extract(predictors, envtrain_init)
envtrain <- data.frame(cbind(pa=pb_train, envtrain)) 
envtrain <- envtrain[sample(nrow(envtrain)),]
head(envtrain)
envtest_init <- rbind(pres_test, backg_test)
pb_test <- c(rep(1, nrow(pres_test)), rep(0, nrow(backg_test))) 
envtest <- extract(predictors, envtest_init)
envtest <- data.frame(cbind(pa=pb_test, envtest))
envtest <- envtest[sample(nrow(envtest)),]
head(envtest)
length(envtrain$pa)
length(envtest$pa)

write.csv(envtrain, file ='data/envtrain_yb.csv', row.names=FALSE)
write.csv(envtest, file ='data/envtest_yb.csv', row.names=FALSE)

## ------------------------------------------------------------------------

# plot pycaret/sci-kit learn modeling results

py_blend <- raster(
  'outputsexp_id=101, yucca_b/blender-baseline/responses.tif')
py_blend <- py_blend*train.layers[[1]]/train.layers[[1]] #remove nas
py_blend  <- calc(py_blend , fun=function(x){ x[x < 0.5] <- NA; return(x)} )

wrong <- extract(py_blend, data.frame(yb_final$lon, yb_final$lat))
wrong <- yb_final[which(is.na(wrong)),]
wrong

py_blend <- raster('
    outputsexp_id=101, yucca_b/blender-baseline/probability_1.0.tif')
py_blend <- py_blend*train.layers[[1]]/train.layers[[1]] #remove nas

png(filename = 'outputsexp_id=101, yucca_b/cont_baseline_yb.png',
pointsize=5, width=2800, height=2000, res=800)
plot(py_blend, main='') #col = 'forest green')
points(yb_final, col= mycol, pch = 16,cex=.3)
#points(wrong, col = 'red', pch = 18,cex=1)
plot(st_geometry(map_states), add = TRUE)
dev.off()

# compute collinearity for recursive ranker function

#mat_1 <-round(cor(extract(predictors, yb_final), method = 'spearman'),7)
#mat_1 <-abs(mat_1)
#write.csv(data.frame(mat_1), 'data/collinearity/yb-baseline.csv',
          #row.names=FALSE)
#write.csv(extract(predictors, yb_final), 'data/collinearity/raw_data_yb.csv',
         #row.names=FALSE)


