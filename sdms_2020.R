#' ---
#' title: "sdms_2020"
#' author: "Daniel Furman"
#' date: "2020"
#' output: html_document
#' ---
#' 
#' Script constructs rforest, extra trees, and lgbm models for 
#' species distribution spatial prediction (sdms).
#' 
#' For an introduction to sdms and dismo:
#' https://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf
#' 
#' 

## ------------------------------------------------------------------------
# Libraries needed for this file, download with install.packages("<command name>") 

library(dismo)
library(raster)
library(maptools)
library(randomForest)
library(lightgbm)
library(extraTrees)

## ------------------------------------------------------------------------

#' First clone the following data repository : 
#' https://github.com/daniel-furman/xvigilis-data-projections
#' 
#' I have unzipped its contents at '/Volumes/HardDrive/'

path <- '/Volumes/HardDrive/' # change to your own path

set.seed(100) #we decided this seed using random number generator

xant <- read.csv(paste(path, 'xvigilis-data-main/xant_complete.csv', 
                      sep = ""),header = TRUE)
xantusia_unique <- unique(xant) #xantusia without duplicates
xantusia_unique <- xantusia_unique[complete.cases(xantusia_unique),] 
e <- extent(min(xantusia_unique$longitude),max(xantusia_unique$longitude),min(
  xantusia_unique$latitude),max(xantusia_unique$latitude))
length_presences <- length(xantusia_unique[,1])

## ------------------------------------------------------------------------

train_layers <- list.files(path=paste(path,'/xvigilis-data-main/train_tifs',
                                      sep = ""), pattern='asc', full.names=TRUE)
train_layers

predictors <- stack(train_layers)
predictors <- predictors[[c(-1,-2,-3,-5,-8,-9,-12,-14,-15,-20,-22)]]
predictors

mask <- raster(train_layers[1])

crs(mask) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

bg2 <- randomPoints(mask, length_presences*3, ext=e, extf = 4) 

colnames(bg2) <- c('lon','lat')

# Extracting values from rasters
#.rs.unloadPackage('tidyr')
presvals <- extract(predictors, xantusia_unique)
absvals <- extract(predictors, bg2)

group <- kfold(xantusia_unique, 5) #divide whole dataset into 5 groups
pres_train <- xantusia_unique[group != 1, ] #make a train set of 80 percent
pres_test <- xantusia_unique[group == 1, ] #make a test set of 20 percent

## ------------------------------------------------------------------------

group <- kfold(bg2, 5) #make background set of 
backg_train <- bg2[group != 1, ] #make back train
backg_test <- bg2[group == 1, ] #make back test

colnames(pres_train) <- c('lon','lat') #make colnames the same
colnames(backg_train) <- c('lon','lat') #make colnames the same

train <- rbind(pres_train, backg_train) #training back and pres binded

pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train))) #col of ones and zeros

envtrain <- extract(predictors, train) #extract raster values for training
envtrain <- data.frame(cbind(pa=pb_train, envtrain) ) #bind 

testpres <- data.frame(extract(predictors, pres_test) ) #make test pres
testbackg <- data.frame(extract(predictors, backg_test) )

testpres <- testpres[complete.cases(testpres),] #remove any nas
testbackg <- testbackg[complete.cases(testbackg),]
envtrain <- envtrain[complete.cases(envtrain),] 

test <- rbind(testpres, testbackg)
pb_test <- c(rep(1, nrow(testpres)), rep(0, nrow(testbackg))) #col of ones and zeros

## ------------------------------------------------------------------------

head(envtrain)
head(testbackg)
head(testpres)

## ------------------------------------------------------------------------

envtrain_corr <- envtrain[,c(-2,-3,-4,-6,-9,-10,-13,-15,-16,-21,-23)] #remove correlates
testpres_corr <- testpres[,c(-1,-2,-3,-5,-8,-9,-12,-14,-15,-20,-22)] #remove correlates
testbackg_corr <- testbackg[,c(-1,-2,-3,-5,-8,-9,-12,-14,-15,-20,-22)] #remove correlates
testbackg_corr
names(envtrain_corr) ##  [1] "pa"        "bclim12"   "bclim14"   "bclim15"   "bclim18"  
##  [6] "bclim19"   "bclim3"    "bclim6"    "bclim7"    "bclim8"   
## [11] "bclim9

write.csv(envtrain_corr, file ='envtrain_corr.csv')
write.csv(testpres_corr, file ='testpres_corr.csv')
write.csv(testbackg_corr, file ='testbackg_corr.csv')



#' 
## ------------------------------------------------------------------------

data(wrld_simpl)
plot(predictors, 1)
#plot(wrld_simpl, add=TRUE)
points(xantusia_unique, col='blue', pch = 16,cex=.4)

plot(predictors, 1)
#plot(wrld_simpl, add=TRUE)
points(bg2, col='red', pch = 1,cex=.4)

#' 
#' 
## ------------------------------------------------------------------------

# fit Rforest

model_corr <- factor(pa) ~
  bclim3   + bclim6 + bclim7 + bclim8 + bclim9   + bclim12  + bclim14 + bclim15 + bclim18 + bclim19
rf_binary <- tuneRF(envtrain_corr[,2:11], factor(envtrain_corr$pa), ntreeTry = 100, doBest = TRUE)

val.pred <- predict(rf_binary, rbind(testpres_corr,testbackg_corr))
tab <- table(observed =  pb_test, predicted = val.pred)

tab # validation confusion matrix

print('class 0 error rforest:')
print(tab[1,2]/(tab[1,1]+tab[1,2])) 
print('class 1 error rforest:') 
print(tab[2,1]/(tab[2,1]+tab[2,2]))

# fit LGBM

dtrain <- lgb.Dataset(
   data = as.matrix(envtrain_corr[,2:11])
   , label = as.numeric(envtrain_corr$pa))
dtrain <- lgb.Dataset.construct(dtrain)
bst <- lightgbm(
   data =  dtrain
   , objective = "binary")
pred <- predict(bst, as.matrix(rbind(testpres_corr,testbackg_corr)))
err <- mean(as.numeric(pred > 0.5) != pb_test)
print(paste("test-error lgbm=", err))
# pbst1<- predict(predictors, bst, ext=ext) #current #doesn't work with bst

# fit etrees

etrees_binary <- extraTrees(envtrain_corr[,2:11], factor(envtrain_corr$pa))

val.pred <- predict(etrees_binary, rbind(testpres_corr,testbackg_corr))
tab <- table(observed =  pb_test, predicted = val.pred)

tab # validation confusion matrix

print('class 0 error etrees:')
print(tab[1,2]/(tab[1,1]+tab[1,2])) 
print('class 1 error etrees:') 
print(tab[2,1]/(tab[2,1]+tab[2,2]))

## ------------------------------------------------------------------------

ext <- extent(-125.0208, -92.00083, 20, 46.9975) #set extent

prf1<- predict(predictors, rf_binary, ext=ext) #current random forest
petrees1<- predict(predictors, etrees_binary, ext=ext) #current extra trees

#png(filename = '/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/rf_current_setseed100.png',
    #pointsize=5, width=2800, height=2000, res=800)
plot(prf1, main='random forest, classification ')
points(bg2, col='red', pch = 16,cex=.2)
points(xantusia_unique, col='black', pch = 16,cex=.2)
plot(wrld_simpl, add=TRUE, border='dark grey')
#dev.off()

plot(petrees1, main='extra trees, classification ')
points(bg2, col='red', pch = 16,cex=.2)
points(xantusia_unique, col='black', pch = 16,cex=.2)
plot(wrld_simpl, add=TRUE, border='dark grey')

crs(prf1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(petrees1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

petrees1 <- calc(petrees1, fun=function(x){ x[x < 1] <- NA; return(x)} )
prf1 <- calc(prf1, fun=function(x){ x[x < 1] <- NA; return(x)} )


prediction <- petrees1 + prf1
prediction <- calc(prediction, fun=function(x){ x[x < 2] <- NA; return(x)} )
prediction <- calc(prediction, fun=function(x){ x[x == 2] <- 1; return(x)} )

#png(filename = '/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/etrees_rforest.png',
#pointsize=5, width=2800, height=2000, res=800)
plot(prediction, main='etrees and rforest', col = 'blue')
points(bg2, col='red', pch = 16,cex=.2)
points(xantusia_unique, col='black', pch = 16,cex=.2)
plot(wrld_simpl, add=TRUE, border='dark grey')
#dev.off()

crs(prediction) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"


print(cellStats(prediction, function(i, ...) sum(!is.na(i)))) 
print(cellStats(petrees1, function(i, ...) sum(!is.na(i)))) 
print(cellStats(prf1, function(i, ...) sum(!is.na(i)))) 

curr <- cellStats(prediction, function(i, ...) sum(!is.na(i)))