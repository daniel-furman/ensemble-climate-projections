#' ---
#' title: "sdms_2020"
#' author: "Daniel Furman"
#' date: "2020"
#' output: html_document
#' ---
#' 
#' For an introduction to sdms and dismo:
#' https://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf
#' 

## ------------------------------------------------------------------------
# Libraries needed for this file, download with install.packages("<command name>") 

library(dismo)
library(raster)
library(maptools)
library(randomForest)

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

colnames(pres_train) = c('lon','lat') #make colnames the same
colnames(backg_train) = c('lon','lat') #make colnames the same

train <- rbind(pres_train, backg_train) #training back and pres binded

pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train))) #col of ones and zeros

envtrain <- extract(predictors, train) #extract raster values for training

envtrain <- data.frame(cbind(pa=pb_train, envtrain) ) #bind 

testpres <- data.frame(extract(predictors, pres_test) ) #make test pres
testbackg <- data.frame(extract(predictors, backg_test) )

testpres = testpres[complete.cases(testpres),] #remove any nas
testbackg = testbackg[complete.cases(testbackg),]
envtrain = envtrain[complete.cases(envtrain),] 

## ------------------------------------------------------------------------

head(envtrain)
head(testbackg)
head(testpres)

## ------------------------------------------------------------------------

envtrain_corr <- envtrain[,c(-2,-3,-4,-6,-9,-10,-13,-15,-16,-21,-23)] #remove correlates
testpres_corr <- testpres[,c(-1,-2,-3,-5,-8,-9,-12,-14,-15,-20,-22)] #remove correlates
testbackg_corr <- testbackg[,c(-1,-2,-3,-5,-8,-9,-12,-14,-15,-20,-22)] #remove correlates
names(envtrain_corr) ##  [1] "pa"        "bclim12"   "bclim14"   "bclim15"   "bclim18"  
##  [6] "bclim19"   "bclim3"    "bclim6"    "bclim7"    "bclim8"   
## [11] "bclim9

#write.csv(envtrain_corr, file ='envtrain_corr.csv')
#write.csv(testpres_corr, file ='testpres_corr.csv')
#rite.csv(testbackg_corr, file ='testbackg_corr.csv')



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

# fit models Rf

model_corr <- factor(pa) ~
  bclim3   + bclim6 + bclim7 + bclim8 + bclim9   + bclim12  + bclim14 + bclim15 + bclim18 + bclim19

rf_binary <- tuneRF(envtrain_corr[,2:11], factor(envtrain_corr$pa), ntreeTry = 500, doBest = TRUE)
rf_binary
varImpPlot(rf_binary)

val.pred <- predict(rf_binary, rbind(test_pres_k,test_backg_k))
tab <- table(observed =  rbind(test_pres_k,test_backg_k)$pa, predicted = val.pred)

tab # validation confusion matrix

print('class 0 validation error:')
print(tab[1,2]/(tab[1,1]+tab[1,2])) 
print('class 1 validation error:') 
print(tab[2,1]/(tab[2,1]+tab[2,2]))

## ------------------------------------------------------------------------

ext <- extent(-125.0208, -92.00083, 20, 46.9975) # set extent

prf1<- predict(predictors, rf_binary, ext=ext) #current 

#png(filename = '/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/rf_current_setseed100.png',
    #pointsize=5, width=2800, height=2000, res=800)
plot(prf1, main='random forest, classification ')
points(bg2, col='red', pch = 16,cex=.2)
points(xantusia_unique, col='black', pch = 16,cex=.2)
plot(wrld_simpl, add=TRUE, border='dark grey')
#dev.off()


crs(prf1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

all <- calc(prf1, fun=function(x){ x[x < 1] <- NA; return(x)} )

curr <- cellStats(all, function(i, ...) sum(!is.na(i)))
print(curr)

#source('sdms_2020.R')


