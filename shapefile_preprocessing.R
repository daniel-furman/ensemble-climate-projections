# ---
# title: "Pre-processing for projections"
# author: "Daniel Furman"
# date: "6/27/2020"
# output: html_document
# ---
  
# This notebook processes CMIP6 TIF files and saves them to file. 
# They are then uploaded to GitHub : 
# https://github.com/daniel-furman/xantusia-data

### Load Climate Files 

# Load packages
library(dismo)
library(rgdal)
library(DescTools)

tag <- 'ssp585_2061-2080' #change to your own file path

files <- list.files(path=paste('/Volumes/HardDrive/',tag, sep = ""), pattern='tif',
full.names=TRUE ) 
print('these are the raw tif files:')

files

# dictionary of model names
#  {'BCC-CSM2-MR', 'CanESM5','CNRM-CM6-1','CNRM-ESM2-1’, 'IPSL-CM6A-LR'
# ,'MIROC-ES2L','MIROC6', and 'MRI-ESM2-0'}


### Stacking Climate Rasters 

#  Now that everything is loaded into R, we  create rasters stacks and save them as
# new tif files. These rasters are now all less than 1.5 mb, and can be uploaded
# to GitHub as a result for collaboration and public storage. 

# Initialize empty model vector
BCC_stack_layers <- stack()
Can_stack_layers <- stack()
CNRMcm6_stack_layers <- stack()
CNRMesm2_stack_layers <- stack()
IPSL_stack_layers <- stack()
MirocES2L_stack_layers <- stack()
Miroc6_stack_layers <- stack()
MRI_stack_layers <- stack()

ext <- extent(-125.0208, -92.00083, 20, 46.9975) # set extent

model_stacks <- c(BCC_stack_layers,Can_stack_layers,CNRMcm6_stack_layers,
                  CNRMesm2_stack_layers,IPSL_stack_layers,MirocES2L_stack_layers,
                  Miroc6_stack_layers ,MRI_stack_layers) #create list of models

# populate raster stacks for each GCM model and their 19 bands
for (i in c(1:8)){ 
  for (l in c(1:19)){
    bclim <- raster(files[[i]], band = l) #for all 19 bands make a raster
    bclim <- crop(bclim, ext) #crop the raster to its extent
    model_stacks[[i]] <- stack(model_stacks[[i]], bclim) #stack the bclim variables 
  }
}

# Rename stacks and layers:

BCC <- model_stacks[[1]]
Can <- model_stacks[[2]]
CNRMcm6 <- model_stacks[[3]] 
CNRMesm2 <- model_stacks[[4]] 
IPSL <- model_stacks[[5]] 
MirocES2L <- model_stacks[[6]]
Miroc6 <- model_stacks[[7]]
MRI <- model_stacks[[8]] 

# Rename the variables as bclimi

for (i in c(1:19)) {
  names(BCC[[i]])<-paste('bclim', i, sep = "")
  names(Can[[i]])<-paste('bclim', i, sep = "")
  names(CNRMcm6[[i]])<-paste('bclim', i, sep = "")
  names(CNRMesm2[[i]])<-paste('bclim', i, sep = "")
  names(IPSL[[i]])<-paste('bclim', i, sep = "")
  names(MirocES2L[[i]])<-paste('bclim', i, sep = "")
  names(Miroc6[[i]])<-paste('bclim', i, sep = "")
  names(MRI[[i]])<-paste('bclim', i, sep = "")
}

BCC <- dropLayer(BCC,c(1,2,4,5,10,11,13,16,17)) #drop correlated variables
Can <- dropLayer(Can,c(1,2,4,5,10,11,13,16,17))
CNRMcm6 <- dropLayer(CNRMcm6,c(1,2,4,5,10,11,13,16,17))
CNRMesm2 <- dropLayer(CNRMesm2,c(1,2,4,5,10,11,13,16,17)) 
IPSL <- dropLayer(IPSL,c(1,2,4,5,10,11,13,16,17)) 
MirocES2L <- dropLayer(MirocES2L,c(1,2,4,5,10,11,13,16,17)) 
Miroc6 <- dropLayer(Miroc6,c(1,2,4,5,10,11,13,16,17)) 
MRI <- dropLayer(MRI,c(1,2,4,5,10,11,13,16,17)) 

# Write Rasters to file and share to github repo


writeRaster(BCC, filename=paste('/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/stacks/',
                                tag, '/BCC/', names(BCC), sep = ""), bylayer=TRUE,format="GTiff")
writeRaster(Can, filename=paste('/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/stacks/',
                                tag, '/Can/', names(Can), sep = ""), bylayer=TRUE,format="GTiff")
writeRaster(CNRMcm6, filename=paste('/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/stacks/',
                                    tag, '/CNRMcm6/', names(CNRMcm6), sep = ""), bylayer=TRUE,format="GTiff")
writeRaster(CNRMesm2 , filename=paste('/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/stacks/',
                                      tag, '/CNRMesm2/', names(CNRMesm2 ), sep = ""), bylayer=TRUE,format="GTiff")
writeRaster(IPSL, filename=paste('/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/stacks/',
                                 tag, '/IPSL/', names(IPSL), sep = ""), bylayer=TRUE,format="GTiff")
writeRaster(MirocES2L, filename=paste('/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/stacks/',
                                      tag, '/MirocES2L/', names(MirocES2L), sep = ""), bylayer=TRUE,format="GTiff")
writeRaster(Miroc6, filename=paste('/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/stacks/',
                                   tag, '/Miroc6/', names(Miroc6), sep = ""), bylayer=TRUE,format="GTiff")
writeRaster(MRI, filename=paste('/Users/danielfurman/Desktop/work/Xantusia_Harvey_Mudd/stacks/',
                                tag, '/MRI/', names(MRI), sep = ""), bylayer=TRUE,format="GTiff")