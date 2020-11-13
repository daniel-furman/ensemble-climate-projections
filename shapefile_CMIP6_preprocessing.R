# ---
# title: "Pre-processing for projections"
# author: "Daniel Furman"
# date: "6/27/2020"
# output: html_document
# ---
  
# This R script processes CMIP6 TIF files and saves them to file. 

### Load Climate Files 

## ------------------------------------------------------------------------

library(dismo)
library(rgdal)
library(DescTools)

tag <- 'ssp370_2081-2100' 
#change to your own file path
files <- list.files(path=paste('/Volumes/HardDrive/raw_CMIP6/',
                    tag, sep = ""), pattern='tif', full.names=TRUE ) 
print('These are the raw tif files:')
files
# dictionary of model names
# {'BCC-CSM2-MR', 'CanESM5','CNRM-CM6-1','CNRM-ESM2-1’, 'IPSL-CM6A-LR'
# ,'MIROC-ES2L','MIROC6', and 'MRI-ESM2-0'}

## ------------------------------------------------------------------------

### Stacking Climate Rasters 
# Now that everything is loaded into R, we  create rasters stacks and save
# them as new tif files. These rasters are now all less than 1.5 mb, and
# can be uploaded to GitHub. 

# Initialize empty model vector
BCC_stack_layers <- stack()
Can_stack_layers <- stack()
CNRMcm6_stack_layers <- stack()
CNRMesm2_stack_layers <- stack()
IPSL_stack_layers <- stack()
MirocES2L_stack_layers <- stack()
Miroc6_stack_layers <- stack()
MRI_stack_layers <- stack()

ext <- extent(-122.767, -109.3192, 31.91825, 38.1922)

model_stacks <- c(BCC_stack_layers,Can_stack_layers,CNRMcm6_stack_layers,
                  CNRMesm2_stack_layers,IPSL_stack_layers,
                  MirocES2L_stack_layers, Miroc6_stack_layers,
                  MRI_stack_layers)

# populate raster stacks for each GCM model and their 19 bands
for (i in c(1:8)){ 
  for (l in c(1:19)){
    bclim <- raster(files[[i]], band = l) 
    bclim <- crop(bclim, ext)
    model_stacks[[i]] <- stack(model_stacks[[i]], bclim)
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

## ------------------------------------------------------------------------

# Write Rasters
writeRaster(BCC, filename=paste(
  'data/CMIP6/', tag, '/BCC/', names(BCC), sep = ""), bylayer=TRUE,
  format="GTiff", overwrite = TRUE)
writeRaster(Can, filename=paste(
  'data/CMIP6/', tag, '/Can/', names(Can), sep = ""), bylayer=TRUE,
  format="GTiff", overwrite = TRUE)
writeRaster(CNRMcm6, filename=paste(
  'data/CMIP6/', tag, '/CNRMcm6/', names(CNRMcm6), sep = ""), bylayer=TRUE,
  format="GTiff", overwrite = TRUE)
writeRaster(CNRMesm2 , filename=paste(
  'data/CMIP6/', tag, '/CNRMesm2/', names(CNRMesm2 ), sep = ""), bylayer=TRUE,
  format="GTiff", overwrite = TRUE)
writeRaster(IPSL, filename=paste(
  'data/CMIP6/', tag, '/IPSL/', names(IPSL), sep = ""), bylayer=TRUE,
  format="GTiff", overwrite = TRUE)
writeRaster(MirocES2L, filename=paste(
  'data/CMIP6/', tag, '/MirocES2L/', names(MirocES2L), sep = ""), bylayer=TRUE,
  format="GTiff", overwrite = TRUE)
writeRaster(Miroc6, filename=paste(
  'data/CMIP6/', tag, '/Miroc6/', names(Miroc6), sep = ""), bylayer=TRUE,
  format="GTiff", overwrite = TRUE)
writeRaster(MRI, filename=paste(
  'data/CMIP6/', tag, '/MRI/', names(MRI), sep = ""), bylayer=TRUE,
  format="GTiff", overwrite = TRUE)
