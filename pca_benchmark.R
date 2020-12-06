# ---
# title: "pca_jacard_analysis"
# author: "Daniel Furman"
# date: "2020"
# ---
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
library(factoextra)
library(plotly)
library(geometry)

## ------------------------------------------------------------------------


train.files <- list.files(
  path='data/train-rasters-2.5m', full.names=TRUE)
train.layers <- stack(train.files)
names(train.layers)

underlying.pca.data <- as.matrix(train.layers)
head(underlying.pca.data)
underlying.pca.data <- underlying.pca.data[
  complete.cases(underlying.pca.data),]
length(underlying.pca.data[,1])
underlying.pca = prcomp(underlying.pca.data, center = T, scale. = T)
underlying.pca
# screeplot
fviz_eig(underlying.pca)

# visualizing variables
fviz_pca_var(underlying.pca,
             col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T)
## ------------------------------------------------------------------------

xv_baseline <- read.csv('data/collinearity/raw_data.csv')

xv.pca <- scale(xv_baseline, underlying.pca$center,
                underlying.pca$scale) %*% underlying.pca$rotation 

yb_baseline <- read.csv('data/collinearity/raw_data_yb.csv')

yb.pca <- scale(yb_baseline, underlying.pca$center,
                underlying.pca$scale) %*% underlying.pca$rotation 

head(xv.pca)
head(yb.pca)

both <- rbind(xv.pca,yb.pca)
label <- c(rep(0, nrow(xv.pca)), rep(1, nrow(yb.pca))) 
both <- data.frame(both)
both$X <- label
head(both) 

pca_vis.1 <- plot_ly(both, x = ~PC1, y = ~PC2, z = ~PC3, color = ~X,
                   colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))

pca_vis.1 

hull_xan1 = convhulln(both[both$X == 0,2:4], output.options = 'FA')
hull_xan1$vol #vol = 92


## ------------------------------------------------------------------------

# The rest of thise file requires more data to complete,
# data that is not included in the GitHub repo. The future climate data
# was too large to share open-spource

names <- c('BCC', 'Can', 'CNRMcm6','CNRMesm2', 'IPSL', 'Miroc6', 
           'MirocES2L', 'MRI')

for (i in names){ 
  train.files <- list.files(
    path=paste('data/CMIP6/ssp370_2061-2080/', i, sep = ''),
    full.names=TRUE, pattern='tif')
  stack = stack(train.files)
  write.csv(extract(stack, acsel),
    paste('data/CMIP6/ssp370_2061-2080/', i, '/xv_pca.csv', sep = ''),
    row.names=FALSE)
  write.csv(extract(stack, yb_final),
    paste('data/CMIP6/ssp370_2061-2080/', i, '/yb_pca.csv', sep = ''),
    row.names=FALSE)
}  
  
for (i in names){ 
  train.files <- list.files(
    path=paste('data/CMIP6/ssp370_2081-2100/', i, sep = ''),
    full.names=TRUE, pattern='tif')
  stack = stack(train.files)
  write.csv(extract(stack, acsel),
            paste('data/CMIP6/ssp370_2081-2100/', i, '/xv_pca.csv', sep = ''),
            row.names=FALSE)
  write.csv(extract(stack, yb_final),
            paste('data/CMIP6/ssp370_2081-2100/', i, '/yb_pca.csv', sep = ''),
            row.names=FALSE)
}  

for (i in names){ 
  train.files <- list.files(
    path=paste('data/CMIP6/ssp370_2041-2060/', i, sep = ''),
    full.names=TRUE, pattern='tif')
  stack = stack(train.files)
  write.csv(extract(stack, acsel),
            paste('data/CMIP6/ssp370_2041-2060/', i, '/xv_pca.csv', sep = ''),
            row.names=FALSE)
  write.csv(extract(stack, yb_final),
            paste('data/CMIP6/ssp370_2041-2060/', i, '/yb_pca.csv', sep = ''),
            row.names=FALSE)
} 

#pause
avg_xv_41 <- read.csv('data/CMIP6/ssp370_2041-2060/average_PCA_xv.csv')
avg_xv_61 <- read.csv('data/CMIP6/ssp370_2061-2080/average_PCA_xv.csv')
avg_xv_81 <- read.csv('data/CMIP6/ssp370_2081-2100/average_PCA_xv.csv')

avg_xv_41 <- avg_xv_41[,2:20]
avg_xv_41 <- scale(avg_xv_41, underlying.pca$center,
                underlying.pca$scale) %*% underlying.pca$rotation 

avg_xv_61 <- avg_xv_61[,2:20]
avg_xv_61 <- scale(avg_xv_61, underlying.pca$center,
                   underlying.pca$scale) %*% underlying.pca$rotation 

avg_xv_81 <- avg_xv_81[,2:20]
avg_xv_81 <- scale(avg_xv_81, underlying.pca$center,
                   underlying.pca$scale) %*% underlying.pca$rotation 

both <- rbind(xv.pca, avg_xv_41, avg_xv_61, avg_xv_81)
label <- c(rep(0, nrow(xv.pca)), rep(1, nrow(avg_xv_41)), rep(2, nrow(avg_xv_61)),
           rep(3, nrow(avg_xv_81))) 
both <- data.frame(both)
both$X <- label
head(both)

which.max(xv.pca[,2])
length(xv.pca[,2])
pca_vis.3 <- plot_ly(both, x = ~PC1, y = ~PC2, z = ~PC3, color = ~X,
                     colors = c( '#000000', '#db675b' , '#d24031', '#952c21' )) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))

pca_vis.3
hull_xv1 = convhulln(xv.pca[,2:4], output.options = 'FA')
hull_xv1$vol #vol = 85

avg_yb_41 <- read.csv('data/CMIP6/ssp370_2041-2060/average_PCA_yb.csv')
avg_yb_61 <- read.csv('data/CMIP6/ssp370_2061-2080/average_PCA_yb.csv')
avg_yb_81 <- read.csv('data/CMIP6/ssp370_2081-2100/average_PCA_yb.csv')

avg_yb_41 <- avg_yb_41[,2:20]
avg_yb_41 <- scale(avg_yb_41, underlying.pca$center,
                   underlying.pca$scale) %*% underlying.pca$rotation 

avg_yb_61 <- avg_yb_61[,2:20]
avg_yb_61 <- scale(avg_yb_61, underlying.pca$center,
                   underlying.pca$scale) %*% underlying.pca$rotation 

avg_yb_81 <- avg_yb_81[,2:20]
avg_yb_81 <- scale(avg_yb_81, underlying.pca$center,
                   underlying.pca$scale) %*% underlying.pca$rotation 

both <- rbind(yb.pca, avg_yb_41, avg_yb_61, avg_yb_81)
label <- c(rep(0, nrow(yb.pca)), rep(1, nrow(avg_yb_41)), rep(2, nrow(avg_yb_61)),
           rep(3, nrow(avg_yb_81))) 
both <- data.frame(both)
both$X <- label
head(both)

pca_vis.4 <- plot_ly(both, x = ~PC1, y = ~PC2, z = ~PC3, color = ~X,
                     colors = c( '#000000', '#db675b' , '#d24031', '#952c21' )) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))

pca_vis.4

hull_yuc1 = convhulln(yb.pca[,2:4], output.options = 'FA')
hull_yuc1$vol #vol = 83
# extra tutorial on jacard similarity statistic

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/clusterIndex.R") 
library(cluster)
y <- matrix(rnorm(5000), 1000, 5, dimnames=list(paste("g", 1:1000, sep=""), paste("t", 1:5, sep="")))
head(y)
?clara
clarax <- clara(y, 49)
clarax$clustering
clV1 <- clarax$clustering
clarax <- clara(y, 50)
clV2 <- clarax$clustering
ci <- cindex(clV1=clV1, clV2=clV2, self=FALSE, minSZ=1, method="jaccard")
ci[2:3]

length(clV1)
length(clV2)


