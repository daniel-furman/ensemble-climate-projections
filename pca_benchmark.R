#' ---
#' title: "Principal Component Analysis for SICB 2021"
#' author: "Daniel Furman"
#' date: "June 26, 2019"
#' output: html_document
#' ---
#' 
#' 
## -----------------------------------------------------------------------------------------------------

library(raster)
library(plot3D)
library(plotly)
library(factoextra)
library(sets)
library(geometry)
library(rgl)
library(tidyr)

## -----------------------------------------------------------------------------------------------------

pca_df_total = read.csv('data/total_pre_pca.csv', header = TRUE) #climate values from all 19 BIOCLIM
pca_clim_total = prcomp(pca_df_total, center = T, scale. = T) # PCA dimension analysis for the 19

# screeplot
fviz_eig(pca_clim_total) # we see that the first two component contain a substantial amount of the total
# visualizing variables
fviz_pca_var(pca_clim_total,
             col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = T) #We see that precip and temp variables do split 

total_pca_complete = read.csv('data/total_pca.csv', header = TRUE)

## -----------------------------------------------------------------------------------------------------

p_total <- plot_ly(total_pca_complete, x = ~PC1, y = ~PC2, z = ~PC3, color = ~X,
  colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                     yaxis = list(title = 'PC2'),
                     zaxis = list(title = 'PC3')))
p_total

## -----------------------------------------------------------------------------------------------------
xan_set1 <-
  total_pca_complete %>%
  filter(X == 0) 
xan_set1 = xan_set1[,2:4]
yucca_set1 <-
  total_pca_complete %>%
  filter(X == 1) 

yucca_set1 = yucca_set1[,2:4]

hull_xan1 = convhulln(xan_set1, output.options = 'FA')
hull_xan1$vol #vol = 531.6
hull_yuc1 = convhulln(yucca_set1, output.options = 'FA')
hull_yuc1$vol #vol = 516.9

xanmesh = to.mesh3d(hull_xan1)
yucmesh = to.mesh3d(hull_yuc1)

wire3d(xanmesh, meshColor = c("faces"), color = "red")
dot3d(xanmesh, meshColor = c("faces"), color = "red")
wire3d(yucmesh, meshColor = c("faces"), color = "blue")
dot3d(yucmesh, meshColor = c("faces"), color = "blue")
