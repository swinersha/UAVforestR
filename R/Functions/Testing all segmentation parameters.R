rm(list=ls())

library(maptools)
library(rgdal)
library(raster)
library(ggplot2)
library(gridExtra)
library(splancs)
library(RColorBrewer)
library(palettetown)
library(rgeos)
library(caret)
library(dplyr)
library(tidyr)
library(lmtest)
library(car)

source("R/build_dtm.R")
source("R/ITC segment update fast_gobble.R")
source("R/gdal_polygonizeR.R")
source("R/dtm_model_predict.R")
source("R/ggmultiplot.R")
source("R/DEMderiv.R")
source('R/min_extent.R')
source('R/poly_rast_vals.R')
source('R/crown_overlap.R')
source("R/rq_lut.R") # The function to run the quantile regression
source("R/Image preprocessing tools.R") # The function to run the quantile regression

# Load the UAV and LiDAR images
uav_dsm<-raster("data/uav_dsm_matched_cropped.tif")
uav_dtm<-raster("data/uav_dtm_matched_cropped.tif")

uav_dsm<-uav_dsm-17.6
uav_dtm<-uav_dtm-17.6
# remove very low values from UAV DSM
uav_dsm[uav_dsm<0]<-NA
uav_dtm[uav_dtm<0]<-NA
# Calculates the canopy height models for the UAV
uav_chm<-uav_dsm-uav_dtm # UAV DTM

e1<-extent(313700, 313900, 9746700, 9746900)
uav_chm_e1<-crop(uav_chm, e1)

# Blurring and finding the sobel edges:
uav_chm_blur<-blur(uav_chm)
uav_chm_sobel<-sobel_edge(uav_chm)

# Make a parameter matrix:
THRESHSeed_vec<-seq(from=0.5, to=0.9, by=0.1)
THRESHCrown_vec<-seq(from=0.5, to=0.9, by=0.1)
SOBELstr_vec<-seq(from=0.5, to=6, by=1)
params<-expand.grid(THRESHSeed=THRESHSeed_vec,
              THRESHCrown=THRESHCrown_vec, 
              SOBELstr=SOBELstr_vec)

# Calculate the 90th percentile of the manual crowns height predicted from radius:
man_lid<-readShapeSpatial("/Users/Tom/Documents/Work/R projects/UAVforestR/data/Manual trees LiDAR UTM.shp")
lut_lid<-rq_lut(x=man_lid$lid_CH_max, y=man_lid$R, log=FALSE)

for(i in 1:nrow(params)){
  print(i); print(Sys.time())
  ut<-itcIMG_fast(uav_chm_blur, im_sobel=uav_chm_sobel, 
                THRESHSeed=params$THRESHSeed[i], THRESHCrown=params$THRESHCrown[i], 
                htod=function(x) htod_lookup(x, lut=lut_lid, tau=90),
                specT=2, 
                SOBELstr=params$SOBELstr[i],
                lm.searchwin=3,
                gobble='off', cantext='off')

  file_name<-paste('Data/ITC trees_params/seed_', params[i,1],
                   '_crown_', params[i,2],
                   '_sobel_', params[i,3], sep='')
  writeSpatialShape(ut, file_name)
}

#------------------------------------------------
#
# Finding the trees with the greatest overlap and 
# pulling out their segmentation parameters
# 
#------------------------------------------------

manual<-readShapeSpatial("/Users/Tom/Documents/Work/R projects/UAVforestR/data/Manual trees UTM.shp")
crs(manual)<-crs(uav_chm)

out<-matrix(0, nrow=nrow(params), ncol=6)
colnames(out)<-c('manual.id', 'mean.ratio', 'sd.ratio', 'mean.overlap', 'sd.overlap', 'n.0s')

for(i in 76:nrow(params)){
  print(i); print(Sys.time())
  file_name<-paste('Data/ITC trees_params/seed_', params[i,1],
                   '_crown_', params[i,2],
                   '_sobel_', params[i,3], '.shp', sep='')
  print(file_name)
  ut<-readShapeSpatial(file_name)
  crs(ut)<-crs(uav_chm)
  # Find the crowns with the greatest overlap:
  matched<-crown_overlap(auto_trees=ut, manual_trees=manual, buffer_by=60)
  
  # Summarises the quality of the match:
  out[i,1]<-i
  out[i,2]<-mean(matched$size_ratio)
  out[i,3]<-sd(matched$size_ratio)
  out[i,4]<-mean(matched$overlap_auto_tree)
  out[i,5]<-sd(matched$overlap_auto_tree)
  out[i,6]<-sum(matched$overlap_auto_tree==0)
  
  # add the auto tree data to the manual trees:
  matched<-matched[!is.na(matched$id_auto_tree),]
  auto<-ut[matched$id_auto_tree,]
  auto@data$shapebound<-rowSums(auto@data[,c('hbound', 'sobelbound')])
  auto@data$spacebound<-rowSums(auto@data[,c('allombound', 'crownbound')])
  auto@data<-auto@data[,c('CH_mean', 'CH_max', 'CR_m', 'shapebound', 'spacebound')]
  matched@data<-cbind(matched@data, auto@data)
  
  par(mfrow=c(1,1), mar=c(0,0,0,0))
  plot(matched)
  plot(auto, add=TRUE, border='red')
  
  # Save the output:
  file_name<-paste('Data/ITC trees_params/Matched/seed_', params[i,1],
                   '_crown_', params[i,2],
                   '_sobel_', params[i,3], sep='')
  writeSpatialShape(matched, file_name)
  print(out[i,])
}

par(mfrow=c(1,2), mar=c(4,4,2,2))
plot(out[,'mean.ratio'], type='l')
points(out[,'mean.overlap'], type='l', col='red')
plot(out[,'sd.ratio'], type='l')
points(out[,'sd.overlap'], type='l', col='red')



out2<-matrix(0, nrow=nrow(params), ncol=6)
colnames(out2)<-c('manual.id', 'mean.ratio', 'sd.ratio', 'mean.overlap', 'sd.overlap', 'n')

for(i in 1:nrow(params)){
  print(i); print(Sys.time())
  file_name<-paste('Data/ITC trees_params/Matched/seed_', params[i,1],
                   '_crown_', params[i,2],
                   '_sobel_', params[i,3], '.shp', sep='')
  print(file_name)
  mt<-readShapeSpatial(file_name)
  mt<-mt[!is.na(mt$shapebound),]
  mt<-mt[mt$shapebound>0.9,]
  out2[i,'manual.id']
  out2[i,'mean.ratio']<-mean(mt$size_ratio)
  out2[i,'sd.ratio']<-sd(mt$size_ratio)
  out2[i,'mean.overlap']<-mean(mt$overlap_au)
  out2[i,'sd.overlap']<-sd(mt$overlap_au)
  out2[i,'n']<-nrow(mt)
}

par(mfrow=c(1,3), mar=c(4,4,2,2))
plot(out2[,'mean.ratio'], type='l')
points(out2[,'mean.overlap'], type='l', col='red')
plot(out2[,'sd.ratio'], type='l')
points(out2[,'sd.overlap'], type='l', col='red')
plot(out2[,'n'], type='l', col='red')



ut2<-readOGR("Data/ITC trees_params/seed_0.5_crown_0.8_sobel_4.5.shp")
ut3<-readOGR("Data/ITC trees_params/seed_0.7_crown_0.7_sobel_3.shp")
ut4<-readOGR("Data/ITC trees_params/seed_0.5_crown_0.5_sobel_3.shp")
crs(ut4)<-crs(ut3)<-crs(ut2)<-crs(uav_chm)

plot(uav_chm_e1)
plot(ut2, add=TRUE)
plot(manual, add=TRUE, col=rgb(0.8, 0, 0, 0.5))

plot(ut2)
plot(manual, add=TRUE, col=rgb(0.8, 0, 0, 0.5))

#plot(ut3, add=TRUE, border='blue')
#plot(ut4, add=TRUE, border='red')
#plot(ut1, add=TRUE, border='red')


# 1. Which parameters give a mean overlap closest to one with the least variance?
# 2. Are there particular parameters that interact to make the best prediction of which 
# trees are most confidently segmented according to overlap?


