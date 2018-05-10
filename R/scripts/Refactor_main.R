rm(list=ls())

library(raster)
library(maptools)
library(rgdal)
library(ggplot2)
library(gridExtra)
library(splancs)
library(RColorBrewer)
library(rgeos)
library(tidyverse)
#library(lmtest)
#library(car)

# Load the UAV and LiDAR images
uav_dsm<-raster("data/raster/uav_dsm_matched_cropped.tif")
uav_dtm<-raster("data/raster/uav_dtm_matched_cropped.tif")
lid_dsm<-raster("data/raster/lid_dsm_matched_cropped.tif")
lid_dtm<-raster("data/raster/lid_dtm_matched_cropped.tif")
lid_chm<-raster("data/raster/lid_chm_matched_cropped.tif")

# Ensure the UAV and LiDAR DEMs are vertically aligned (numbers from manual alignment):
uav_dsm<-uav_dsm-17.6
uav_dtm<-uav_dtm-17.6
# remove very low values from UAV DSM
uav_dsm[uav_dsm<0]<-NA
uav_dtm[uav_dtm<0]<-NA
# Calculates the canopy height models for the UAV
uav_chm<-uav_dsm-uav_dtm # UAV DTM
# Calculate the 90th percentile of the manual crowns height predicted from radius:
man_lid<-readOGR("data/shape/Manual trees LiDAR UTM.shp")
man_lid$R<-sqrt(man_lid$Area/pi) # Calculate a pseudoradius
man_lid$lid_CH_mean<-raster::extract(lid_chm, man_lid, fun=mean)
man_lid$lid_CH_max<-raster::extract(lid_chm, man_lid, fun=max)
lut_lid<-rq_lut(x=man_lid$lid_CH_max, y=man_lid$R, log=TRUE)

ggplot(man_lid@data, aes(x=lid_CH_max, y=R)) +
  geom_point() +
  geom_smooth(method="lm") +
  xlim(0, 50) + ylim(0,10)

# Load the UAV and LiDAR images
uav_chm_blur<-raster("data/raster/uav_chm_blur_e1.tif")
uav_chm_sobel<-raster("data/raster/uav_chm_sobel_e1.tif")

crowns_sp<-itcIMG_fast(uav_chm_blur,
            uav_chm_sobel,
            THRESHSeed = 0.4,
            THRESHCrown = 0.4,
            tau = 90,
            specT=2,
            SOBELstr=4,
            lm.searchwin = NULL,
            pypath = "/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_polygonize.py"
)


plot(uav_chm_blur, asp=1)
plot(crowns_sp, add=TRUE)

ggplot(crowns_sp@data, aes(x = treeHeights_max, sobelbound)) +geom_point() +ylim(0,1)
ggplot(crowns_sp@data, aes(x = treeHeights_max, hbound_mean)) +geom_point() +ylim(0,1)
ggplot(crowns_sp@data, aes(x = treeHeights_max, hbound_max)) +geom_point() +ylim(0,1)
ggplot(crowns_sp@data, aes(x = treeHeights_max, allombound)) +geom_point() +ylim(0,1)
ggplot(crowns_sp@data, aes(x = treeHeights_max, crownbound)) +geom_point() +ylim(0,1)
ggplot(crowns_sp@data, aes(x = treeHeights_max, tallerbound)) +geom_point() +ylim(0,1)

# Make a parameter matrix:
THRESHSeed_vec<-seq(from=0.3, to=0.9, by=0.2)
THRESHCrown_vec<-seq(from=0.3, to=0.9, by=0.2)
SOBELstr_vec<-seq(from=1.5, to=6, by=2)
params<-expand.grid(THRESHSeed=THRESHSeed_vec,
                    THRESHCrown=THRESHCrown_vec,
                    SOBELstr=SOBELstr_vec)
