#
#
# This is a script for identifying which las cells contain the LiDAR data needed
# To compare the UAV predictions out of set
#
#

library(lidR)
library(sp)
library(rgdal)
library(raster)
library(rgeos)

ctg1 <- catalog("/Volumes/HARAPAN/Remote_Sensing/Original/2014_LiDAR_Bioclime/LAS1")
ctg2 <- catalog("/Volumes/HARAPAN/Remote_Sensing/Original/2014_LiDAR_Bioclime/LAS2")

boundary<-readOGR("data/shape/Boundary_outer.shp")

crs(ctg1)<-crs(ctg2)<-crs(boundary)

plot(ctg1)
plot(ctg2)


