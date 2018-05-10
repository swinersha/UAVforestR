rm(list=ls())

library(raster)
library(maptools)
library(rgdal)
library(ggplot2)
library(rgeos)
library(foreach)
library(doParallel)
#library(gridExtra)
#library(splancs)
#library(RColorBrewer)
#library(tidyverse)
#library(UAVforestR)
#library(lmtest)
#library(car)

# Load R source files
R_source_files<-list.files(
  path = "R",
  pattern = "*.R$",
  full.names = TRUE
)
sapply(R_source_files, function(x) source(x, local = FALSE,  echo = FALSE))

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
man_lid$Area<-as.numeric(as.character(man_lid$Area))
man_lid$R<-sqrt(man_lid$Area/pi) # Calculate a pseudoradius
man_lid$lid_CH_mean<-raster::extract(lid_chm, man_lid, fun=mean)
man_lid$lid_CH_max<-raster::extract(lid_chm, man_lid, fun=max)
lut_lid<-rq_lut(x=man_lid$lid_CH_max, y=man_lid$R, log=TRUE)

man_uav<-readOGR("data/shape/Manual trees UTM.shp")
man_uav$Area<-as.numeric(as.character(man_uav$Area))
man_uav$R<-sqrt(man_uav$Area/pi) # Calculate a pseudoradius
man_uav$uav_CH_mean<-raster::extract(uav_chm, man_uav, fun=mean)
man_uav$uav_CH_max<-raster::extract(uav_chm, man_uav, fun=max)
lut_uav<-rq_lut(x=man_uav$uav_CH_max, y=man_uav$R, log=TRUE)



ggplot(man_lid@data, aes(x=lid_CH_max, y=R)) +
  geom_point() +
  geom_smooth(method="lm") +
  xlim(0, 50) + ylim(0,10)

# Load the UAV and LiDAR images
uav_chm_blur<-raster("data/raster/uav_chm_blur.tif")
uav_chm_sobel<-raster("data/raster/uav_chm_sobel.tif")

# lid_chm_blur<-blur(lid_chm)
# lid_chm_sobel<-sobel_edge(lid_chm)
# writeRaster(lid_chm_blur, "data/raster/lid_chm_blur.tif")
# writeRaster(lid_chm_sobel, "data/raster/lid_chm_sobel.tif")

lid_chm_blur<-raster("data/raster/lid_chm_blur.tif")
lid_chm_sobel<-raster("data/raster/uav_chm_sobel.tif")


# Make a parameter matrix:
THRESHSeed_vec<-seq(from=0.3, to=0.9, by=0.1)
THRESHCrown_vec<-seq(from=0.3, to=0.9, by=0.1)
SOBELstr_vec<-seq(from=1.5, to=6, by=1)
params<-expand.grid(THRESHSeed=THRESHSeed_vec,
                    THRESHCrown=THRESHCrown_vec,
                    SOBELstr=SOBELstr_vec)

htod<-function(x, tau) htod_lookup(x, lut=lut_uav, tau)

cl<-parallel::makeCluster(25)
registerDoParallel(cl)

foreach(i = 1:nrow(params),
        .packages = c("raster", "rgeos", "sp"),
        .inorder = TRUE) %dopar% {
          itc<-itcIMG_fast(uav_chm_blur,
                      uav_chm_sobel,
                      THRESHSeed=params$THRESHSeed[i], THRESHCrown=params$THRESHCrown[i],
                      tau = 90,
                      specT=0,
                      SOBELstr=params$SOBELstr[i],
                      lm.searchwin = 3,
                      #pypath = "gdal_polygonize.py"
                      pypath = "/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_polygonize.py"
          )
          # file_name<-paste('data/shape/ITC_trees_params_cost/seed_', params[i,1],
          #                  '_crown_', params[i,2],
          #                  '_sobel_', params[i,3], sep='')
          file_name<-paste('data/shape/ITC_trees_params_cost_uav/seed_', params[i,1],
                           '_crown_', params[i,2],
                           '_sobel_', params[i,3],
                           '_specT_0', sep='')
          writeOGR(itc,
                   dsn = paste(file_name, '.shp', sep=''),
                   layer = basename(file_name),
                   drive = 'ESRI Shapefile')
        }
stopCluster(cl)

