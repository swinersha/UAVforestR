# UAVforestR
Measuring forest quality using UAV imagary and structure from motion data

## Raster preparation
[Raster preparation]("R/Raster preparation.R")
Should be run first. It is used to align, trim and clean up the UAV and LiDAR
DEMs.

## Testing all segmentation parameters
[Testing all segmentation parameters]("R/Testting all segmentation parameters.R")
This is an engine for running through a list of reasonable parameters
for ITC segment update fast_gobble.R and calculating their performance 
based on a comparison with manually segmented tree crowns. It uses
the allometry calculated from the LiDAR derived manually segmented crowns.

It should enable the selection of a useful set of parameters specifically 
for identifying trees which can be segmented with a high level of confidence.