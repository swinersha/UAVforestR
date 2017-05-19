# A function for estimating a DTM from the trees segmented from the UAV CHM
# A raster is output 
#
#
# Dr Tom Swinfield
# 17.11.16

#trees=conf.trees_all; model=model; dsm=uav_dsm; chm=uav_chm; predict_grid_by=5; k=100

img_xy_extract<-function(x, y, img){
  xy<-data.frame(x=x, y=y)
  xy<-SpatialPoints(xy, proj4string = crs(img))
  z<-raster::extract(img, xy)
  return(z)
}

#x=conf_auto$maxx; y=conf_auto$maxy; z.est=conf_auto$CH50; img=uav_dsm

dtm_extract<-function(x, y, z.est, img){
  xy<-data.frame(x=x, y=y)
  xy_spatial<-SpatialPoints(xy, proj4string = crs(img))
  z<-raster::extract(img, xy_spatial)
  dtm.est<-z-z.est
  xyz<-cbind(xy, z=dtm.est)
  return(xyz)
}

#x=conf_auto_uav$maxx; y=conf_auto_uav$maxy; z=conf_auto_uav$CH50
#img=uav_dsm;
#predict_grid_by=25;k=100


surface_predict<-function(xyz, img, predict_grid_by, k, type='response'){
  fm.xyz <-mgcv::gam(z ~ s(x, y, k=k), data=xyz) # fit the 2D spline.
  # summary(xyz)
  grid<-img_grid(img, grid_size=predict_grid_by) # produce a grid for prediction.
  if(!any(type %in% c('response', 's.e.')))
    stop("dtm_model_predict: predict must be set to response or s.e.")
  if(type == 'response')
    xyz<-mgcv::predict.gam(fm.xyz, newdata=grid) # predict
  if(type == 's.e.')
    xyz<-mgcv::predict.gam(fm.xyz, newdata=grid, se.fit=TRUE)$se.fit

  xyz<-cbind(grid, xyz)
  sp::coordinates(xyz)=~x+y # convert to gridded spatial object
  sp::proj4string(xyz)<-raster::crs(img)
  sp::gridded(xyz) <-TRUE
  xyz <- raster::raster(xyz) # convert to raster
  xyz<-resample(xyz, img)
  xyz<-mask(xyz, img)
  return(xyz)
}

dtm_predict<-function(x, y, z, img, predict_grid_by, k, type='response'){
  xyz<-dtm_extract(x=x, y=y, z.est=z, img=img)
  # Create a raster with the predicted values
  
  xyz<-surface_predict(xyz, img=img,
            predict_grid_by=predict_grid_by,
            k=k, type=type)
  
  return(xyz)
}


dtm_model_predict<-function(trees, model, dsm, chm, predict_grid_by, k, type='response')
{
  # Extracts the DSM values for the confidently segmented trees:
  #treepos_xy<-data.frame(x=trees$maxx, y=trees$maxy)
  #treepos<-SpatialPoints(treepos_xy, proj4string = crs(dsm))
  #treepos_dsm<-raster::extract(dsm, treepos)
  #treepos_uavchm<-raster::extract(chm, treepos)
  #treepos_lidchm<-raster::extract(lid_chm, treepos)
  
  # Creates the prediction dataframe from the segmented trees:
  segtree<-data.frame(x=trees$maxx, y=trees$maxy, h_uav=trees$CH_max, CR_m=trees$CR_m, sobel=log(trees$sobel+1), hbound=trees$hbound+trees$sobelbound)
  segtree<-segtree[complete.cases(segtree),]
  
  # Makes the prediction  
  segtree$pred<-predict(model, newdata=segtree)
  
  # Plots out the GAM predictions:
  #par(mfrow=c(1,2), mar=c(4,4,2,2))
  #plot(x=segtree$h_uav, y=segtree$pred); abline(0,1, col='red')
  #hist(segtree$h_uav-segtree$pred)
  
  # Extracts the DSM values for the confidently segmented trees:
  treepos_xy<-data.frame(x=segtree$x, y=segtree$y)
  treepos<-SpatialPoints(treepos_xy, proj4string = crs(dsm))
  treepos_dsm<-raster::extract(dsm, treepos)
  treepos_uavchm<-raster::extract(chm, treepos)
  #treepos_lidchm<-raster::extract(lid_chm, treepos)
  treepos_dtm<-treepos_dsm-segtree$pred
  treepos_dtm<-data.frame(treepos_xy, z=treepos_dtm, w=segtree$h_uav) # You can try playing arond with what you weight by ...
  Esc.treepos_dtm<<-treepos_dtm
  # I have tried h_uav, hbound and sobel, they are all pretty similar 
  # but sobel looks best
  
  # Estimating the ground 
  cat("\nFitting GAM\n")
  #require(mgcv)
  est.dtm <-mgcv::gam(z ~ s(x, y, k=k), data=treepos_dtm, weights=w/mean(w)) # fit the 2D spline.
  #summary(est.dtm)
  # Create a raster with the predicted values
  grid<-img_grid(dsm, grid_size=predict_grid_by) # produce a grid for prediction.
  if(!any(type %in% c('response', 's.e.')))
    stop("dtm_model_predict: predict must be set to response or s.e.")
  if(type == 'response')
    dtm_pred<-mgcv::predict.gam(est.dtm, newdata=grid) # predict
  if(type == 's.e.')
    dtm_pred<-mgcv::predict.gam(est.dtm, newdata=grid, se.fit=TRUE)$se.fit
  
  dtm_pred<-cbind(grid, dtm_pred)
  sp::coordinates(dtm_pred)=~x+y # convert to gridded spatial object
  sp::proj4string(dtm_pred)<-raster::crs(chm)
  sp::gridded(dtm_pred) <-TRUE
  dtm_pred <- raster::raster(dtm_pred) # convert to raster
  dtm_pred<-resample(dtm_pred, chm)
  dtm_pred<-mask(dtm_pred, chm)
  return(dtm_pred)
}
