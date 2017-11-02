
#' Creates polygon boundaries for grid cells
#'
#' @param x x coordinates as numeric vector
#' @param y y coordinates as numeric vector
#' @param size the grid size as a number
#' @return A matrix of x and y coordinates
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 16-10-10
#'


grid.poly<-function(x, y, size)
{
  size<-size/2 # because the xy positions are in the centre, the size is halfed.
  x<- x + c(-size, -size, size, size, -size)
  y<- y + c(-size, size, size, -size, -size)
  xy<-cbind(x=x, y=y)
  return(xy)
}

#' Converts an ordinary polygon to a SpatialPolygon
#'
#' @param poly A matrix describing a polygon
#' @param ID The ID to be handed to sp::Polygons
#' @return A Polygons object
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 16-10-10
#'
#'
sp.poly<-function(poly, ID=NULL)
{
  poly<-Polygon(data.frame(poly))
  poly<-Polygons(list(poly), ID=ID)
  return(poly)
}


#' Create a SpatialPolygons object for all xy combinations
#'
#' @param xy a matrix of xy coordinates
#' @param size The desired size of the polygons
#' @param ID The ID to be handed to sp::Polygons
#' @param proj4string The CRS of the polygons
#' @return A spatialPolygons object
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 16-10-10
#'
#'

polygrid<-function(xy, size, ID=NULL, proj4string)
{
  if(is.null(ID))
    ID<-paste('cell', 1:nrow(xy), sep='')

  xy<-cbind(xy, ID)
  pg<-apply(xy, 1, function(CELL)
    {
      poly<-grid.poly(x=as.numeric(CELL['x']), y=as.numeric(CELL['y']), size=size)
      poly<-sp.poly(poly, ID=CELL['ID'])
    }
  )
  pg<-SpatialPolygons(pg, 1:length(pg), proj4string = proj4string)
  return(pg)
}


#' A helper function which builds a grid of sample points.
#'
#' @param imagery A raster image
#' @param grid_size The size of the grid to make
#' @return A matrix of xy coordinates
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 16-10-10
#'

img_grid<-function(imagery, grid_size)
{
  minx_grid<-extent(imagery)@xmin
  miny_grid<-extent(imagery)@ymin
  maxx_grid<-extent(imagery)@xmax
  maxy_grid<-extent(imagery)@ymax

  n_x<-(maxx_grid-minx_grid) %/% grid_size
  n_y<-(maxy_grid-miny_grid) %/% grid_size

  x_grid<-seq(from=minx_grid, by=grid_size, length=n_x)
  y_grid<-seq(from=miny_grid, by=grid_size, length=n_y)

  xy<-expand.grid(x=x_grid, y=y_grid)
}


#' Build a DTM by extracting the minimum values within grid cells
#' These are then averaged using a 2D GAM.
#'
#' @param imagery A raster image
#' @param grid_by The size of the grid to make to extract values
#' @param k The number of smoothing knots for the GAM
#' @param predict_grid_by The resolution of the raster to predict for
#' @return A matrix of xy coordinates
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 16-10-10
#'

build_dtm<-function(imagery, grid_by, k, predict_grid_by, ...)
{
  grid_by<-round(grid_by/res(imagery)[1])
  row_ind<-0:((nrow(imagery) %/% grid_by))*grid_by +1
  col_ind<-0:((ncol(imagery) %/% grid_by))*grid_by +1
  row_ind<-row_ind[row_ind<nrow(imagery)] # excludes any values that are out of range
  col_ind<-col_ind[col_ind<ncol(imagery)] # ...
  xy_grid<-expand.grid(row_ind, col_ind)

  min_vals_x<-vector('integer', length=nrow(xy_grid))
  min_vals_y<-vector('integer', length=nrow(xy_grid))
  min_vals<-vector('numeric', length=nrow(xy_grid))
#min_vals<-matrix(0, nrow=nrow(xy_grid), ncol=3)
  cat("Extracting gridded values\n")
  for(i in 1:nrow(xy_grid))
  {
    prog.seq<-seq(0,1, by=0.1)
    progress<-round(nrow(xy_grid) * prog.seq)  %in% i
    if(any(progress))
      cat(paste((which(progress)-1)*10), '%   ', sep='')
    # Ensure the grid doesn't exceed the size of the raster:
    if((xy_grid[i,2] + grid_by) > dim(imagery)[2])
      grid_by_x<-(ncol(imagery)+1 %% grid_by) -1
    else
      grid_by_x<-grid_by
    if((xy_grid[i,1] + grid_by) > dim(imagery)[1])
      grid_by_y<-(nrow(imagery)+1 %% grid_by) -1
    else
      grid_by_y<-grid_by

      # Pull out the values for the block
    blk<-getValuesBlock(imagery, row=xy_grid[i,1], nrows=grid_by_y, col=xy_grid[i,2], ncols=grid_by_x, format='matrix')
    if(all(is.na(blk))) # Make sure there are some values:
    {
      min_vals_x[i]<-NA
      min_vals_y[i]<-NA
      min_vals[i]<-NA
    }
    else # Find the lowest value and it's position:
    {
      min_vals[i]<-min(blk, na.rm=TRUE)
      min_pos<-which(blk==min_vals[i], arr.ind=TRUE)
      min_vals_y[i]<-min_pos[1,1]+xy_grid[i,1]
      min_vals_x[i]<-min_pos[1,2]+xy_grid[i,2]
    }
  }

  min_vals_y<-extent(imagery)@ymax - (min_vals_y* res(imagery)[1]) # These resolution values could be the wrong way round!!!
  min_vals_x<-extent(imagery)@xmin + (min_vals_x* res(imagery)[2])
  min_vals<-data.frame(x=min_vals_x, y=min_vals_y, z=min_vals)

  cat("\nFitting GAM\n")
  m.dtm <- mgcv::gam(z ~ s(x, y, k=k), data=min_vals) # fit the 2D spline.

  # Create a raster with the predicted values
  grid<-img_grid(imagery, grid_size=predict_grid_by) # produce a grid for prediction.
  dtm_pred<-mgcv::predict.gam(m.dtm, newdata=grid) # predict
  dtm_pred<-cbind(grid, dtm_pred)
  sp::coordinates(dtm_pred)=~x+y # convert to gridded spatial object
  sp::proj4string(dtm_pred)<-raster::crs(imagery)
  sp::gridded(dtm_pred) <-TRUE
  dtm_pred <- raster::raster(dtm_pred) # convert to raster
  dtm_pred<-resample(dtm_pred, imagery, ...)
  return(dtm_pred)
}

#' A function to build a DTM by extracting the minimum values within grid cells
#' These are then averaged using a 2D GAM.
#'
#' @param imagery A raster image
#' @param grid_by The size of the grid to make to extract values
#' @param type Either "raster" or "grid" the former returns a raster, the latter a matrix
#' @param fun A function to summarise the points in the grid
#' @return Either a raster or matrix depending on the type argument entered
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 16-10-10
#'


grid_fun<-function(imagery, grid_by, type='raster', fun)
{
  grid_by_map<-grid_by # stored for conversion back to map scale
  grid_by<-round(grid_by/res(imagery)[1])
  row_ind<-0:((nrow(imagery) %/% grid_by))*grid_by +1
  col_ind<-0:((ncol(imagery) %/% grid_by))*grid_by +1
  row_ind<-row_ind[row_ind<nrow(imagery)] # excludes any values that are out of range
  col_ind<-col_ind[col_ind<ncol(imagery)] # ...

  xy_grid<-expand.grid(row_ind, col_ind)

  y<-vector('numeric', length=nrow(xy_grid))
  #min_vals<-matrix(0, nrow=nrow(xy_grid), ncol=3)
  cat("Extracting gridded values\n")
  for(i in 1:nrow(xy_grid))
  {
    prog.seq<-seq(0,1, by=0.1)
    progress<-round(nrow(xy_grid) * prog.seq)  %in% i
    if(any(progress))
      cat(paste((which(progress)-1)*10), '%   ', sep='')
    # Ensure the grid doesn't exceed the size of the raster:
    if((xy_grid[i,2] + grid_by) > dim(imagery)[2])
      grid_by_x<-(ncol(imagery)+1 %% grid_by) -1
    else
      grid_by_x<-grid_by
    if((xy_grid[i,1] + grid_by) > dim(imagery)[1])
      grid_by_y<-(nrow(imagery)+1 %% grid_by) -1
    else
      grid_by_y<-grid_by

    # Pull out the values for the block
    blk<-getValuesBlock(imagery, row=xy_grid[i,1], nrows=grid_by_y, col=xy_grid[i,2], ncols=grid_by_x, format='matrix')
    if(all(is.na(blk))) # Make sure there are some values:
    {
      y[i]<-NA
    }
    else # Find the lowest value and it's position:
    {
      y[i]<-fun(blk)
    }
  }

  # Create a raster with the predicted values
  # Sets xy_grid to the map scale:
  xy_grid<-expand.grid(1:length(row_ind), 1:length(col_ind))
  colnames(xy_grid)<-c('y', 'x')
  xy_grid$x<-xy_grid$x*grid_by_map+extent(imagery)@xmin
  xy_grid$y<-extent(imagery)@ymax-xy_grid$y*grid_by_map
  xy_grid<-data.frame(x=xy_grid$x, y=xy_grid$y, y)
  sp::coordinates(xy_grid)=~x+y # convert to gridded spatial object
  sp::proj4string(xy_grid)<-raster::crs(imagery)
  sp::gridded(xy_grid) <-TRUE
  if(type=='grid')
    return(xy_grid)
  xy_raster <- raster::raster(xy_grid) # convert to raster
  if(type=='raster')
    return(xy_raster)
}

