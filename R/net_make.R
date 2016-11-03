
#rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}
pythag <- function(b, c) {sqrt(b^2 +c^2)}
pairdist <- function(x1, x2, y1, y2) {pythag(x1-x2, y1-y2)}

# A function to calculate the pairwise distances for adjacent pixels
# Pixels must have already been ordered by distance from the local maximum.
line_pairdist <- function(xy) 
{
  xy_dist<-vector('numeric', length=nrow(xy))
  i<-2
  for(i in 2:length(xy_dist))
    xy_dist[i]<-pairdist(xy[i-1,1], xy[i,1], xy[i-1,2], xy[i,2])
  return(xy_dist)
}


# makes the sample lines for the net:
net_lines<-function(nlines, distmax)
{  
  angles<-0:(nlines-1) * (360/nlines)
  x<-round(sin(deg2rad(angles))*distmax)
  y<-round(cos(deg2rad(angles))*distmax)
  xy<-cbind(x,y)
  lines<-lapply(1:nrow(xy), function(i) rbind(c(0,0), c(xy[i,1], xy[i,2])))
  lines <- lapply(lines, function(x) raster::spLines(x))
  return(lines)
}


#nlines<-12; distmax<-20; im_res<-c(0.33333,0.33333)
net_make<-function(nlines, distmax, im_res)
{
  lines<-net_lines(nlines, distmax)
  
  # makes the raster template for the net:
  rowmax<-round(distmax/im_res[2] - im_res[2]/2) -1
  colmax<-round(distmax/im_res[1] - im_res[1]/2) -1
  ncols<-colmax * 2 + 1 # the number either side of the centre plus the centre
  nrows<-rowmax * 2 + 1 # ... and again for the rows
  xmn<- -((ncols * im_res[1]) / 2); xmx <-(ncols * im_res[1]) / 2
  ymn<- -((nrows * im_res[2]) / 2); ymx <-(nrows * im_res[2]) / 2
  r <- raster::raster(ncols=ncols, nrows=nrows, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, resolution = im_res)
  values(r) <- 1:(nrow(r)*ncol(r))

# Makes integer coordinates for each cell position:
  cart_x<--colmax:colmax
  cart_y<--rowmax:rowmax
#i<-3
# Finds the coordinates of the cells on each line:
  net<-lapply(1:nlines, function(i) 
    {
      line<-lines[[i]]
      cellnums<-raster::extract(r, line, cellnumbers=TRUE)[[1]][,'cell']
      cellind<-arrayInd(cellnums, .dim=c(ncols, nrows))
      x<-cart_x[cellind[,1]]
      y<-cart_y[cellind[,2]]
      xy<-cbind(x, y)
      xy_scale<-cbind(x=x*im_res[2], y=y*im_res[1])
      dist<-sqrt(rowSums(xy_scale^2))
      order.ind<-order(dist)
      x<-x[order.ind]
      y<-y[order.ind]
      dist<-dist[order.ind]
      xy<-cbind(x,y)
      pair.dist<-line_pairdist(xy)
      xy<-data.frame(line=rep(i, nrow(xy)), x=x, y=y, dist=dist, pairdist=pair.dist)
      return(xy)
    }
  )
  par(mfrow=c(1,1), mar=c(4,4,2,2)); plot(0, type='n', xlim=range(cart_x), ylim=range(cart_y, asp=1))
  lapply(net, function(X) points(X[,'x'], X[,'y'], type='l'))
  
  net<-do.call(rbind, net)
  return(net)
}
