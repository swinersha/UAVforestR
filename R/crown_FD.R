
sp.crowns.features<-lapply(1:length(sp.crowns), function (i)
{

poly<-crowns<-sp.crowns



# A function which takes a polygon describing a tree crown alongside a chm and orthomosaic raster
# it returns a load of colour related image features for the crown
# If plot is set to 'on' the raster data will be displayed as a plot
crown_FD<-function(om, chm, poly, plot='on')
{
  om_crown<-spTransform(poly, crs(om)) # Transform the crown polygon to the correct crs.
  chm_crown<-spTransform(crown, crs(chm))
  
  om<-brick(om.tile, package='raster')
  om.crown <- crop(om, om_crown, snap='out') # crops the om
  om.crown.mask <- mask(om.crown, om_crown, inverse=FALSE) # sets values outside of the polygon to NA
  chm.crown <- crop(chm, chm_crown, snap='out') # crops the chm
  chm.crown.mask<- mask(chm.crown, chm_crown, inverse=FALSE) # sets values outside of the polygon to NA
  
  if(plot=='on')
  {
    par(mfrow=c(1,2))
    plotRGB(om.crown.mask)
    plot(crown, add=TRUE)
    plot(chm.crown.mask, axes=FALSE, legend=FALSE, box=FALSE)
    plot(chm_crown, add=TRUE)
  }
  Red <- subset(om.crown.mask,1)[,]
  Green <- subset(om.crown.mask,2)[,]
  Blue <- subset(om.crown.mask,3)[,]

  black_ind<-((Red+Green+Blue)/3)<80
  white_ind<-((Red+Green+Blue)/3)>200
  rgb_ind<-!(black_ind | white_ind)
  pBlack<-sum(black_ind, na.rm=TRUE)/length(Red)
  pWhite<-sum(white_ind, na.rm=TRUE)/length(Red)

  mRed<-mean(Red[rgb_ind], na.rm=T)
  mGreen<-mean(Green[rgb_ind], na.rm=T)
  mBlue<-mean(Blue[rgb_ind], na.rm=T)
  sdRed<-sd(Red[rgb_ind], na.rm=T)
  sdGreen<-sd(Green[rgb_ind], na.rm=T)
  sdBlue<-sd(Blue[rgb_ind], na.rm=T)
# Consider using image filters/convolutions and clustering here.
  return(data.frame(tree.id=i, mRed, mGreen, mBlue, sdRed, sdGreen, sdBlue, pBlack, pWhite))
}


source("R/imgPoly_match.R")
om.path<-'/Users/Tom/Documents/Work/RSPB/HRF/Drone/Missions/16 Block B/OM 4 cm cc blocked/OM 4 cm cc blocked'
chm.path<-'/Users/Tom/Documents/Work/RSPB/HRF/Drone/Missions/16 Block B/DSM 2'
chunk<-2
om.pattern=paste('Chunk', chunk, 'OM', '.*.tif$')
chm.pattern=paste('Chunk', chunk, 'CHM', '.*.tif$')

  crowns_om<-imgPoly_match2(om.path, pattern=om.pattern, crowns)
  crowns_chm<-imgPoly_match2(chm.path, pattern=chm.pattern, crowns)

  # You need to write a small checking function to make sure that the CHM only returns one match
  # ... or a way to handle multiple cases... perhaps choose the prefered chunk
  # For 1 tree:
  i<-30
  crown<-crowns[i,]
  crown_om<-brick(crowns_om[[i]], package='raster')
  crown_chm<-raster::raster(crowns_chm[[i]])
  crown_FD(crown_om, crown_chm, crown, plot='on')
  if(length(om.tile)==1)
    return(data.frame(tree.id=i, mRed, mGreen, mBlue, sdRed, sdGreen, sdBlue, pBlack, pWhite))
  else
    return(data.frame(tree.id=i, mRed=NA, mGreen=NA, mBlue=NA, sdRed=NA, sdGreen=NA, sdBlue=NA, pBlack=NA, pWhite=NA))
})
sp.crowns.features<-do.call(rbind, sp.crowns.features)

# Then you need to write a function which will allow the user to specify classifcations:

user.class<-select.list(classes, title='Choose the species')
