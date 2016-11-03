# A function to find the minimum extent of two rasters
# 161007

min_extent<-function(r1, r2)
{
  ex1<-extent(r1) 
  ex2<-extent(r2)
  ex<-ex1
  ex@xmin<-max(ex1@xmin, ex2@xmin)
  ex@xmax<-min(ex1@xmax, ex2@xmax)
  ex@ymin<-max(ex1@ymin, ex2@ymin)
  ex@ymax<-min(ex1@ymax, ex2@ymax)
  return(ex)
}