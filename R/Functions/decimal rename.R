#-----------------------------------------------
# 
# A function that renames shape files when the .5 has dropped off the end...
#
# I know that this is a bit of a hack but wahey!
#
# Tom Swinfield
# 17-02-16
#
#-----------------------------------------------

setwd('/Users/Tom/Documents/Work/R projects/UAVforestR/data/ITC trees_params/Matched')
all<-list.files()

shx<-grep('sobel_[0-9].shx', all)
shp<-grep('sobel_[0-9].shp', all)

all[shx]<-gsub('.shx', '.5.shx', all[shx])
all[shp]<-gsub('.shp', '.5.shp', all[shp])

file.rename(from=list.files(), to=all)

setwd('/Users/Tom/Documents/Work/R projects/UAVforestR')
