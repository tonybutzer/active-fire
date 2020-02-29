l8fire_detect.R <-
function(b7,b6,b5,b4,b3,b2,fname="out.tiff")
{
require(raster)
require(rgdal)
dyn.load("l8fire_detection.dll")
out <- band7
out[,]<- 0
jcol<- dim(band7)[2]
irow<- dim(band7)[1]
bs<- .C("l8fire_detection",jcol=as.integer(jcol),irow=as.integer(irow),data7=as.double(b7[,]),data6=as.double(b6[,]),data5=as.double(b5[,]),data4=as.double(b4[,]),data3=as.double(b3[,]),data2=as.double(b2[,]),out=as.integer(out[,]))
dyn.unload("l8fire_detection.dll")
out[,]<-bs$out
print(unique(bs$out))
print(fname)
writeRaster(out,fname,format="GTiff",overwrite=T,datatype="INT1U")
}
