# Ocean productivity data from Oregon State
# www.science.oregonstate.edu/ocean.productivity
# Downloaded monthly NPP rasters from the updated CbPM model for 2011
# on 17 July 2013

setwd('//bioark.bio.unc.edu/hurlbertlab/gis/oceanproductivity')
analysis_dir = '//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/sebastes/'
library(R.utils)

files = list.files(pattern='\\.gz$')

for (f in files) {
  gunzip(f)
}

hdffiles = list.files(pattern='\\.hdf$')

# This uses a java script available from http://asascience.com/software/downloads/
# to convert hdf viles to netcdf, however it does not work when the script
# is on my local machine while the files are on Bioark, so I actually
# converted them manually in the shell rather than with this loop
for (h in hdffiles) {
  shell(paste("java -jar modis_hdf2nc.jar ", h, sep=""))
}

# Moved .nc files back to Bioark



#############
library(ncdf)

ncfiles = list.files(pattern='\\.nc$')
for (m in 1:12) {
  assign(paste("ncdf",m,sep=""), open.ncdf(ncfiles[m]))
  getvar = get.var.ncdf(get(paste("ncdf",m,sep="")))
  temp = t(raster(getvar, xmn = -90, ymn = -180, xmx = 90, ymx = 180))
  temp[temp == -9999] <- NA
  assign(paste("npp", m, sep=""), temp)
}

mean.npp = mosaic(npp1, npp2, npp3, npp4, npp5, npp6, npp7, npp8, npp9, npp10, npp11, npp12, fun = mean)
plot(mean.npp)

# Aggregate down to 1 degree resolution, and crop to NorthEastern Pacific
mean.npp.1dg = aggregate(mean.npp, fact = 6, fun = mean)
xmin = -152
xmax = -105
ymin = 22
ymax = 60
NEP = extent(xmin, xmax, ymin, ymax)
NPP.NEP = crop(mean.npp.1dg, NEP)  
plot(NPP.NEP)  

# Pull out eastern most values at each latitude within cropped raster
NPP.NEPmat = as.matrix(NPP.NEP)
coast.npp = unlist(apply(NPP.NEPmat, 1, function(x) x[min(which(is.na(x)))-1]))
# Because Alaska coastline goes south again as you go west, the above command
# does not provide a value for 60 degrees latitude. Add it manually.
coast.npp = c(NPP.NEPmat[1, max(which(!is.na(NPP.NEplotPmat[1,])))], coast.npp)
lat = ymax:(ymin+1)
plot(lat, coast.npp)

res = 2 # resolution in degrees
mean.npp.2d = aggregate(mean.npp, fact = 6*res, fun = mean)
NPP2.NEP = crop(mean.npp.2d, NEP)
coast.npp2d = unlist(apply(as.matrix(NPP2.NEP), 1, function(x) x[min(which(is.na(x)))-1]))
lat2d = seq(ymax, ymin+res, -res)
lat2d. = lat2d[2:length(lat2d)] #no NPP value for first latitude row
plot(lat2d., coast.npp2d, type = 'l')


npp.1deg = data.frame(latitude=lat, npp=coast.npp)
npp.2deg = data.frame(latitude=lat2d., npp=coast.npp2d)
write.csv(npp.1deg,paste(analysis_dir,'/PacificNPPvsLatitude_1dg.csv',sep=''), row.names=F)
write.csv(npp.2deg,paste(analysis_dir,'/PacificNPPvsLatitude_2dg.csv',sep=''), row.names=F)

# Moving average function over n=5 latitudinal bins (+2, + 1, 0, -1, -2)
ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
n = 5
npp.1deg.ma = ma(npp.1deg$npp, n)
# fill in MA values for endpoints
npp.1deg.ma[1] = mean(npp.1deg$npp[1:floor(n/2)])
npp.1deg.ma[2] = mean(npp.1deg$npp[1:4])
npp.1deg.ma[length(npp.1deg.ma)] = mean(npp.1deg$npp[(length(npp.1deg.ma)-2):length(npp.1deg.ma)])
npp.1deg.ma[length(npp.1deg.ma)-1] = mean(npp.1deg$npp[(length(npp.1deg.ma)-3):length(npp.1deg.ma)])

npp.1d.ma = data.frame(latitude = npp.1deg$lat, npp = npp.1deg.ma)
write.csv(npp.1d.ma, paste(analysis_dir,'/PacificNPPvsLatitude_1dg_MA.csv',sep = ''), row.names=F)
plot(npp.1d.ma$latitude, npp.1d.ma$npp, type='l', col='green')
