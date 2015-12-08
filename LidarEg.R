library(raster)
#library(proj4)
source("readLAS.R")

liddat <- readLAS("2013_SJER_AOP_point_cloud_unclassified.las")
liddat <- data.frame(liddat)
resol <- 10 # meters
ras <- raster(ncols=(ceiling(max(liddat$x))-floor(min(liddat$x)))/2,
              nrows=(ceiling(max(liddat$y))-floor(min(liddat$y)))/2,
              xmn=floor(min(liddat$x)),
              xmx=ceiling(max(liddat$x)),
              ymn=floor(min(liddat$y)),
              ymx=ceiling(max(liddat$y)),
              resolution=resol)
sp1 <- SpatialPoints(liddat[,c("x","y")])
r1 <- rasterize(sp1, ras, liddat$z, fun=min)
plot(r1)