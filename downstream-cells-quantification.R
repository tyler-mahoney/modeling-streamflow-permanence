library(sp)
library(raster)
library(dynatopmodel)
library(topmodel)


FR.HRU <- raster('SpatialInputData/FallingRockHRUs.tif')
FR.fdr <- raster('SpatialInputData/fr_fdr')
FR.dem <- raster('SpatialInputData/fr1meterDEM.tif')
FR.soils <- raster('SpatialInputData/FRSoils1mMask.tif')
FR.fac.dinf <- raster('SpatialInputData/fr_fac_dinf.tif')
TWI.dinf <- raster('C:/Users/david/OneDrive/Documents/ArcGIS/Projects/FallingRock/twi_dinf.tif')

HRU.matrix <- as.matrix(FR.HRU)
fdr.matrix <- as.matrix(FR.fdr)
dem.matrix <- as.matrix(FR.dem)
soils.matrix <- as.matrix(FR.soils)
TWI.dinf.matrix <- as.matrix(TWI.dinf)
TWI.dinf.matrix[TWI.dinf.matrix==0] <- NA

twi <- build_layers(FR.dem)
twi <- twi$atb

normalize.matrix <- function(in.matrix) {
  out.matrix <- (in.matrix-min(in.matrix,na.rm=T))/(max(in.matrix,na.rm=T)-min(in.matrix,na.rm=T))
  return(out.matrix)
}

twi.dist <- normalize.matrix(as.matrix(twi$atb))
twi.hist <- hist(twi.dist,plot = F)

twi.dist.2 <- as.matrix(FR.TWI.calc)
twi.dist.2[twi.dist.2 == 0] <- NA
twi.dist.2 <- normalize.matrix(twi.dist.2)
twi.hist.2 <- hist(twi.dist.2, plot = F)

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")


plot(twi.hist,col=c3)
plot(twi.hist.dinf,col=c1,add=T)



fdr.downstream.matrix <- matrix(nrow=nrow(HRU.matrix),ncol=ncol(HRU.matrix))
 

for (i in 1:nrow(HRU.matrix)) {
  for (j in 1:ncol(HRU.matrix)) {
    if (is.na(fdr.matrix[i,j])) {
      fdr.downstream.matrix[i,j] <- NA
    } else {
      if (fdr.matrix[i,j]==1) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i,j+1]
      } else if (fdr.matrix[i,j]==2) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j+1]
      } else if (fdr.matrix[i,j]==4) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j]
      } else if (fdr.matrix[i,j]==8) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j-1]
      } else if (fdr.matrix[i,j]==16) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i,j-1]
      } else if (fdr.matrix[i,j]==32) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j-1]
      } else if (fdr.matrix[i,j]==64) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j]
      } else if (fdr.matrix[i,j]==128) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j+1]
      } else {
        fdr.downstream.matrix[i,j] <- NA
      }
    }
  }
}

try.downstream.weighting.matrix <- matrix(nrow=length(unique(FR.HRU)),ncol=length(unique(FR.HRU)))
(start <- Sys.time())
for (i in 1:length(unique(FR.HRU))) {
  print(unique(FR.HRU)[i])
  for (j in 1:length(unique(FR.HRU))) {
    try.downstream.weighting.matrix[i,j] <- length(which((HRU.matrix==unique(FR.HRU)[i])&(fdr.downstream.matrix==unique(FR.HRU)[j])))
  }
}
(end <- Sys.time())

options("scipen"=100, "digits"=5)

rownames(downstream.weighting.matrix) <- unique(FR.HRU)
colnames(downstream.weighting.matrix) <- unique(FR.HRU)

weighting.matrix <- downstream.weighting.matrix/rowSums(downstream.weighting.matrix)

## Run Kmeans
set.seed(123)
plot(x=TWI.dinf.matrix,y=soils.matrix)

Fr.kmeans.df <- data.frame('TWI'=as.vector(TWI.dinf.matrix),'soils'= as.vector(soils.matrix))
Fr.kmeans <- kmeans(Fr.kmeans.df,centers = 20, iter.max = 5, nstart = 10)
