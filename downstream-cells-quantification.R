library(sp)
library(raster)
library(dynatopmodel)
library(topmodel)
library(ecbtools)
library(cluster)
library(factoextra)
library(tidyverse)
library(NbClust)
library(basicClEval)
library(ggplot2)

FR.HRU <- raster('SpatialInputData/FallingRockHRUs.tif')
FR.fdr <- raster('SpatialInputData/fr_fdr')
FR.dem <- raster('SpatialInputData/fr1meterDEM.tif')
FR.soils <- raster('SpatialInputData/FRSoils1mMask.tif')
FR.fac.dinf <- raster('SpatialInputData/fr_fac_dinf.tif')
TWI.dinf <- raster('C:/Users/david/OneDrive/Documents/ArcGIS/Projects/FallingRock/twi_dinf.tif')
TWI.dinf[TWI.dinf==0] <- NA

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

Fr.elev.bands <- cut(FR.dem, 4)
TWI.dinf.cut <- cut(TWI.dinf,4)
image(TWI.dinf.cut)
Fr.stack <- stack(TWI.dinf.cut,FR.soils, Fr.elev.bands) 

Fr.kmeans.2 <- raster.kmeans(Fr.stack, k = 20, iter.max=10, nstart = 10, geo = T, geo.weight = 1)
writeRaster(Fr.kmeans.2,'C:/Users/david/OneDrive/Documents/ArcGIS/Projects/FallingRock/HRUs_21JULY22.tif',overwrite=T)
x11()
image(Fr.kmeans.2)

Fr.kmeans.hist <- hist(as.vector(Fr.kmeans.2))
x11()
plot(Fr.kmeans.hist)
set.seed(123)
Fr.kmeans.df <- data.frame('elev.bands'=as.vector(Fr.elev.bands),'TWI'=as.vector(TWI.dinf.cut),'soils'=as.vector(FR.soils))
Fr.kmeans.df <- Fr.kmeans.df[complete.cases(Fr.kmeans.df),]
Fr.kmeans <- kmeans(Fr.kmeans.df,centers = 10, iter.max = 10, nstart = 10)
Fr.kmeans.wcss <- wcss(Fr.kmeans.df,Fr.kmeans$cluster)

x11()
hist(Fr.kmeans$cluster)

wss <- function(k) {
  kmeans(Fr.kmeans.df,k,nstart=10)$tot.withinss
}

# Compute and plot wss for k = 1 to k = 40
k.values <- 1:40

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values,wss)
x11()

plot(k.values,wss_values,type='b',pch=19,frame=F,xlab='Number of clusters',ylab='Total within-cluster sum of sq')

Fr.HRU.vector <- data.frame('hru' = as.vector(Fr.kmeans.2))
Fr.TWI.vector <- data.frame('TWI' = as.vector(TWI.dinf.cut))
Fr.TWI.vector.nc <- data.frame('TWI.nc' = as.vector(TWI.dinf))
Fr.soils.vector <- data.frame('soils' = as.vector(FR.soils))
Fr.elev.bands.vector <- data.frame('Elev.bands' = as.vector(Fr.elev.bands))
Fr.elev.bands.vector.nc <- data.frame('Elev.bands.nc'=as.vector(FR.dem))
Fr.rasters.df <- data.frame(Fr.HRU.vector,Fr.TWI.vector,Fr.soils.vector,Fr.elev.bands.vector,Fr.elev.bands.vector.nc, Fr.TWI.vector.nc)
Fr.rasters.df.clean <- Fr.rasters.df[complete.cases(Fr.rasters.df),]

## Get the total ss within each cluster

wss.raster <- function(k) {
  Fr.kmeans.raster <- raster.kmeans(Fr.stack, k, iter.max=10, nstart = 10, geo = T, geo.weight = 1)
  Fr.HRU.vector <- data.frame('hru' = as.vector(Fr.kmeans.raster))
  Fr.TWI.vector <- data.frame('TWI' = as.vector(TWI.dinf.cut))
  Fr.soils.vector <- data.frame('soils' = as.vector(FR.soils))
  Fr.elev.bands.vector <- data.frame('Elev.bands' = as.vector(Fr.elev.bands))
  Fr.rasters.df <- data.frame(Fr.HRU.vector,Fr.TWI.vector,Fr.soils.vector,Fr.elev.bands.vector)
  Fr.rasters.df.clean <- Fr.rasters.df[complete.cases(Fr.rasters.df),]
  within.ss <- wcss(Fr.rasters.df.clean[,-1],Fr.rasters.df.clean[,1])
  return(within.ss$WCSS)
}

set.seed(123)
# Compute and plot wss for k = 1 to k = 40
k.values <- 1:40

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values,wss.raster)
x11()

plot(k.values,wss_values,type='b',pch=19,frame=F,xlab='Number of clusters',ylab='Total within-cluster sum of sq')

fviz_cluster(Fr.rasters.df.clean[,1],Fr.rasters.df.clean[,-1])


library(ggplot2)
x11()
ggplot(data=Fr.rasters.df.clean, aes(x=TWI,y=Elev.bands,col=hru))+geom_point()

## Overlay the stream network atop the FR Hru raster
Fr.stream <- raster('SpatialInputData/FR_stream_network.tif')
Fr.HRU <- Fr.kmeans.2
Fr.HRU.ov <- overlay(Fr.stream,Fr.HRU,fun=sum)
