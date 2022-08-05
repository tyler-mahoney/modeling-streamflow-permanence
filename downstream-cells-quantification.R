library(sp)
library(raster)
library(dynatopmodel)
library(topmodel)
library(ecbtools) # remotes::install_github("ozjimbob/ecbtools")
library(cluster)
library(factoextra)
library(tidyverse)
library(NbClust)
library(basicClEval) # remotes::install_github("sarafrr/basicClEval")
library(ggplot2)
library(mapview)

#FR.HRU <- raster('SpatialInputData/FallingRockHRUs.tif')
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




## Run Kmeans
set.seed(123)

Fr.elev.bands <- cut(FR.dem, 4)
TWI.dinf.cut <- cut(TWI.dinf,4)
image(TWI.dinf.cut)
Fr.stack <- stack(TWI.dinf.cut,FR.soils, Fr.elev.bands) 

Fr.kmeans.2 <- raster.kmeans(Fr.stack, k = 20, iter.max=10, nstart = 10, geo = T, geo.weight = 1)
#writeRaster(Fr.kmeans.2,'C:/Users/david/OneDrive/Documents/ArcGIS/Projects/FallingRock/HRUs_21JULY22.tif',overwrite=T)
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

# This is where it would be good to reorder the HRUs. 

mean.twi <- data.frame(matrix(nrow=length(unique(Fr.kmeans.2)),ncol=5))
names(mean.twi) <- c('HRU.name','atb.bar.unsort','id.new','atb.bar','id')
mean.twi$HRU.name <- unique(Fr.kmeans.2)


Fr.hillslope.hru.matrix <- as.matrix(Fr.kmeans.2)
TWI.dinf.matrix <- as.matrix(TWI.dinf)

for (i in 1:length(unique(Fr.kmeans.2))) {
  
  mean.twi$atb.bar.unsort[i] <- mean(TWI.dinf.matrix[which(Fr.hillslope.hru.matrix==unique(Fr.kmeans.2)[i])])
}

mean.twi$atb.bar <- sort(mean.twi$atb.bar.unsort,decreasing=F)
mean.twi$id.new <- match(mean.twi$atb.bar.unsort,mean.twi$atb.bar)
mean.twi$id <- 1:length(unique(Fr.kmeans.2))

reclassify.matrix <- data.frame(matrix(nrow=length(unique(Fr.kmeans.2)),ncol=2))
names(reclassify.matrix) <- c('is','becomes')
reclassify.matrix$is <- mean.twi$HRU.name
reclassify.matrix$becomes <- mean.twi$id.new


Fr.kmeans.resample <- reclassify(Fr.kmeans.2,reclassify.matrix)

## Overlay the stream network atop the FR Hru raster
Fr.stream <- raster('SpatialInputData/FR_stream_network.tif')
Fr.stream[Fr.stream==0] <- NA
Fr.HRU <- Fr.kmeans.resample*10 + max(max(unique(Fr.stream)),100)
Fr.HRU[Fr.stream>0] <- Fr.stream[Fr.stream>0]
mapview(Fr.HRU)


#writeRaster(Fr.HRU,'C:/Users/david/OneDrive/Documents/ArcGIS/Projects/FallingRock/HRUs_26JULY22.tif',overwrite=T)

HRU.matrix <- as.matrix(Fr.HRU)

fdr.downstream.matrix <- matrix(nrow=nrow(HRU.matrix),ncol=ncol(HRU.matrix))
FR.fdr <- raster('SpatialInputData/fr_fdr')
fdr.matrix <- as.matrix(FR.fdr)

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

try.downstream.weighting.matrix <- matrix(nrow=length(unique(Fr.HRU)),ncol=length(unique(Fr.HRU)))
(start <- Sys.time())
for (i in 1:length(unique(Fr.HRU))) {
  print(unique(Fr.HRU)[i])
  for (j in 1:length(unique(Fr.HRU))) {
    try.downstream.weighting.matrix[i,j] <- length(which((HRU.matrix==unique(Fr.HRU)[i])&(fdr.downstream.matrix==unique(Fr.HRU)[j])))
  }
}
(end <- Sys.time())

options("scipen"=100, "digits"=5)

rownames(try.downstream.weighting.matrix) <- unique(Fr.HRU)
colnames(try.downstream.weighting.matrix) <- unique(Fr.HRU)

weighting.matrix <- try.downstream.weighting.matrix/rowSums(try.downstream.weighting.matrix)
weighting.matrix.df <- data.frame(weighting.matrix[,-1])
weighting.matrix.df$sum <- rowSums(weighting.matrix.df[1:53])

library(dynatopmodel)

# Calculate the mean TWI for each HRU - relabel the HRUs such that high valued TWI is the highest number HRU and low valued TWI is the lowest HRU

mean.twi <- data.frame(ncol=1,nrow=length(unique(Fr.HRU)))

for (i in 1:length(unique(Fr.HRU))) {
  mean.twi <- mean(TWI.dinf[Fr.HRU==unique(Fr.HRU)[i]])
}


# Create the groups matrix
cell.size <- as.numeric(readline('Enter the resolution of the raster in m (for 10-m raster enter 10): '))

area_cell <- function(raster.matrix,raster_ID) {
  area_c <- sum(raster.matrix==raster_ID, na.rm=T)
  return(area_c)
}

mean_raster <- function(raster.matrix, HRU.matrix, raster_ID) {
  mean_rast <- mean(raster.matrix[HRU.matrix==raster_ID],na.rm=T)
  return(mean_rast)
}

Fr.slope <- terrain(FR.dem,opt='slope',units='tangent')



groups <- data.frame(matrix(nrow=length(unique(Fr.HRU)),ncol=26))
names(groups) <- c('id','tag','chan.no','order','area_pc','area','sbar','atb.bar','gauge.id','catch.id','srz_max','ln_t0',
                   'm','srz0','td','vchan','vof','k0','CD','sd_max','pe_fact','vof_fact','rain_fact','mann.n','S0','ex_max')
groups$id <- unique(Fr.HRU)
groups$tag <- unique(Fr.HRU)
groups$chan.no[1:length(unique(Fr.stream))] <- unique(Fr.stream)
groups$order <- 1:length(unique(Fr.HRU))
groups$area_pc <- sapply(X=groups$id,FUN=area_cell,raster.matrix=HRU.matrix)/sum(!is.na(HRU.matrix))*100
groups$area <- groups$area_pc*sum(!is.na(HRU.matrix))*cell.size^2/1000/1000/100
groups$sbar <- zonal(x=Fr.slope,z=Fr.HRU,fun='mean')[,2]
groups$atb.bar[(1+length(unique(Fr.stream))):length(unique(Fr.HRU))] <- zonal(x=TWI.dinf,z=Fr.HRU,fun='mean')[(1+length(unique(Fr.stream))):length(unique(Fr.HRU)),2]
groups$gauge.id <- 1
groups$catch.id <- 1
groups$srz_max <- 0.1
groups$ln_t0 <- 7
groups$m <- 0.01
groups$srz0 <- 0
groups$td <- 1
groups$vchan[1:length(unique(Fr.stream))] <- 1000
groups$vof[(1+length(unique(Fr.stream))):length(unique(Fr.HRU))] <- 100
groups$k0 <- 1e+8
groups$CD <- 0.1
groups$sd_max <- 0.5
groups$pe_fact <- 1
groups$vof_fact <- 1
groups$rain_fact <- 1
groups$mann.n <- 0.01
groups$S0 <- 0.1
groups$ex_max <- 1












bromp.hru <- brompton$disc$hru
bromp.hru[brompton$chans$chans] <- 1
bromp.hru.matrix <- as.matrix(bromp.hru)
bromp.slope <- terrain(brompton$dem,opt='slope',units='tangent')
bromp.slope.matrix <- as.matrix(bromp.slope)
bromp.layers <- build_layers(brompton$dem)
bromp.twi.matrix <- as.matrix(bromp.layers$atb)

mean.slope.bromp <- matrix(nrow=length(unique(bromp.hru)),ncol=1)
for (i in 1:length(unique(bromp.hru))) {
  mean.slope.bromp[i] <- mean(bromp.slope.matrix[bromp.hru.matrix==unique(bromp.hru)[i]],na.rm=T)
}

try.zonal <- zonal(x=bromp.layers$atb,z=bromp.hru,fun="mean")

mean.twi.bromp <- matrix(nrow=length(unique(bromp.hru)),ncol=1)
for (i in 1:length(unique(bromp.hru))) {
  mean.twi.bromp[i] <- mean(bromp.twi.matrix[bromp.hru.matrix==unique(bromp.hru)[i]],na.rm=T)
}
