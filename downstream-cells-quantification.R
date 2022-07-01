library(sp)
library(raster)

FR.HRU <- raster('SpatialInputData/FallingRockHRUs.tif')
FR.fdr <- raster('SpatialInputData/fr_fdr')

HRU.matrix <- as.matrix(FR.HRU)
fdr.matrix <- as.matrix(FR.fdr)

fdr.downstream.matrix <- matrix(nrow=nrow(HRU.matrix),ncol=ncol(HRU.matrix))



for (i in 1:nrow(HRU.matrix)) {
  for (j in 1:ncol(HRU.matrix)) {
    if (is.na(fdr.matrix[i,j])) {
      fdr.downstream.matrix[i,j] <- NA
    } else {
      if (fdr.matrix[i,j]==1) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i,j+1]
      } else if (fdr.matrix[i,j]==2) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j+1]
      } else if (fdr.matrix[i,j]==4) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j]
      } else if (fdr.matrix[i,j]==8) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j-1]
      } else if (fdr.matrix[i,j]==16) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i,j-1]
      } else if (fdr.matrix[i,j]==32) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j-1]
      } else if (fdr.matrix[i,j]==64) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j]
      } else if (fdr.matrix[i,j]==128) {
        fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j+1]
      } else {
        fdr.downstream.matrix[i,j] <- NA
      }
    }
  }
}


for (i in 1:20) {
  for (j in 1:20) {
    fdr.downstream.matrix[i,j] <- i+j
  }
}

fdr.matrix[1000,1000]==64


i=1
j=1

if (is.na(fdr.matrix[i,j])) {
  fdr.downstream.matrix[i,j] <- NA
} else {
  if (fdr.matrix[i,j]==1) {
    fdr.downstream.matrix[i,j] <- HRU.matrix[i,j+1]
  } else if (fdr.matrix[i,j]==2) {
    fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j+1]
  } else if (fdr.matrix[i,j]==4) {
    fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j]
  } else if (fdr.matrix[i,j]==8) {
    fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j-1]
  } else if (fdr.matrix[i,j]==16) {
    fdr.downstream.matrix[i,j] <- HRU.matrix[i,j-1]
  } else if (fdr.matrix[i,j]==32) {
    fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j-1]
  } else if (fdr.matrix[i,j]==64) {
    fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j]
  } else if (fdr.matrix[i,j]==128) {
    fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j+1]
  } else {
    fdr.downstream.matrix[i,j] <- NA
  }
}
  


if (fdr.matrix[i,j]==1) {
  fdr.downstream.matrix[i,j] <- HRU.matrix[i,j+1]
} else if (fdr.matrix[i,j]==2) {
  fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j+1]
} else if (fdr.matrix[i,j]==4) {
  fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j]
} else if (fdr.matrix[i,j]==8) {
  fdr.downstream.matrix[i,j] <- HRU.matrix[i-1,j-1]
} else if (fdr.matrix[i,j]==16) {
  fdr.downstream.matrix[i,j] <- HRU.matrix[i,j-1]
} else if (fdr.matrix[i,j]==32) {
  fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j-1]
} else if (fdr.matrix[i,j]==64) {
  fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j]
} else if (fdr.matrix[i,j]==128) {
  fdr.downstream.matrix[i,j] <- HRU.matrix[i+1,j+1]
} else {
  fdr.downstream.matrix[i,j] <- NA
}