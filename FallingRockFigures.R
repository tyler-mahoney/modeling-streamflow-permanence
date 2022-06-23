## Create figures for reaches
reach.figures <- function(routed.flow.instant= explicit.routing$routed.flow.instant.m3_hr, explicit.reach.table=explicit.reach.table,
                          routed.flow.instant.mm_s=explicit.routing$routed.flow.instant.mm_s) {
  library(fitdistrplus)
  library(hydroTSM)
  Q <- data.frame(explicit.routing$routed.flow.instant.mm_s)
  r.names <- names(Q)
  #r.name <- 'X3794'
  
  # Histograms of log flow
  x11()
  par(mfrow=c(3,11))
  par(mar=c(5,2,1,1))
  for (r.name in 1:length(r.names)) {
    Q.in <- log(as.numeric(Q[,r.name]))
    main.title <- paste0('Reach ', r.name, ' Order ' , explicit.reach.table[r.name,5])
    hist(Q.in,
         main=main.title,
         xlab ='ln(Q [m^3/hr])')
    
  }
  
  # weibull distribution of log flow
  for (r.name in 27:length(r.names)) {
    infile <- paste('C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/Figures/weibull_pdfs/',
                      'Weibull_reach',r.name,'.jpg',sep='')
    jpeg(file=infile)
    Q.in <- log(as.numeric((Q[,r.name])))
    if (any(Q.in<=0)) {
      Q.in <- na.omit(as.numeric(Q.in-min(Q.in)+.01))
      }
   # Q.in.2 <- na.omit(na.exclude(Q.in))
    
    main.title <- paste0('Reach ', r.name, ' Order ' , explicit.reach.table[r.name,5])
    fit.weibull <- fitdist(Q.in[Q.in>=0],distr='weibull', method='mle',lower=c(0,0)) # Weibull distributionby shifting data such that negative values are included...
    plot(fit.weibull)
    dev.off()
  }
  
  #exp distribution
  for (r.name in 1:length(r.names)) {
    infile <- paste('C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/FR_Spatial_Function_23FEB2021/Routing/pdfs/',
                    'reach',r.name,'.jpg',sep='')
    jpeg(file=infile)
    Q.in <- (as.numeric((Q[,r.name])))
    main.title <- paste0('Reach ', r.name, ' Order ' , explicit.reach.table[r.name,5])
    fit.exp <- fitdist(Q.in,'exp')
    plot(fit.exp)
    dev.off()
  }
  
  
  # exp distribution
  for (r.name in 1:length(r.names)) {
    infile <- paste('C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/FR_Spatial_Function_23FEB2021/Routing/pdfs/',
                    'reach',r.name,'.jpg',sep='')
    jpeg(file=infile)
    Q.in <- (as.numeric((Q[,r.name])))
    main.title <- paste0('Reach ', r.name, ' Order ' , explicit.reach.table[r.name,5])
    fit.exp <- fitdist(Q.in,'exp')
    plot(fit.exp)
    dev.off()
  }
  
  # ecdf plots
  x11()
  par(mfrow=c(3,11))
  par(mar=c(5,2,1,1))
  #ggplot(Q.plot, aes(Q)) + stat_ecdf(geom = 'step')
  for (r.name in 1:length(r.names)) {
    Q.in.2 <- ifelse(Q[,r.name]<=0, 0.001, Q[,r.name] ) 
    Q.in.2 <- log(as.numeric(Q.in.2))
    
    Q.in <- ecdf(as.numeric(Q.in.2))
    main.title <- paste0('Reach ', r.name, ' Order ' , explicit.reach.table[r.name,5])
    plot(Q.in,
         main=main.title,
         xlab ='ln(Q [m^3/hr])')
  }
  
  # FDC
  x11()
  par(mfrow=c(3,11))
  par(mar=c(5,2,1,1))
  sim.timestamp <- time(result.optim.run$run$Qsim[1:nrow(result.optim.run$run$fluxes$qbf)])
  #ggplot(Q.plot, aes(Q)) + stat_ecdf(geom = 'step')
  for (r.name in 1:length(r.names)) {
    Q.in <- zoo(Q[,r.name],sim.timestamp)
    flow.dur <- fdc(Q.in, lQ.thr=0.7, hQ.thr=0.2, plot=TRUE, log="y",
                    main = paste0("FDC reach ",r.name), xlab = "% time exceeded", ylab = 'Q [m3/hr]', yat=c(0.1, .1, 1),
                    xat=c(0.01, 0.025, .05), col=palette("default")[1:NCOL(Q.in)], pch=1:NCOL(Q.in), lwd=rep(1, NCOL(Q.in)), 
                    lty=1:NCOL(Q.in), cex=0.4, cex.axis=1.2, cex.lab=1.2, leg.txt=NULL,
                    leg.cex=1, leg.pos="topright", verbose=TRUE, thr.shw=TRUE, new=TRUE)
  }
  
  
  x11()
  # Combined ecdf
  for (r.name in 1:length(r.names)) {
    Q.in.2 <- ifelse(Q[,r.name]<=0, 0.01, Q[,r.name] ) 
    Q.in.2 <- log(as.numeric(Q.in.2))
    
    ecdf.in <- ecdf(as.numeric(Q.in.2))
    main.title <- paste0('Reach ', r.name, ' Order ' , explicit.reach.table[r.name,5])
    plot(ecdf.in,
         xlab ='ln(Q [m^3/hr])')
    par()
  }
  
  
  # Combined FDC
  x11()
  #par(mfrow=c(3,11))
  #par(mar=c(5,2,1,1))
  sim.timestamp <- time(result.optim.run$run$Qsim[1:nrow(result.optim.run$run$fluxes$qbf)])
  #ggplot(Q.plot, aes(Q)) + stat_ecdf(geom = 'step')
  for (r.name in 1:length(r.names)) {
    if (r.name==26) {
      r.name <- 25
    }
    Q.in <- zoo(Q[,r.name],sim.timestamp)
    order <- explicit.reach.table[r.name,5]
    if (order ==4){
      col='black'; lwd = 3
    } else if (order ==3) {
      col='blue'; lwd =2
    } else if (order == 2 ) {
      col ='green'; lwd =1;
    } else {
      col = 'red'; lwd =.25
    }
    flow.dur <- fdc.zoo(Q.in, lQ.thr=0.7, hQ.thr=0.2, plot=TRUE, log="y",
                    main = paste0("FDC reach ",r.name), xlab = "% time exceeded", ylab = 'Q [m3/hr]', yat=c(0.1, .1, 1),
                    xat=c(0.01, 0.025, .05), col=col, pch=NA, lwd=1, 
                    lty=1:NCOL(Q.in), cex=0.4, cex.axis=1.2, cex.lab=1.2, leg.txt=NULL,
                    leg.cex=1, leg.pos="topright", verbose=TRUE, thr.shw=TRUE, new=TRUE)
    par(new=T)
  }
  
  
  
  Flow.duration.reaches <- matrix(0,nrow = length(Q[,1]),ncol = length(r.names))
  Q.in.rank <- matrix(0,nrow = length(Q[,1]),ncol = length(r.names))
  colnames(Flow.duration.reaches) <- r.names
  colnames(Q.in.rank) <- r.names
  x11()
  
  ## Calculate FDC
  for (r.name in 1:length(r.names)) {
    if (r.name==26) {
      r.name <- 25
    }
    Q.in <- Q[,r.name]
    order <- explicit.reach.table[r.name,5]
    Q.in.sort <- data.frame(flow.asc=sort(Q.in,decreasing =T))
    rank <- 1:length(Q.in)
    Q.in.sort$rank <- rank
    Q.in.sort$Prob.Q.in <- 100*(Q.in.sort$rank/(length(Q.in)+1))
    Flow.duration.reaches[,paste0('X',as.character(explicit.reach.table$link_no[r.name]))] <- Q.in.sort$Prob.Q.in
    Q.in.rank[,paste0('X',as.character(explicit.reach.table$link_no[r.name]))] <- Q.in.sort$flow.asc
    order <- explicit.reach.table[r.name,5]
    if (order ==4){
      col='black'; lwd = 4
    } else if (order ==3) {
      col='#E69F00'; lwd =3
    } else if (order == 2 ) {
      col ='#0072B2'; lwd =2;
    } else {
      col = '#CC79A7'; lwd =.1
    }
    
    plot(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
         log='y',ylab='Q [mm/s]',xlab='% Exceedance',col=col,lwd=lwd)
    
    par(new=T)
  }
  # PLot this such that the flows are all on the same axis
  x11()
  # Plot the FDCs
  for (r.name in 1:length(r.names)) {
    order <- explicit.reach.table[r.name,5]
    if (order ==4){
      col='black'; lwd = 3
    } else if (order ==3) {
      col='#0072B2'; lwd = 2
    } else if (order == 2 ) {
      col ='#E69F00'; lwd =1;
    } else {
      col = '#CC79A7'; lwd =.1
    }
    
    if (r.name==1) { 
      #plot(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
      #     log='y',ylab='Q [m3/hr]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-03,1e05))
      plot(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
           log='y',ylab='Q [mm/s]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-10,1e-1),
           cex =1.5,cex.axis=1.5,cex.lab=1.5)
      
    } else {
      #lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
       #    ylab='Q [m3/hr]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-03,1e05))
      lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
            ylab='Q [mm/s]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-10,1e-1))
    }
    
  }
  
  legend(65,max(Q.in.rank)*0.5, legend = c('Order 4','Order 3','Order 2','Order 1'),
         col=c('black','#0072B2','#E69F00','#CC79A7'),lwd=c(4,3,2,.1), title='Stream Order',
         cex =1.5)
  
  # Normalize the flows and plot this to see differences in the FDCs
  normalize.Q <- function(x) {(x-min(x))/(max(x)-min(x))}
  flow.norm <- apply(Q.in.rank,2,normalize.Q)
  
  flow.norm <- flow.norm[1:(nrow(flow.norm)-1),]
  flow.dur.norm <- Flow.duration.reaches[1:(nrow(Flow.duration.reaches)-1),]
  # PLot this such that the flows are all on the same axis
  x11()
  # Plot the FDCs
  for (r.name in 1:length(r.names)) {
    order <- explicit.reach.table[r.name,5]
    if (order ==4){
      col='black'; lwd = 3
    } else if (order ==3) {
      col='#E69F00'; lwd = 2
    } else if (order == 2 ) {
      col ='#0072B2'; lwd =1;
    } else {
      col = '#CC79A7'; lwd =.1
    }
    if(r.name!=26) {
    if (r.name==1) { 
      #plot(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
      #     log='y',ylab='Q [m3/hr]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-03,1e05))
      plot(x=as.numeric(flow.dur.norm[,r.name]),y=as.numeric(flow.norm[,r.name]), type ='l',
          log ='y', ylab='Q [normalized]',xlab='% Exceedance',col=col,lwd=lwd)
      
    } else {
      #lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
      #     ylab='Q [m3/hr]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-03,1e05))
      lines(x=as.numeric(flow.dur.norm[,r.name]),y=as.numeric(flow.norm[,r.name]), type ='l',
            ylab='Q [normalized]',xlab='% Exceedance',col=col,lwd=lwd)
    }
    }
  }
  
  legend(70,max(Q.in.rank)*0.5, legend = c('Order 4','Order 3','Order 2','Order 1'),
         col=c('black','#E69F00','#0072B2','#CC79A7'),lwd=c(4,3,2,.1), title='Stream Order')
  
  # Plot the average of each stream order 
  stream.orders <- explicit.reach.table[,5]
  (which(stream.orders==1))
  
  order.1 <- (Q.in.rank[,which(stream.orders==1)])
  order.1 <- rowMeans(order.1[,-10])
  order.2 <- rowMeans(Q.in.rank[,which(stream.orders==2)])
  order.3 <- rowMeans(Q.in.rank[,which(stream.orders==3)])
  order.4 <- rowMeans(Q.in.rank[,which(stream.orders==4)])
  
  
  x11()
  plot(x=as.numeric(Flow.duration.reaches[,1]),y=as.numeric(order.4), type ='l',
       log='y',ylab='Q [mm/s]',xlab='% Exceedance',col='black',lwd=2,ylim=c(1e-10,1e-1))
  
  lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(order.3), type ='l',
        ylab='Q [mm/s]',xlab='% Exceedance',col='#E69F00',lwd=2,ylim=c(1e-10,1e-1))
  
  lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(order.2), type ='l',
        ylab='Q [mm/s]',xlab='% Exceedance',col='#0072B2',lwd=2,ylim=c(1e-10,1e-1))
  
  lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(order.1), type ='l',
        ylab='Q [mm/s]',xlab='% Exceedance',col='#CC79A7',lwd=2,ylim=c(1e-10,1e-1))
  
  legend(70,max(Q.in.rank)*0.5, legend = c('Order 4','Order 3','Order 2','Order 1'),
         col=c('black','#E69F00','#0072B2','#CC79A7'),lwd=c(4,3,2,1), title='Stream Order') 
    
  
  # Plot fdc for explicit reaches
  
  
  example.reaches <- c('X3802','X3474','X2738','X58')
  
  # Plot the FDCs
  
  for (r.name in example.reaches) {
    x11()
    par(mar=c(4,6,4,4))
    name.shorter <- substr(r.name,2,nchar(r.name))
    order <- explicit.reach.table[which(explicit.reach.table$link_no == name.shorter),5]
    if (order ==4){
      col='black'; lwd = 5
    } else if (order ==3) {
      col='#0072B2'; lwd = 4
    } else if (order == 2 ) {
      col ='#E69F00'; lwd =3;
    } else {
      col = '#CC79A7'; lwd =3
    }
    
    main.title.FDC <- paste0('Reach ', name.shorter, ' Flow Exceedance')
    
      #plot(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
      #     log='y',ylab='Q [m3/hr]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-03,1e05))
      plot(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
           log='y',ylab='Q [mm/s]',main=main.title.FDC, col=col,lwd=lwd,xlab='',
           ylim=c(1e-10,1e-1),cex.lab=2, cex.axis=2,cex=2, cex.main=2)
    
    
  }

 
 
  
  x11()
  # Combined ecdf
  for (r.name in 1:length(r.names)) {
    Q.in.2 <- ifelse(Q[,r.name]<=0, 0.000001, Q[,r.name] ) 
    Q.in.2 <- log(as.numeric(Q.in.2))
    order <- explicit.reach.table[r.name,5]
    if (order ==4){
      col='black'; lwd = 3
    } else if (order ==3) {
      col='#E69F00'; lwd = 2
    } else if (order == 2 ) {
      col ='#0072B2'; lwd =1;
    } else {
      col = '#CC79A7'; lwd =.1
    }
    ecdf.in <- ecdf(as.numeric(Q.in.2))
    main.title <- paste0('Reach ', r.name, ' Order ' , explicit.reach.table[r.name,5])
    if (r.name!=26) {
    if(r.name ==1) {
      plot(ecdf.in,
           xlab ='ln(Q [mm/s])', col=col,lwd=lwd,xlim=c(-23,-5), main='Empirical CDF',
           ylab='Frequency')
    } else {
      lines(ecdf.in, col = col,lwd=lwd,)
    }
    }
  }
  
  legend(-22,.8, legend = c('Order 4','Order 3','Order 2','Order 1'),
         col=c('black','#E69F00','#0072B2','#CC79A7'),lwd=c(4,3,2,1), title='Stream Order') 
  
  
  
  
  ## Read in the 1/0 timeseries
  logger.dir <- 'C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/3 DATA/2-PROJECT1-ROBINSON-FOREST/4-LOGGERDATA'
  files <- list.files(logger.dir,pattern='*clean.csv')
  files <- files[1:4]
  # plot the logger data
  x11()
  par(mfrow=c(4,1))
  par(mar=c(2,5,1,1))
  for (file.name in files) {
    short.name <- substr(file.name,start=1,stop=(nchar(file.name)-10))
    TS.logger <- read.csv(paste0(logger.dir,'/',file.name))
    #TS.logger.clean <- TS.logger %>% select(RoundTime,ApproxState)
    TS.logger.clean <- data.frame(Date=TS.logger$RoundTime,State=TS.logger$ApproxState)
    
    TS.logger.clean$State[which(TS.logger.clean$State=='Missing data')] <- 3
    TS.logger.clean$Date <- (ymd_hms(TS.logger.clean$Date))
    #sum(is.finite(TS.logger.clean$State))
    TS.logger.clean$Date2 <- as.Date(TS.logger.clean$Date)
    TS.logger.clean$State <- as.numeric(TS.logger.clean$State)
    
    TS.zoo <- zoo(TS.logger.clean$State,TS.logger.clean$Date)
    #TS.logger.clean <- data.frame(TS.logger.clean)
    #time(TS.logger.clean) <- time(TS.logger)
    
    #plot(TS.logger.clean$State)
    #image(TS.logger.clean$ApproxState, col=c('red','blue','yellow'))
    ones <- rep(1,nrow(TS.logger.clean))
    colors = c("black", "#0072B2",'white')
    #on <- sum(TS.logger.clean$State==1)
    #off <- sum(TS.logger.clean$State==0)
    #Total.logger <- on+off
    #percent.on <- round(on/Total.logger*100)
    plot(TS.logger.clean$Date2,ones, type="h", col=colors[TS.zoo +1], ylim=c(0,1), 
         xaxt='n',yaxt='n',lwd=2,ylab=paste0(short.name), cex=2, cex.lab=1.5,cex.axis=1.5)
    ytick <- c(0,1)
    axis(side=2,at=ytick)
    axis.Date(1,at=TS.logger.clean$Date2,format = '%m-%y',tick=F,cex.axis=1.5)
  }
  
  legend(1,10, legend = c('Dry','Wet','No data'),
         col=c('black','#0072B2','white'),lwd=c(2,2,2), title='Legend') 
  
  
  
  ### Determine which simulations have overland flow...
  overland.files.dir <-  'C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/Calibration/CalibrationTest10/fluxes_stores/'
  setwd(overland.files.dir)
  files <- list.files(pattern='fluxqof.*\\.csv',recursive=T)
  
  for (file.name in files) {
    try.file <- read.csv(file.name)
    try.file <- try.file[,2:ncol(try.file)]
    if (sum(try.file)>0){
      print(file.name)
      which(try.file>0)
    }
  }
  
  try.file <- read.csv('fluxqof191.csv')
  try.file <- try.file[,2:ncol(try.file)]
  
  
  runoff.in <- try.file
  
  ##### Create a movie of the dynamics of the HRU scale storages...
  
  FR.HRU <- FR.Spatial$disc$hru  
  x11()
  plot(FR.HRU)
  groups <- FR.Spatial$disc$groups
  HRU.list <- as.numeric(groups$id)
  
  #replace.values <- data.frame(id=HRU.list,qof=t(runoff.in[5787,]))
  #FR.Runoff.plot <- subs(FR.HRU,replace.values)
  #x11()
  #plot(FR.Runoff.plot)
  
  event.runoff <- seq(4470,4490)
  event.runoff.2 <- seq(5784,5795)
  
  dates.ts <- time(result.optim.run$sim)
  times.event <- dates.ts[5784:5795]
  
  for (i in event.runoff.2) {
    
    print(i)
    replace.values <- data.frame(id=HRU.list,qof=t(runoff.in[i,]))
    FR.Runoff.plot <- subs(FR.HRU,replace.values)
    if(i==event.runoff.2[1]) {
      FR.runoff.stack <- FR.Runoff.plot
    } else {
      FR.runoff.stack <- stack(FR.runoff.stack,FR.Runoff.plot)
    }
    
  }
  
  #stream.net <- 'C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/1 GIS ANALYSIS/Shapes/FR4000StreamNet.shp'
  
  library(lattice)
  library(viridisLite)
  nlayers(FR.runoff.stack)
  x11()
  levelplot(FR.runoff.stack,names.attr=as.character(times.event),layout=c(4,3),scales=list(tck=c(1,1),draw=FALSE),xlim=c(5678500,5684000), ylim=c(3704500,3709000),
            colorkey=list(title='Runoff [m/hr]'), col.regions=viridis(100))
  x11()
  levelplot(FR.runoff.stack$X5792,names.attr=as.character(times.event[9]),layout=c(1,1),scales=list(tck=c(1,1),draw=FALSE),xlim=c(5678500,5684000), ylim=c(3704500,3709000),
            colorkey=list(title='Runoff [m/hr]'), col.regions=viridis(100))
  
  save.map <- 'C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/Figures/WetnessMaps'
  writeRaster(FR.runoff.stack$X5791,paste0(save.map,'/Runoff5791.tif'))
  save.map <- 'C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/Figures/WetnessMaps'
  writeRaster(FR.runoff.stack$X5792,paste0(save.map,'/Runoff5792.tif'))
  
  #####
  # PLot this such that the flows are all on the same axis
  x11()
  # Plot the FDCs
  
  for (r.name in 1:length(r.names)) {
    r.name=25
    order <- explicit.reach.table[r.name,5]
    if (order ==4){
      col='black'; lwd = 3
    } else if (order ==3) {
      col='#E69F00'; lwd = 2
    } else if (order == 2 ) {
      col ='#0072B2'; lwd =1;
    } else {
      col = '#CC79A7'; lwd =.1
    }
    
    
      #lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
      #    ylab='Q [m3/hr]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-03,1e05))
      plot(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',log='y',
            ylab='Q [mm/s]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-10,1e-1))
    
  }
  
  legend(70,max(Q.in.rank)*0.5, legend = c('Order 4','Order 3','Order 2','Order 1'),
         col=c('black','#E69F00','#0072B2','#CC79A7'),lwd=c(4,3,2,.1), title='Stream Order')
  
  
  ## Read in the 1/0 timeseries and compare to simulation
  logger.dir <- 'C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/3 DATA/2-PROJECT1-ROBINSON-FOREST/4-LOGGERDATA'
  files <- list.files(logger.dir,pattern='*clean.csv')
  # plot the logger data
  
  
    file.name <- files[4]
    short.name <- substr(file.name,start=1,stop=(nchar(file.name)-10))
    TS.logger <- data.frame(read.csv(paste0(logger.dir,'/',file.name)))
    #TS.logger.clean <- TS.logger %>% select(RoundTime,ApproxState)
    TS.logger.clean <- data.frame(Date=TS.logger$RoundTime,State=TS.logger$ApproxState)
    TS.logger.clean$State[which(TS.logger.clean$State=='Missing data')] <- 3
    TS.logger.clean$Date <- (ymd_hms(TS.logger.clean$Date))
    #sum(is.finite(TS.logger.clean$State))
    TS.logger.clean$Date2 <- as.Date(TS.logger.clean$Date)
    TS.logger.clean$State <- as.numeric(TS.logger.clean$State)
    
    TS.zoo <- zoo(TS.logger.clean$State,TS.logger.clean$Date)
    #TS.logger.clean <- data.frame(TS.logger.clean)
    #time(TS.logger.clean) <- time(TS.logger)
    
    #plot(TS.logger.clean$State)
    #image(TS.logger.clean$ApproxState, col=c('red','blue','yellow'))
    on <- sum(TS.logger.clean$State==1)
    off <- sum(TS.logger.clean$State==0)
    Total.logger <- on+off
    percent.on <- on/Total.logger
    percent.on
    
    ones <- rep(1,nrow(TS.logger.clean))
    colors = c("black", "#0072B2",'white')
    plot(TS.logger.clean$Date2,ones, type="h", col=colors[TS.zoo +1], ylim=c(0,1), xlab='date',
         xaxt='n',yaxt='n',lwd=5,ylab=short.name)
    ytick <- c(0,1)
    axis(side=2,at=ytick)
    axis.Date(1,at=TS.logger.clean$Date2,format = '%m-%y',tick=F)
    
    
    
  
  legend(1,10, legend = c('Dry','Wet','No data'),
         col=c('black','#0072B2','white'),lwd=c(2,2,2), title='Legend') 
  
  order.1.on <- ifelse(explicit.routing$routed.flow.instant.mm_s$`810`>2.804104e-06,1,0)
  order.1.on <- data.frame(order.1.on)
  time.sim <- time(result.optim.run$run$fluxes$qin)
  order.1.on$Date <- time.sim
  
  sim.on.off <- order.1.on[2218:14323,]
  
  logger.TS <- data.frame(state=TS.zoo,Date=TS.logger.clean$Date)
  try.match <- left_join(order.1.on,logger.TS,by='Date')
  try.match$state <- as.numeric(try.match$state)
  try.match$difference <- try.match$state-try.match$order.1.on
  try.match$difference <- ifelse(try.match$difference==3|try.match$difference==2,NA,try.match$difference)
  
  zero.difference <- sum(try.match$difference==0,na.rm=T)
  missed.on <- sum(try.match$difference==1,na.rm=T)
  missed.off <- sum(try.match$difference==-1,na.rm=T)
  correct.state.percent <- zero.difference/(zero.difference+missed.on+missed.off)*100
  correct.state.percent
  
  correct.on.off <- try.match[2218:14323,]
  correct.on.off$order.1.on <- ifelse(correct.on.off$state==3,3,correct.on.off$order.1.on)
  correct.on.off$Date2 <- as.Date(correct.on.off$Date)
  correct.on.off.zoo <- zoo(correct.on.off,correct.on.off$Date)
  logger.zoo <- zoo(correct.on.off$state,correct.on.off$Date)
  sim.onoff.zoo <- zoo(correct.on.off$order.1.on,correct.on.off$Date)
  
  
  
  x11()
  par(mfrow=c(2,1))
  
  ones <- rep(1,nrow(correct.on.off))
  colors = c("black", "#0072B2",'white')
  plot(correct.on.off$Date2,ones, type="h", col=colors[logger.zoo +1], ylim=c(0,1), xlab='date',
       xaxt='n',yaxt='n',lwd=1,ylab=short.name)
  ytick <- c(0,1)
  axis(side=2,at=ytick)
  axis.Date(1,at=correct.on.off$Date2,format = '%m-%y',tick=F)
  
  ones <- rep(1,nrow(correct.on.off))
  colors = c("black", "#0072B2",'white')
  plot(correct.on.off$Date2,ones, type="h", col=colors[sim.onoff.zoo +1], ylim=c(0,1), xlab='date',
       xaxt='n',yaxt='n',lwd=1,ylab='Simulation')
  ytick <- c(0,1)
  axis(side=2,at=ytick)
  axis.Date(1,at=correct.on.off$Date2,format = '%m-%y',tick=F)
  
  
  
  ## calculate FD of internal storages...
  hru.names <- groups$id
  Flow.duration.hru.qbf <- matrix(0,nrow = nrow(result.optim.run$run$fluxes$qbf),ncol = ncol(result.optim.run$run$fluxes$qbf))
  unsort.Flow.duration.hru.qbf <- matrix(0,nrow = nrow(result.optim.run$run$fluxes$qbf),ncol = ncol(result.optim.run$run$fluxes$qbf))
  Flow.duration.hru.sd <- matrix(0,nrow=nrow(result.optim.run$run$storages$sd),ncol=ncol(result.optim.run$run$storages$sd))
  qbf.in.rank <- matrix(0,nrow = nrow(result.optim.run$run$fluxes$qbf),ncol = ncol(result.optim.run$run$fluxes$qbf))
  sd.in.rank <- matrix(0,nrow=nrow(result.optim.run$run$storages$sd),ncol=ncol(result.optim.run$run$storages$sd))
  colnames(Flow.duration.hru.qbf) <- hru.names
  colnames(unsort.Flow.duration.hru.qbf) <- hru.names
  colnames(qbf.in.rank) <- hru.names
  colnames(Flow.duration.hru.sd) <- hru.names
  colnames(sd.in.rank) <- hru.names
  
  qbf <- matrix(result.optim.run$run$fluxes$qin,nrow=nrow(result.optim.run$run$fluxes$qbf),ncol = ncol(result.optim.run$run$fluxes$qbf))
  qbf.2 <- matrix(result.optim.run$run$fluxes$qbf,nrow=nrow(result.optim.run$run$fluxes$qbf),ncol = ncol(result.optim.run$run$fluxes$qbf))
  
  #sd <- matrix(result.optim.run$run$storages$sd,nrow=nrow(result.optim.run$run$storages$sd),ncol=ncol(result.optim.run$run$storages$sd))
  #sd.2 <- result.optim.run$run$storages$sd
  ## Calculate FDC
  for (hru.name in 1:length(hru.names)) {
    if (hru.name >1) {
    qbf.in <- qbf[,hru.name]
    qbf.in.sort <- data.frame(flow.asc=sort(qbf.in,decreasing =T))
    #sd.in <- sd[,hru.name]
    #sd.in.sort <- data.frame(sd.asc=sort(sd.in,decreasing=T))
    
    
    rank.qbf <- 1:length(qbf.in)
    qbf.in.sort$rank <- rank.qbf
    #rank.sd <- 1:length(sd.in)
    #sd.in.sort$rank <- rank.sd
    
    qbf.in.sort$Prob.qbf.in <- 100*(qbf.in.sort$rank/(length(qbf.in)+1))
   # sd.in.sort$Prob.sd.in <- 100*(sd.in.sort$rank/length(sd.in)+1)
    
    
    
    Flow.duration.hru.qbf[,hru.name] <- qbf.in.sort$Prob.qbf.in 
    qbf.in.rank[,hru.name] <- qbf.in.sort$flow.asc
    
    #Flow.duration.hru.sd[,hru.name] <- qbf.in.sort$Prob.qbf.in 
    #sd.in.rank[,hru.name] <- sd.in.sort$flow.asc
    
    
    plot(x=as.numeric(Flow.duration.hru.qbf[,hru.name]),y=as.numeric(qbf.in.rank[,hru.name]), type ='l',
         log='y',ylab='Q [mm/s]',xlab='% Exceedance',col='black',lwd=2)
    par(new=TRUE)
    #plot(x=as.numeric(Flow.duration.hru.sd[,hru.name]),y=as.numeric(sd.in.rank[,hru.name]), type ='l',
    #     log='y',ylab='Q [mm/s]',xlab='% Exceedance',col='black',lwd=2)
    
    foo <- data.frame(qbf=qbf.in,rank=rank(1/qbf.in))
    
    foo$probability <- foo$rank/(1+nrow(foo))*100
    unsort.Flow.duration.hru.qbf[,hru.name] <- foo$probability
    }
  }
  
  
  FR.HRU <- FR.Spatial$disc$hru  
  groups <- FR.Spatial$disc$groups
  HRU.list <- as.numeric(groups$id)
  
  #replace.values <- data.frame(id=HRU.list,qof=t(runoff.in[5787,]))
  #FR.Runoff.plot <- subs(FR.HRU,replace.values)
  #x11()
  #plot(FR.Runoff.plot)
  
  event.runoff <- 2121
  

    replace.values.qin <- data.frame(id=HRU.list,qin=(qbf[event.runoff,]))
    FR.qbf.plot <- subs(FR.HRU,replace.values.qin)
    
    replace.values.qbf <- data.frame(id=HRU.list,qbf=(qbf.2[event.runoff,]))
    FR.qbf2.plot <- subs(FR.HRU,replace.values.qbf)
    
    replace.values.sd <- data.frame(id=HRU.list,sd=(sd[event.runoff,]))
    FR.sd.plot <- subs(FR.HRU,replace.values.sd)
  x11()
  plot(FR.qbf.plot)
  
  x11()
  plot(FR.sd.plot)
  
  save.map <- 'C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/Figures/WetnessMaps'
  writeRaster(FR.qbf.plot,paste0(save.map,'/qbf_plot.tif'))
  writeRaster(FR.sd.plot,paste0(save.map,'/sd_plot.tif'))
  writeRaster(FR.qbf2.plot,paste0(save.map,'/qbf2_plot.tif'))
  getwd()
}


# plot the FDC

setwd('C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/Calibration/CalibrationTest13/fluxes_stores')
qin.run <- read.csv(paste0('fluxqin',qin.file,'.csv'))                           # Read qin 
qin.run <- xts(qin.run[,-1],time.qbf)                                            # Convert to xts
qbf.run <- read.csv(paste0('fluxqbf',qin.file,'.csv'))                           # Read qbf
qbf.run <- xts(qbf.run[,-1],time.qbf)                                            # Convert to xts
ae.run <- read.csv(paste0('fluxae',qin.file,'.csv'))                             # Read ae
ae.run <- xts(ae.run[,-1],time.qbf)                                              # Convert to xts
rain.run <- read.csv(paste0('fluxrain',qin.file,'.csv'))                         # Read rain
rain.run <- xts(rain.run[,-1],time.qbf)                                          # Convert to xts
qof.run <- read.csv(paste0('fluxqof',qin.file,'.csv'))                           # Read qof
qof.run <- xts(qof.run[,-1],time.qbf)                                            # Convert to xts
Qsim.run <- read.csv(paste0('qsim',qin.file,'.csv'))                             # Read Qsim
Qsim.run <- xts(Qsim.run[,-1],time.qbf)                                          # Convert to xts
q_specific.run <- read.csv(paste0('qsim_specific',qin.file,'.csv'))
q_specific.run <- data.frame(q_specific.run[,-1]*24000)
compile.fluxes <- list(qin = qin.run, qbf = qbf.run, ae = ae.run,                # Compile fluxes 
                       rain = rain.run, qof = qof.run)                           # Into list
compile.run <- list(fluxes=compile.fluxes, Qsim = Qsim.run)                      # Compile fluxes and Qsim into list

# Run routing function -- note TOPMODEL puts out data in m/hr for reach-specific flow in runoff
explicit.routing.run <- explicit.routing.instant(FR.Spatial,                     # Read in FR.Spatial and the compiled run
                                                 run=compile.run,
                                                 model.timestep = model.timestep)# Run the explicit routing code for the simulation
Q <- data.frame(explicit.routing.run$routed.flow.instant.m3_hr)                   # Assign the routed flow instant (mm/s) to a variables
r.names <- names(Q)                                                              # Get the reach names

# Convert Q from m/hr to mm/day
#Q <- Q*24*60*60                                                                # This is important for some of the empirical calculations we'll later need to do

#Q.compare.test <- data.frame(sim_routing=Q$X3802, sim_dtop=q_specific.run,diff = (Q$X3802-q_specific.run)/q_specific.run*100)

# Calculate the fdc of each reach 
Flow.duration.reaches <- matrix(0,nrow = length(Q[,1]),ncol = length(r.names))   # initialize the flow duration matrix
Q.in.rank <- matrix(0,nrow = length(Q[,1]),ncol = length(r.names))               # Initialize the Q in rank matrix
colnames(Flow.duration.reaches) <- r.names                                       # Assign names to the flow duration matrix
colnames(Q.in.rank) <- r.names                                                   # Assign names to the Q in rank matrix

# Run the FDC code for the reaches
for (r.name in 1:length(r.names)) {                                              # For each of the reaches
  if (r.name != 26) {                                                            # Skip reach 26 because it isn't working properly
    Q.in <- Q[,r.name]                                                           # Assing the Q in for a reach
    order <- explicit.reach.table[r.name,5]                                      # Get the stream order
    Q.in.sort <- data.frame(flow.asc=sort(Q.in,decreasing =T))                   # Rank/sort the Q in
    rank <- 1:length(Q.in)                                                       # Get the rank
    Q.in.sort$rank <- rank                                                       # assign the rank for the sorted data dataframe
    Q.in.sort$Prob.Q.in <- 100*(Q.in.sort$rank/(length(Q.in)+1))                 # Calculate the probablility
    Flow.duration.reaches[,paste0('X',as.character(explicit.reach.table$link_no[r.name]))] <- Q.in.sort$Prob.Q.in                                # Record the probability 
    Q.in.rank[,paste0('X',as.character(explicit.reach.table$link_no[r.name]))] <- Q.in.sort$flow.asc                                             # Record the ascending flow
    order <- explicit.reach.table[r.name,5]                                      # Determine the order 
  }
}
Flow.duration.reaches <- data.frame(Flow.duration.reaches)                       # Write flow duration as a dataframe
Q.in.rank <- data.frame(Q.in.rank)      
x11()
for (r.name in 1:ncol(Q.in.rank)) {
  order <- explicit.reach.table[r.name,5]
  if (order ==4){
    col='black'; lwd = 3
  } else if (order ==3) {
    col='#0072B2'; lwd = 2
  } else if (order == 2 ) {
    col ='#E69F00'; lwd =1;
  } else {
    col = '#CC79A7'; lwd =.1
  }
  
  if (r.name==1) { 
    #plot(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
    #     log='y',ylab='Q [m3/hr]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-03,1e05))
    plot(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
         log='y',ylab='Q [mm/s]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(0.001,1000),
         cex =1.5,cex.axis=1.5,cex.lab=1.5)
    
  } else {
    #lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
    #    ylab='Q [m3/hr]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(1e-03,1e05))
    lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(Q.in.rank[,r.name]), type ='l',
          ylab='Q [mm/s]',xlab='% Exceedance',col=col,lwd=lwd,ylim=c(0.001,1000))
  }
}


# Plot the average of each stream order 
stream.orders <- explicit.reach.table[,5]


order.1 <- (Q.in.rank[,which(stream.orders==1)])
order.1 <- rowMeans(order.1[,-10])
order.2 <- rowMeans(Q.in.rank[,which(stream.orders==2)])
order.3 <- rowMeans(Q.in.rank[,which(stream.orders==3)])
order.4 <- rowMeans(Q.in.rank[,which(stream.orders==4)])


x11()
plot(x=as.numeric(Flow.duration.reaches[,1]),y=as.numeric(order.4), type ='l',
     log='y',ylab='Q [m^3/d]',xlab='% Exceedance',col='black',lwd=2,ylim=c(1e-2,2000))

lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(order.3), type ='l',
      ylab='Q m^3/d]',xlab='% Exceedance',col='#E69F00',lwd=2,ylim=c(1e-2,2000))

lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(order.2), type ='l',
      ylab='Q [m^3/d]',xlab='% Exceedance',col='#0072B2',lwd=2,ylim=c(1e-2,2000))

lines(x=as.numeric(Flow.duration.reaches[,r.name]),y=as.numeric(order.1), type ='l',
      ylab='Q [m^3/d]',xlab='% Exceedance',col='#CC79A7',lwd=2,ylim=c(1e-2,2000))

legend(70,max(Q.in.rank)*0.25, legend = c('Order 4','Order 3','Order 2','Order 1'),
       col=c('black','#E69F00','#0072B2','#CC79A7'),lwd=c(4,3,2,1), title='Stream Order') 


fdc.plot

outlet.flow <- data.frame(Flow.duration=Flow.duration.reaches$X3802,
                          Q.in=Q.in.rank$X3802)

# Plot the average probability of flow presence for each reach
fdc.plot <- ggplot(data=outlet.flow, aes(x=Flow.duration,y=Q.in))+geom_line()+scale_y_continuous(trans='log10')
fdc.plot