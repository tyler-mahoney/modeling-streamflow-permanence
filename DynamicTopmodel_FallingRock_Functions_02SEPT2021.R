# Function to read in the results of a run, run explicit routing code, output the fdc of reaches, estimate the 'threshold flow',
# estimate the time that the threshold flow is exceeded, determine the timing of the 'on' and 'off' of the sensor...
# Inputs: qbf/Qin for each HRU for a behavioral simulation, FR.Spatial configuration, sensor inputs
Headwater.evaluation.dynamic <- function(qin.file, FR.Spatial, TS.logger.clean.all, TS.zoo.all, logger.out.percent, Reach.identifiers, Reach.identifiers.second,model.timestep, result.optim.run, time.qbf, ...) {
  # Read in explicit HRU results for fluxes and storages from the behavioral sets
  # NOTE WE MUST SET THE DIRECTORY TO THE PROPER CALIBRATION TEST FOLDER, AS OF 7/6/2021 IT IS CalibrationTest13
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
  compile.fluxes <- list(qin = qin.run, qbf = qbf.run, ae = ae.run,                # Compile fluxes 
                         rain = rain.run, qof = qof.run)                           # Into list
  compile.run <- list(fluxes=compile.fluxes, Qsim = Qsim.run)                      # Compile fluxes and Qsim into list
  
  # Run routing function -- note TOPMODEL puts out data in m/hr for reach-specific flow in runoff
  explicit.routing.run <- explicit.routing.instant(FR.Spatial,                     # Read in FR.Spatial and the compiled run
                                                   run=compile.run,
                                                   model.timestep = model.timestep)# Run the explicit routing code for the simulation
  Q <- data.frame(explicit.routing.run$routed.flow.instant.mm_s)                   # Assign the routed flow instant (mm/s) to a variables
  r.names <- names(Q)                                                              # Get the reach names
  
  # Convert Q to mm/day
  Q <- Q*24000                                                                  # This is important for some of the empirical calculations we'll later need to do
  
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
  Q.in.rank <- data.frame(Q.in.rank)                                               # Write flow duration as a dataframe  
  
  # Write a function that returns the flow associated with the flow threshold percentage desired by the user
  nearest <- function(want,have) {                                                 # Function's name is nearest, two arguments in
    near <- which(abs(want-have)==min(abs(want-have)))                             # Return the percentage that's nearest to the desired percentage threshold
    return(near)                                                                   # Return the percentage
  }
  
  x=Sys.time()
  ## Run the analysis for 1000 flow thresholds, output the optimum flow threshold. 
  min.threshold <- 0.05                                                           # This was the minimum transmissivity scaling exponent measured from the Prancevic and Kirchner 2018 Paper
  max.threshold <- 10                                                           # This was the maximum transmissivity scaling exponent measured from the Prancevic and Kirchner 2018 Paper
  gamma.range <- seq(from=log(min.threshold), to = log(max.threshold),                       # Creates a sequency of 1000 thresholds for which we'll calculate the % correct
                     by = (log(max.threshold)-log(min.threshold))/500)                 # Continued 
  
  # Format the data 
  names.loggers <- names(TS.logger.clean.all)
  
  # Initialize the lists for each of the loggers. 
  logger.TS <- list()                                                              # Initialize the logger Timeseries list
  logger.timeseries <- list()
  #sim.obs.difference <- list()
  confusion.matrix.thresh <- list()
  out.thresh.df <- list()
  reach.slope <- data.frame(matrix(nrow=4))
  reach.contributing.area <- data.frame(matrix(nrow=4))
  
  for (name.logger in 1:length(names.loggers)) {
    # Format the time series properly
    time.sim <- time(result.optim.run$run$fluxes$qin)                                # Get the time stamp for the qsim
    time.sim <- data.frame(Date=time.sim)                                            # Set the col name to 'Date'
    logger.TS[[name.logger]] <- data.frame(state=TS.zoo.all[[name.logger]],
                                           Date=TS.logger.clean.all[[name.logger]]$Date)                  # Get the logger time series
    
    # Join the logger.TS to the time.sim 
    logger.timeseries[[name.logger]] <- left_join(time.sim,logger.TS[[name.logger]],by='Date')                                    # Join the logger data to the time series
    logger.timeseries[[name.logger]]$state <- ifelse(logger.timeseries[[name.logger]]$state==3,NA,logger.timeseries[[name.logger]]$state)        # Set NA values
    logger.timeseries[[name.logger]] <- logger.timeseries[[name.logger]][complete.cases(logger.timeseries[[name.logger]]),]                      # Remove rows with NA values
    
    # initialize the matrix that records the difference in simulated and observed
    sim.obs.difference[[name.logger]] <- data.frame(
      matrix(nrow=length(logger.timeseries[[name.logger]]$Date),ncol=length(gamma.range)+1))
    
    # Get the contributing area and slope for the reach 
    reach.name.dynamic <- which(                                                 # Determine which row the reach name of the logger belongs to
      FR.Spatial$explicit.ChanTable$link_no==                                    # Compare the list of reach numbers to the current logger being evaluated
        as.numeric(sub('.','',Reach.identifiers[name.logger])))                  # Define the reach name as a numeric
    reach.slope[name.logger,1] <-                                                  # Determine the slope of the reach (ft/ft) or (m/m)
      FR.Spatial$explicit.ChanTable$stream_slope_ftpft[reach.name.dynamic]       # Using the which function from above
    reach.contributing.area[name.logger,1] <-                                      # Determine the contributing area (m^2) using the which function from aboce
      FR.Spatial$explicit.ChanTable$US_area_m2[reach.name.dynamic]               # Using the function from above
    
    # Convert contributing area to km2
    reach.contributing.area[name.logger,1] <- 
      reach.contributing.area[name.logger,1]*1e-6
    
    # Run the analysis to simulate the 1/0 time series for the four sensors over the 1000 transmissivity scaling exponents
    for (gamma.in in 1:length(gamma.range)) {                                              
      # Write the simulated 1/0 time series for the average threshold
      
      # Determine the subsurface flow capacity of the valley below the reach. This equation is proposed in Godsney and Kirchner (2018) and Prancevic and Kirchner (2014)
      transmissivity <- exp(gamma.range[gamma.in])
      thresh <- transmissivity*reach.slope[name.logger,1]
      
      #percent.thresh <- Flow.duration.reaches[which(r.names==Reach.identifiers[name.logger])][nearest(Q.in.rank[which(r.names==Reach.identifiers[name.logger])],         # Look at the flow duration for reach 810 (associated with FC4) and 
      #                                                                                                thresh),1]                            # Look at the flow duration for reach 810 (associated with FC4) and 
      explicit.r.names <- names(explicit.routing.run$routed.flow.instant.m3_hr)
      sim.810.on <- 
        ifelse(explicit.routing.run$routed.flow.instant.m3_hr[which(explicit.r.names==Reach.identifiers.second[name.logger])]>thresh,           # USing the iterative threshold
               1,0)                                                               # Set the values greater than the threshold to one, otherwise set to zero 
      sim.810.on <- data.frame(sim.810.on)                                           # Convert to a data frame
      time.sim <- time(result.optim.run$run$fluxes$qin)                              # Get the time stamp for the qsim
      sim.810.on$Date <- time.sim                                                    # Set the date of the data frame to the time.sim 
      
      # Compare the logger 1/0 to the simulated 1/0
      sim.obs.match <- left_join(logger.timeseries[[name.logger]],sim.810.on,by='Date')                           # Join the logger time series to the simulated time series by the date
      sim.obs.match$state <- as.numeric(sim.obs.match$state)                                   # Write the state as a numeric
      sim.obs.match[,3] <- as.numeric(sim.obs.match[,3])
      sim.obs.match$difference <- sim.obs.match$state-sim.obs.match[,3]                     # Subtract simulated state from the observed state. state is observed, order.1.on is simulated; if the difference is 0 then the sim and obs match. If it's 1 then the obs = 1 and sim = 0. If it's -1 then obs = 0 and sim = 1
      sim.obs.match$difference <- ifelse(sim.obs.match$difference==3|
                                           sim.obs.match$difference==2,  # If the difference is 3 or 2, then this should just be set to NA. I set observed NA values to 3 before, not sure why but it makes the code work.
                                         NA,sim.obs.match$difference)                          # If it's 3 or 2 then it will just be NA, otherwise it will be the correct difference.
      
      # Sum the number of correct, missed on, missed off values 
      zero.difference <- sum(sim.obs.match$difference==0,na.rm=T)                          # Zero values mean sim and obs state match; this is the sum of zero difference time steps
      true.on <- sum(ifelse(sim.obs.match$state==1 & sim.obs.match[,3]==1, # See if both obs and model equal 1
                            1,0),na.rm=T)                                                    # Calculate the true positives
      true.off <- sum(ifelse(sim.obs.match$state==0 & sim.obs.match[,3]==0,# See if both obs and model equal 0
                             1,0),na.rm=T)                                                   # Calculate the true negatives
      missed.on <- sum(sim.obs.match$difference==1,na.rm=T)                                # 1 means the obs = 1 and sim = 0; so the model missed the sensor being 'on' or 'wet' since 1 - 0 = 1
      missed.off <- sum(sim.obs.match$difference==-1,na.rm=T)                              # -1 means the obs = 0 and the sim = 1; so the model missed the sensor being 'off' or 'dry' since 0 - 1 = -1
      
      correct.state.percent.thresh <-                                                         # The correct percent is the total number of times when the difference is zero
        zero.difference/(zero.difference+missed.on+missed.off)*100                     # divided by the total number of timesteps when the sensor has data 
      
      flow.permanence.thresh <- sum(sim.obs.match[,3]==1)/length(sim.obs.match[,3])     # Flow permanence is the sum of time steps when the sensor is 'on' or 'wet' divided by the total number of time steps 
      
      if(gamma.in==1) {
        confusion.matrix.thresh[[name.logger]] <- data.frame(true.on=true.on,                         # Calc the confusion matrix
                                                             true.off=true.off,                      # True off  
                                                             missed.on=missed.on,                    # missed on
                                                             missed.off=missed.off)                  # missed off
      } else { 
        confusion.matrix.join.thresh <- data.frame(true.on=true.on,                    # Join to the previous matrix      
                                                   true.off=true.off,                       # true off
                                                   missed.on=missed.on,                     # missed on
                                                   missed.off=missed.off)                   # missed off
        confusion.matrix.thresh[[name.logger]] <- rbind(confusion.matrix.thresh[[name.logger]],                       # Join to the previous
                                                        confusion.matrix.join.thresh)                  # Just created join matrix
      }
      
      # Calculate the percent thresh 
      
      if (gamma.in==1) {
        out.thresh.df[[name.logger]] <- data.frame(percent.correct=correct.state.percent.thresh,
                                                   flow.thresh=thresh, sum.correct=sum(true.on,true.off), # percent.thresh=percent.thresh,
                                                   sum.incorrect=sum(missed.on,missed.off),
                                                   transmissivity = transmissivity)
      } else {
        out.thresh.bind <- data.frame(percent.correct=correct.state.percent.thresh,
                                      flow.thresh=thresh, sum.correct=sum(true.on,true.off), # percent.thresh=percent.thresh,
                                      sum.incorrect=sum(missed.on,missed.off),
                                      transmissivity = transmissivity)
        out.thresh.df[[name.logger]] <- rbind(out.thresh.df[[name.logger]],out.thresh.bind)
      }
      
      # Write out the difference into a big matrix.
      #if (gamma.in==1) {
      #  sim.obs.difference[[name.logger]][,1] <- sim.obs.match$Date
      #  sim.obs.difference[[name.logger]][,gamma.in+1] <- sim.obs.match$difference
      #} else {
      #  sim.obs.difference[[name.logger]][,gamma.in+1] <- sim.obs.match$difference
      #}
    }
  }
  Sys.time()-x
  # Calculate the respective error for each flow thresh for each sensor
  residual.on.off <- matrix(nrow=length(gamma.range),ncol=length(names.loggers))
  residual.on.off[,1] <- confusion.matrix.thresh[[1]][,3]+
    confusion.matrix.thresh[[1]][,4]
  residual.on.off[,2] <- confusion.matrix.thresh[[2]][,3]+
    confusion.matrix.thresh[[2]][,4]
  residual.on.off[,3] <- confusion.matrix.thresh[[3]][,3]+
    confusion.matrix.thresh[[3]][,4]
  residual.on.off[,4] <- confusion.matrix.thresh[[4]][,3]+
    confusion.matrix.thresh[[4]][,4]
  
  # Calculate the total residual error by summing the error for each sensor
  rowSums.residual.on.off <- data.frame(sum.error=rowSums(residual.on.off),transmissivity=exp(gamma.range), 
                                        thresh.1 = out.thresh.df[[1]]$flow.thresh,
                                        thresh.2 = out.thresh.df[[2]]$flow.thresh,
                                        thresh.3 = out.thresh.df[[3]]$flow.thresh,
                                        thresh.4 = out.thresh.df[[4]]$flow.thresh)
  
  # calculate which threshold produces the minimum error accross all four sensors
  thresh.best <- (rowSums.residual.on.off$transmissivity[which.min(rowSums.residual.on.off$sum.error)])
  
  best.residual.on.off <- rowSums.residual.on.off[which(rowSums.residual.on.off$transmissivity==thresh.best),]
  # Output results for all thresholds into a list
  flow.thresh.best.results <- 
    list(out.thresh.data=out.thresh.df, #sim.obs.difference=sim.obs.difference, 
         confusion.matrix=confusion.matrix.thresh,rowSums.residual.on.off=rowSums.residual.on.off)
  
  # 
  for (name.logger in 1:length(names.loggers)) {
    if (name.logger ==1) {
      save.out.thresh.best <- data.frame(out.thresh.df[[name.logger]][which(out.thresh.df[[name.logger]]$flow.thresh==thresh.best),])
    }else {
      save.out.thresh.best.bind <- out.thresh.df[[name.logger]][which(out.thresh.df[[name.logger]]$flow.thresh==thresh.best),]
      save.out.thresh.best <- rbind(save.out.thresh.best,save.out.thresh.best.bind)
    }
  }
  
  # Calculate the percent of time that the network is on and the percent of time the network is off
  # Also calcualte the time series of on/off 
  Q.sub.c.best <- thresh.best*FR.Spatial$explicit.ChanTable$stream_slope_ftpft
  total.network.on.off <- data.frame(matrix(nrow=nrow(explicit.routing.run$routed.flow.instant.m3_hr),ncol=ncol(explicit.routing.run$routed.flow.instant.m3_hr)))
  colnames(total.network.on.off) <- r.names
  for (reach in 1:length(Q.sub.c.best)) {
    total.network.on.off[,reach] <- ifelse(explicit.routing.run$routed.flow.instant.m3_hr[,reach]>Q.sub.c.best[reach],1,0)
  }
  
  flow.network.thresh <- explicit.routing.run$routed.flow.instant.m3_hr
  colnames(flow.network.thresh) <- r.names
  for (reach in 1:length(Q.sub.c.best)) {
    flow.network.thresh[flow.network.thresh[,reach] < (Q.sub.c.best[reach]),reach] <- 0
  }
  
  
  # Calculate the percent of time that each reach is considered to be on
  reach.names.list <- colnames(total.network.on.off)
  
  percent.on.network <- c()
  for (name.iter in 1:length(reach.names.list)) {
    percent.on.network[name.iter] <- 
      sum(total.network.on.off[,name.iter])/length(total.network.on.off[,name.iter])
  }
  
  #best.percent.thresh <- out.thresh.df[[4]]$percent.thresh[which(out.thresh.df[[4]]$flow.thresh==(thresh.best))]
  
  # Write the simulated 1/0 time series for the best threshold
  #order.1.on <- 
  #  ifelse(explicit.routing.run$routed.flow.instant.m3_hr$`810`>thresh.best, # USing the averaged threshold
  #         1,0)                                                                    # Set the values greater than the threshold to one, otherwise set to zero 
  #order.1.on <- data.frame(order.1.on)                                             # Convert to a data frame
  #time.sim <- time(result.optim.run$run$fluxes$qin)                                # Get the time stamp for the qsim
  #order.1.on$Date <- time.sim                                                      # Set the date of the data frame to the time.sim 
  
  # Compare the logger 1/0 to the simulated 1/0
  #final.match <- left_join(order.1.on,logger.timeseries[[4]],by='Date')                           # Join the logger time series to the simulated time series by the date
  #final.match$state <- as.numeric(final.match$state)                                   # Write the state as a numeric
  #final.match$difference <- final.match$state-final.match$order.1.on                     # Subtract simulated state from the observed state. state is observed, order.1.on is simulated; if the difference is 0 then the sim and obs match. If it's 1 then the obs = 1 and sim = 0. If it's -1 then obs = 0 and sim = 1
  #final.match$difference <- ifelse(final.match$difference==3|final.match$difference==2,  # If the difference is 3 or 2, then this should just be set to NA. I set observed NA values to 3 before, not sure why but it makes the code work.
  #                                 NA,final.match$difference) 
  
  # Sum the number of correct, missed on, missed off values for 2003 
  #zero.difference.yr2003 <- sum(final.match$difference[2218:4473]==0,na.rm=T)        # This is the time when the sensor is actively collecting data between 2003-2004; sum the values equal to 0
  #missed.on.yr2003 <- sum(final.match$difference[2218:4473]==1,na.rm=T)              # Sum the values equal to 1
  #missed.off.yr2003 <- sum(final.match$difference[2218:4473]==-1,na.rm=T)            # Sum the values equal to 0
  #correct.state.percent.yr2003 <- zero.difference.yr2003/                          # Calculate correct percent for 2003
  #  (zero.difference.yr2003+missed.on.yr2003+missed.off.yr2003)*100                # number of time steps equal to zero divided by total number of timesteps for 2003
  #flow.permanence.yr2003.sim <-                                                    # Simulated 2003 'on' divided by total 
  #  sum(final.match$order.1.on[2218:4473]==1)/length(final.match$order.1.on[2218:4473])# 2003 time steps
  #flow.permanence.yr2003 <-                                                        # observed 2003 fraction 'on' 
  #  sum(final.match$state[2218:4473]==1,na.rm=T)/                                    # Total 1 divided by 
  #  ((sum(final.match$state[2218:4473]==1,na.rm=T))+                                 # Total 1 plus
  #     (sum(final.match$state[2218:4473]==0,na.rm=T)))                               # Total 0
  
  # Sum the number of correct, missed on, missed off values for 2005
  #zero.difference.yr2005 <- sum(final.match$difference[10135:14323]==0,na.rm=T)      # This is the time when the sensor is actively collecting data between 2005-2006; sum the values equal to 0
  #missed.on.yr2005 <- sum(final.match$difference[10135:14323]==1,na.rm=T)            # Sum the values equal to 1
  #missed.off.yr2005 <- sum(final.match$difference[10135:14323]==-1,na.rm=T)          # Sum the values equal to -1
  #correct.state.percent.yr2005 <- zero.difference.yr2005/                            # Calculate correct percent for 2005
  #  (zero.difference.yr2005+missed.on.yr2005+missed.off.yr2005)*100                  # number of time steps equal to zero divided by total number of time steps
  #flow.permanence.yr2005.sim <-                                                      # Simulated 2005 'on' divided by 
  #  sum(final.match$order.1.on[10135:14323]==1)/                                     # total 2005 simulated
  #  length(final.match$order.1.on[10135:14323])                                      # timesteps
  #flow.permanence.yr2005 <-                                                          # Observed 2005 fraction 'on'  
  #  sum(final.match$state[10135:14323]==1,na.rm=T)/                                  # total 1 divided by 
  #  ((sum(final.match$state[10135:14323]==1,na.rm=T))+                               # Total 1 plus
  #     (sum(final.match$state[10135:14323]==0,na.rm=T)))                             # Total 0
  
  # Calculate statistics for the flow exceeded, correct state percent for the total period, 2003, and 2005, flow permanence
  #stats.on.off <- data.frame(                                          # Write flow exceeded and correct.state.percent
  #  correct.percent.2003=correct.state.percent.yr2003, correct.percent.2005 = correct.state.percent.yr2005,               # Write the correct state for 2003 and 2005
  #  flow.permanence.2003.sim=flow.permanence.yr2003.sim,                           # Write the obs and sim flow permanences
  #  flow.permanence.2005.sim=flow.permanence.yr2005.sim,flow.permanence.2003=flow.permanence.yr2003,                      # Write the obs and sim flow permanences
  #  flow.permanence.2005=flow.permanence.yr2005) #best.percent.exceeded = best.percent.threshold,                                                                          # Write the obs and sim flow permanences
  
  # Compile everything as a list
  out.headwater.evaluation <- list(                                       # Compile a list of everything we want to output
    flow.thresh.best.results=flow.thresh.best.results, # stats.on.off=stats.on.off,
    thresh.best=thresh.best, #best.percent.thresh=best.percent.thresh,
    total.network.on.off=total.network.on.off,
    flow.network.thresh=flow.network.thresh, percent.on.network=percent.on.network, 
    save.out.thresh.best= save.out.thresh.best,
    best.residual.on.off=best.residual.on.off,
    Q.sub.c.best=Q.sub.c.best
  )                                                                                      # Continued
  
  return(out.headwater.evaluation)                                                 # return the list and exit the function
}



#################################################################
############ FUNCTIONS FOR DYNATOPMODEL CALIBRATION #############
#################################################################

### Wrapper function for running the dynatopmodel: 
Dynatophydromod <- function(param.values,                                        # Parameter vector containing values of dynatopmodel parameters
                            inner.timestep,                                      # Number of inner timesteps
                            rain,                                                # Timeseries of precipitation at dt to run the model (continuous data)
                            PET,                                                 # Timeseries of PET at dt needed to run the model (continuous data)
                            obs,                                                 # Timeseries of observed Q data to compare the simulated results 
                            disc,                                                # Spatial discretization (weighting matrix) and 'groups' parameter matrix 
                            RoutingTable=NULL,                                   # Routing table to discretize stream network
                            date,                                                # Optional
                            dates.df,                                            # Date frame of dates used for warm up and calibration periods 
                            model.drty,                                          # need the model directory
                            ...
) {
  
  # Run the simulation
  simLump <- runPSOCalibrationTopmodel(param.values=param.values,                         # Parameter vector containing values of dynatopmodel parameters
                                       inner.timesteps=inner.timestep,                    # Number of inner timesteps
                                       rains = rain,                                      # Timeseries of precipitation at dt to run the model (continuous data)
                                       PETs = PET,                                        # Timeseries of PET at dt needed to run the model (continuous data)
                                       obss = obs,                                        # Timeseries of observed Q data to compare the simulated results 
                                       discs = disc,                                      # Spatial discretization (weighting matrix) and 'groups' parameter matrix 
                                       RoutingTables = NULL,                              # Routing table to discretize stream network
                                       dates.dfs = dates.df)                              # Initial and Final date for the calibration period)                                 
  
  # Process model results
  
  obs <- window(obs,start=dates.df$calibration.initial,
                end = dates.df$calibration.final)
  qsim <- window(simLump$sim, start = dates.df$calibration.initial,
                 end = dates.df$calibration.final)
  if (length(obs)<length(qsim)) {
    qsim <- qsim[1:length(obs)]
  } else {
    obs <- obs[1:length(qsim)]
  }
  
  n <- length(qsim)                                                         # Removing the warming up period
  gof <- KGE(sim=as.numeric(log(qsim+1e-07)), obs=as.numeric(log(obs+1e-07), method="2012"))                             # Caclculate KGE of the log of the ts. added a small amount to avoid zeros. 
  # gof <- KGE(sim=as.numeric(log(qsim+.0001)), obs=as.numeric(log(obs+.0001), method="2012"))
  
  # If the gof is greater than 0.3, save the output of the run - NOTE THIS SAVES A TIME SERIES THAT INCLUDES THE WARM UP PERIOD
  if (gof>0.3) {
    save.out <- paste0(model.drty,'/fluxes_stores/')
    fluxes <- simLump$run$fluxes
    names.fluxes <- names(fluxes)
    for (flux in 1:length(fluxes)) {
      flux.name <- names.fluxes[flux]
      for (iter in 1:600) { # note: 600 is arbitrary, right now i'm only running the model 500 times. if a file name is greater than 250 then there's a problem somewhere. 
        out.flux.file <- paste0(save.out,'flux',flux.name,iter,'.csv')
        if (!file.exists(out.flux.file)) {
          write.csv(fluxes[[flux]],out.flux.file)
          break
        }
      }
    }
    stores <- simLump$run$storages
    names.stores <- names(stores)
    for (store in 1:length(stores)) {
      store.name <- names.stores[store]
      for (iter in 1:600) { 
        out.store.file <- paste0(save.out,'store',store.name,iter,'.csv')
        if (!file.exists(out.store.file)) {
          write.csv(stores[[store]],out.store.file)
          break
        }
      }
    }
    
    for (iter in 1:600) {  
      out.q.file <- paste0(save.out,'qsim',iter,'.csv')
      if (!file.exists(out.q.file)) {
        #write.csv(qsim,out.q.file)
        write.csv(simLump$run$Qsim[1:length(simLump$run$fluxes$qbf[,1])],paste0(save.out,'Qsim',iter,'.csv'))
        write.csv(simLump$run$qsim[1:length(simLump$run$fluxes$qbf[,1])],paste0(save.out,'qsim_specific',iter,'.csv'))
        write.csv(gof,paste0(save.out,'gof_log_kge',iter,'.csv'))
        
        break
      }
    }
  }
  # Creating the output of the R function
  nelements <- 2                                                            # Elements to include in the return variable
  out <- vector("list", nelements)                                          # Output as a vector
  out[[1]] <- gof                                                           # write the gof
  out[[2]] <- qsim                                                          # Write the simulated data
  print(gof)
  # Mandatory names for the elements of the output
  names(out)[1:nelements] <- c("GoF", "sim")                                # Must be the names
  return(out)                                                               # Return the output from the function
} # 'Dynatopmodel End' end



### Function for temporal dynatop function 
# This function calls the spatail data for FR and only runs the temporal series to create a time series
# It is important to go ahead and run the Spatial parts of Dynamictopmodel beforehand
runPSOCalibrationTopmodel <- function(param.values,               # Parameter vector containing values of dynatopmodel parameters
                                      inner.timesteps,            # Number of inner timesteps
                                      rains,                      # Timeseries of precipitation at dt to run the model (continuous data)
                                      PETs,                       # Timeseries of PET at dt needed to run the model (continuous data)
                                      obss,                       # Timeseries of observed Q data to compare the simulated results 
                                      discs,                      # Spatial discretization (weighting matrix) and 'groups' parameter matrix 
                                      RoutingTables=NULL,              # Routing table to discretize stream network
                                      dates.dfs                   # Important dates for the period
) {   
  
  # Read in the time step and number of parameters
  ntt= inner.timesteps                                            # Number of inner timesteps needed to run the numerical model for the baseflow storage
  param = param.values                                            # Parameters to run the model
  
  # Read in the precip, PET, and qobs 
  rain.in <- as.xts(rains)                                        # The previously specified precip data for the calibration period. Additionally, it has been disaggregated.
  pe.in <- as.xts(PETs)                                           # The previously specified PET data for the calibration period. It has already been disaggregated.
  qobs.in <- as.xts(obss)                                         # The previously specified Q.obs
  
  # Read in the existing spatial information
  disc.in <- discs                                                # Load in the discretization data from the Spatial component of the model
  groups <- disc.in$groups                                        # Groups is a data frame within the discretization list that contains parameter values for the HRUs
  RoutingTable.in <- NULL                                # This is the routing table generated from the Spatial component of the model related to the length function of the stream network
  weights <- disc.in$weights                                      # Weights is the weighting matrix showing the percentage of an HRU that flows into another HRU (or the stream)
  
  # Assign the randomly generated parameter to a variable
  v.of.in <- param[1]                                             # The overland flow velocity
  m.in <- exp(param[2])                                                # m - the flow recession parameter
  srz.max.in <- param[3]                                          # The maximum saturated root zone parameter
  srz.0.in <- param[4]                                            # The initial saturation in the root zone
  v.chan.in <- param[5]                                           # The parameter for the velocity in the channel
  natlog.T0.in <- param[6]                                        # The maximum saturation
  sd.max.in <- param[7]                                           # The maximum saturation deficit parameter at which point baseflow is disconnected
  td.in <- param[8]                                               # A timing parameter
  mann.n.in <- param[9]                                           # The manning n
  S0.in <- param[10]                                              # Slope parameter
  CD.in <- param[11]                                              # CD parameter - unused
  K0.in <- param[12]                                              # Initial conductivity parameter - unused
  
  m1.in <- exp(param[13])                                              # m for soil 1
  m2.in <- exp(param[14])                                              # m for soil 2
  m3.in <- exp(param[15])                                              # m for soil 3
  
  # Assign the variables to the groups data frame
  groups$m[groups$soils==1] <- m1.in          #m.in                                      # Assignment to 'groups' for the m parameter           
  groups$m[groups$soils==2] <- m2.in          #m.in
  groups$m[groups$soils==3] <- m3.in          #m.in
  groups$td <-  td.in                                             # For the td parameter
  groups$ln_t0 <- natlog.T0.in                                    # For the T0 parameter
  groups$srz_max <- srz.max.in                                    # For the Srz max parameter
  groups$srz0 <- srz.0.in                                         # For the srz0 parameter
  groups$vchan <- v.chan.in                                       # For the V chan parameter
  groups$vof <- v.of.in                                           # For the V overland flow parameter
  groups$k0 <- K0.in                                              # For the k0 parameter
  groups$CD <- CD.in                                              # For the CD parameter
  groups$sd_max <- sd.max.in                                      # For the sd max parameter
  groups$mann.n <- mann.n.in                                      # For the manning n parameter
  groups$S0 <- S0.in                                              # For the s0 parameter
  
  
  # Run Dynamic Topmodel 
  run <- run.dtm(groups=groups,                                   # Set the groups matrix with the parameter values 
                 weights=weights,                                 # The weights matrix containing the percentage of HRUs that flow to other HRUs
                 rain=rain.in,                                    # Read in the precipitation timeseries
                 pe=pe.in,                                        # Read in the PET timeseries
                 qobs=qobs.in,                                    # Read in the observed Q timeseries
                 qt0=0.0001,                                      # The inital flow assumed to initiate the model
                 routing=RoutingTable.in,                         # Read in the routing table
                 graphics.show=F,                                 # Turn of graphics
                 ntt=ntt,                                         # The inner time steps for the routing model
                 max.q=3.5
                 # The max q to plot
  )                                 # This function Runs Dynatopmodel given the parameter values, weighting matrix, precipitation, PET, and observed stream flow for a given watershed.
  
  
  # Ensure Same Length of TS and Calculate GOF 
  
  q.sim <- zoo(window(run$qsim,start=dates.dfs$warmup.initial,    # Reduce the simulation to the calibration period - NOTE due to routing the simulation
                      end=calibration.final))                         # will last longer than the input period 
  
  
  names(q.sim) <- 'sim'                                           # Name the output
  
  out.model <- list(sim=q.sim,run=run)
  ######NOTE I MAY NEED A DATES IN THE QSIM TO SAVE AS A ZOO OBJECT
  return(out.model)                                                     # Return the list from the function
}


Forcing.input <- function(dt=24,dates.df=dates.df,watershed.area=watershed.area,...) {
  
  # Read in the data
  data.directory <- paste('C:/Users/david/OneDrive/Desktop/EPA/',        # Define the directory string
                          'EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/',                        # I'm breaking up the paste for readability
                          '3 DATA/2-PROJECT1-ROBINSON-FOREST/',         # This isn't the working directory where the PSO items will be stored
                          '3-UK-ROBINSONFOREST-DATA/TIME_SERIES/',      # This is just the WD for reading in input files
                          sep='')  
  setwd(data.directory)                                                 # Set the working directory
  
  Met.data.in <- read_csv('WeatherData_15min/CWS_Rain_clean_hourly.csv')                                    # Read in the data with forcing info
  Q.data.in <- read_csv('Streamflow_15 min/Qobs2000_2015.csv')                              # Read in the observed Q data for falling rock (from Chris/Science base)
  PET.data.in <- read_csv('WATER_pub_2013May2_2.csv')                   # Read in the PETdata
  Temp.data.in <- read_csv('Air_Temp.csv')                              # Read in the temp data
  Day.data.in <- read_csv('FR_precip.csv')                              # Read in the daylight hours
  
  #if (precip.in == 'FR_precip.csv') {                                   # This code will read in data from one of three TS; prioritizing the gage at the top of the basin (could alternatively prioritize an average of the two gages at the top and bottom of the basin), then goes to the bottom, then to the average from TW
  #  Met.data.in <- Met.data.in %>% 
  #    mutate(precip_mm = ifelse(is.na(Met.data.in$precip_mm_FR_H),
  #                              ifelse(is.na(Met.data.in$precip_mm_FR_B),
  #                                     Met.data.in$precip_mm_TW,
  #                                     Met.data.in$precip_mm_FR_B),
  #                              Met.data.in$precip_mm_FR_H))
  #}
  
  
  
  # Convert the dates from character to posixct using lubridate 
  dates.met <- ymd_hms(Met.data.in$fullDate)                # Hr data                  # Converts met date data from character to posixct
  dates.q <- ymd_hms(Q.data.in$roundTime)                   # 15 min data                   # Converts q date data from charcter to posixct
  dates.PET <- mdy_hm(PET.data.in$date2)                    # Daily data
  dates.temp <- mdy_hm(Temp.data.in$date2)                  # Daily data - also have 15 min/hourly temp data if needed 
  dates.day <- mdy_hm(Day.data.in$date2)                    # Daily
  
  # Convert to y-m-d format
  dates.2.temp <- strftime(dates.temp, format="%Y-%m-%d")
  dates.2.day <- strftime(dates.day, format="%Y-%m-%d")
  
  # Convert the data to zoo objecets 
  P.zoo <- zoo(Met.data.in$precip.mm,dates.met)                           # Precip as zoo
  PET.zoo <- zoo(PET.data.in$pet_mm,dates.PET)                            # PET as zoo
  Qobs.zoo <- zoo(Q.data.in$ApproxStreamflow_cfs,dates.q)                 # Mean Qobs as zoo
  #Q.max <- Q.data.in$maxdaily_cfs*.0283168466592                          # Max Q as cms
  #Qobs.max.zoo <- zoo(Q.max,dates.q)                                      # Max Qobs as zoo
  Temp.min.zoo <- zoo((Temp.data.in$DailyMin_degF-32)*(5/9),dates.temp)   # in daily temp
  Temp.max.zoo <- zoo((Temp.data.in$DailyMax_degF-32)*(5/9),dates.temp)   # Max daily temp
  day.zoo <- zoo(Day.data.in$daylight_hr,dates.day)                       # Daylight hours
  
  # Window to desired time
  Tmax=window(Temp.max.zoo,start=dates.day[1],
              end=dates.day[length(dates.day)])
  Tmin=window(Temp.min.zoo,start=dates.day[1],
              end=dates.day[length(dates.day)])
  
  for (i in 1:length(Tmax)) {
    if (is.na(Tmax[i])) {
      Tmax[i]=Tmax[i-1]
    }
    
    if (is.na(Tmin[i])) {
      Tmin[i]=Tmin[i-1]
    }
  }
  
  # Calculate the PET
  PET.data <- list(Tmax=Tmax,
                   Tmin=Tmin,
                   n=day.zoo)  
  
  PET.new <- Hamon(data=PET.data)
  
  # Convert the time series form mm/day to m/hr - the requisite of dynatopmodel
  P.mph.full <- P.zoo*(1/1000)                                            # Manipulate the precip data (m/hr) from Chris Barton is previously in mm/hr
  PET.mph.full <- PET.zoo*(1/1000)*(1/24)                                 # Manipulate the PET data (m/h) from ScienceBase
  Qobs.mph.full <- Qobs.zoo*60*60*(1/watershed.area)*(1/35.3147)          # Manipulate the Q data such that it is in the proper units (m/hr) - originally in cfs
  PET.mph.full.new <- PET.new*(1/1000)*(1/24)
  
  # Disaggregate the time series to the desired time step - NOTE: this is just applying an equal dissgregation method to run the model. There are other disggregation techniques that could be applied (e.g., see the tempdisagg function) - ALSO - we have the max flow of the day that perhaps could help with this
  if (dt != 24) {
    dt=dt
    P.calib <- zoo(aggregate_xts(as.xts(P.mph.full),dt=dt))                          # Disaggregate P. Note - in aggregate_obs it appears that whether is.rate is TRUE or FALSE, the output is the same. Also - the dates aren't retained
    PET.calib <- zoo(aggregate_xts(as.xts(PET.mph.full),dt=dt))                      # Disaggregate PET    
    Qobs.calib <- zoo(aggregate_xts(as.xts(Qobs.mph.full),dt=dt))                    # Disaggregate the Qobs
    PET.new.calib <- zoo(aggregate_xts(as.xts(PET.mph.full.new),dt=dt))
    # NOTE: Need to implement the tempdisagg function and compare
  }
  
  # Select the simulation using the input dates 
  P.calib <- window(P.calib,start=dates.df$warmup.initial,             # Reduce the P time series to the desired start
                    end=dates.df$calibration.final)                       # and end dates
  PET.calib <- window(PET.calib,start=dates.df$warmup.initial,         # Reduce the PET time series to the desired start
                      end=dates.df$calibration.final)                     # and end dates 
  Qobs.calib <- window(Qobs.calib,start=dates.df$warmup.initial,       # Reduce the Qobs time series to the desired start 
                       end=dates.df$calibration.final)                    # and end dates
  PET.new.calib <- window(PET.new.calib,start=dates.df$warmup.initial,       # Reduce the Qobs time series to the desired start 
                          end=dates.df$calibration.final)
  
  
  
  # Write the time series to a data frame to be returned. 
  out.input.data <- list(P=P.calib,PET=PET.new.calib,Qobs=Qobs.calib)         # Combine the time series together
  return(out.input.data)                                                  # Return the variable
  
  # Optionally plot the time series 
  #par(mfrow=c(3,1))                                                      # Generate space for three plots
  #plot(P.calib, xlab = 'Date Time', 
  #     col = 'lightblue', main='Precipitation', ylab ='P (m/hr)')        # Plot the precipitation
  #grid()                                                                 # Include a grid on the plot
  #plot(PET.calib, xlab = 'Date Time', col = 'darkgreen',
  #     main = 'Potential Evapotransipiration', ylab ='PET (m/hr)')       # Plot the PET
  #grid()                                                                 # Include a grid on the plot
  #plot(Qobs.calib, xlab = 'Datet Time', col = 'blue',
  #     main = 'Observed Streamflow', ylab = 'Q (mm/hr)')                 # Plot the observed Q
  #grid()                                                                 # Include a grid on the plot
}


Hamon <- function(data=PET.data) {
  Ta <- (data$Tmax + data$Tmin)/2
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax/(data$Tmax + 237.3))
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin/(data$Tmin + 237.3))
  vas <- (vs_Tmax + vs_Tmin)/2
  ET_Hamon.Daily <- 0.55 * 25.4 * (data$n/12)^2 * (216.7 * vas * 10/(Ta + 273.3))/100
  ET.Daily <- ET_Hamon.Daily
  
  return(ET.Daily)
}





### Post processing and routing of Q after calibration...

explicit.routing.instant <- function(FR.Spatial,run,model.timestep) {
  
  ############# READ IN REACH INPUTS AND CLEAN MATRICES #################
  num.reach <- nrow(FR.Spatial$explicit.ChanTable)
  explicit.disc <- FR.Spatial$explicit.disc                           # Read in the explicit disc
  explicit.weights.full <- FR.Spatial$explicit.disc$weights
  rownames(explicit.weights.full)[rownames(explicit.weights.full) == 'R25.6666666666667'] <- 'R26' # Replace R25.6 with R26
  colnames(explicit.weights.full)[colnames(explicit.weights.full) == 'R25.6666666666667'] <- 'R26' # Replace R25.6 with R26
  explicit.weights.full[,17] <- explicit.weights.full[,17]+                   # Correct the matrix for the 17 half reach... I'm assuming that 17.5 is just an extension of 17??
    explicit.weights.full[,18]
  explicit.weights.full[17,] <- explicit.weights.full[17,]+                   # Correct the matrix for the 17 half reach... I'm assuming that 17.5 is just an extension of 17??
    explicit.weights.full[18,]
  explicit.weights.full <- explicit.weights.full[,-18]   
  explicit.weights.full <- explicit.weights.full[-18,] 
  explicit.groups <- FR.Spatial$explicit.disc$groups[-34,]
  groups <- FR.Spatial$disc$groups
  
  # Initialize input matrices
  lumped.chan.inputs <- matrix(0,ncol= ncol(run$fluxes$qbf), 
                               nrow=nrow(run$fluxes$qbf))
  lumped.chan.overland <- matrix(0,ncol= ncol(run$fluxes$qbf), 
                                 nrow=nrow(run$fluxes$qbf))
  lumped.qin <- run$fluxes$qin
  explicit.chan.inputs <- matrix(0,ncol=length(explicit.groups$area),
                                 nrow=nrow(run$fluxes$qbf))
  explicit.chan.inputs.sa <- matrix(0,ncol=length(explicit.groups$area),
                                    nrow=nrow(run$fluxes$qbf))
  explicit.chan.inputs.qbf <- matrix(0,ncol=length(explicit.groups$area),
                                     nrow=nrow(run$fluxes$qbf))
  explicit.chan.inputs.qbf.sa <- matrix(0,ncol=length(explicit.groups$area),
                                        nrow=nrow(run$fluxes$qbf))
  
  evap <- data.frame(matrix(0,ncol=1,nrow=nrow(run$fluxes$ae)))
  names(evap) <- 'ae'
  a.chan.explicit <- explicit.groups$area
  a.chan <- groups$area
  catch.area <- sum(a.chan)
  Q_out <- matrix(0,ncol=1,nrow=nrow(run$fluxes$qin))
  AE_out <- matrix(0,ncol=1,nrow=nrow(run$fluxes$qin))
  total_in <- matrix(0,ncol=1,nrow=nrow(run$fluxes$qin))
  total_out <- matrix(0,ncol=1,nrow=nrow(run$fluxes$qin))
  sd.gain <- matrix(0,ncol=1,nrow=nrow(run$fluxes$qin))
  in.rain.init <- as.vector(run$fluxes$rain[,1]) 
  in.ae.init <- as.vector(run$fluxes$ae[,1])    # NOTE ASSUMES THAT WE ARE USING A SINGULAR AE AND P FOR THE ENTIRE CATCHMENT 
  in.qbf.init <- data.frame(run$fluxes$qbf)
  
  for (iter in 1:nrow(run$fluxes$qbf)) {
    in.rain <- in.rain.init[iter]                               # Rain inputs
    in.ae <- in.ae.init[iter]                                   # Actual Evapotranspiration inputs
    
    in.qbf <- as.numeric(in.qbf.init[iter,])
    
    #evap[iter,'ae'] <- in.ae[2]
    pex <- in.rain-in.ae
    pex[pex<0] <- 0
    
    flows.explicit.qbf <- rep(0,nrow(explicit.weights.full))  
    
    
    flows.explicit.qbf[34:length(flows.explicit.qbf)] <- in.qbf[3:length(in.qbf)] 
    
    explicit.chan.inputs.qbf[iter,] <- as.vector((flows.explicit.qbf*explicit.groups$area) %*% explicit.weights.full)
    explicit.chan.inputs.qbf[iter,] <- explicit.chan.inputs.qbf[iter,] +
      (pex[1])*a.chan.explicit
    
    correction <- as.numeric(run$Qsim[iter])/sum(explicit.chan.inputs.qbf[iter,1:33])
    
    explicit.chan.inputs.qbf[iter,] <- explicit.chan.inputs.qbf[iter,]*correction
    
    
  }
  #plot(x=as.numeric(run$Qsim[1:length(run$fluxes$qin[,1])]),y=as.numeric(rowSums(explicit.chan.inputs.qbf[,1:33])))
  
  if(length(run$Qsim)==12502) {
    perc.diff.Qsim.explicit <- as.vector((run$Qsim-rowSums(explicit.chan.inputs.qbf[1:length(run$Qsim),1:33]))/run$Qsim*100)
    
    check.balances <- data.frame(explicit=rowSums(explicit.chan.inputs.qbf[1:length(run$Qsim),1:33]),
                                 Qsim=run$Qsim,
                                 diff.explicit.Qsim = rowSums(explicit.chan.inputs.qbf[1:length(run$Qsim),1:33]) -
                                   run$Qsim,
                                 qin=lumped.qin[1:length(run$Qsim),1],
                                 diff.qin.Qsim = lumped.qin[1:length(run$Qsim),1]-
                                   run$Qsim,
                                 perc.diff.Qsim.explicit=perc.diff.Qsim.explicit)
    
  } else {
    perc.diff.Qsim.explicit <- as.vector((run$Qsim[1:nrow(run$fluxes$qbf)]-rowSums(explicit.chan.inputs.qbf[,1:33]))/run$Qsim[1:nrow(run$fluxes$qbf)]*100)
    check.balances <- data.frame(explicit=rowSums(explicit.chan.inputs.qbf[,1:33]),
                                 Qsim=run$Qsim[1:nrow(run$fluxes$qbf)],
                                 diff.explicit.Qsim = rowSums(explicit.chan.inputs.qbf[,1:33]) -
                                   run$Qsim[1:nrow(run$fluxes$qbf)],
                                 qin=lumped.qin[,1],
                                 diff.qin.Qsim = lumped.qin[,1]-
                                   run$Qsim[1:nrow(run$fluxes$qbf)],
                                 perc.diff.Qsim.explicit=perc.diff.Qsim.explicit)
  }
  
  
  ############# ROUTE WATER FROM REACHES ############
  
  ####### Routing with no dispersion/lag #######
  explicit.reach.table <- FR.Spatial$explicit.ChanTable           # Write the table of reach information 
  reach.name <- explicit.reach.table$link_no                      # Write the reach names/link numbers
  explicit.chan.inputs <- data.frame(explicit.chan.inputs.qbf)
  reach.name <- c(reach.name,groups$id[3:length(groups$id)])      # Need to check if there are duplicate names
  duplicate.reaches.hru <- sum(duplicated(reach.name))
  # if duplcicated.reaches.hru > 0, then output a warning
  
  names(explicit.chan.inputs) <- reach.name                       # Assign reach names to columns in the inflow matrix
  if (length(run$Qsim) == 12502) {
    explicit.chan.inputs <- explicit.chan.inputs[,1:33]
  } else {
    explicit.chan.inputs <- explicit.chan.inputs[,1:33]
  }
  
  routed.flow.instant <- data.frame(matrix(nrow=nrow(explicit.chan.inputs),  # Initialize a dataframe of the same size of the channel inputs
                                           ncol=ncol(explicit.chan.inputs))) # Continued - for instant routed 
  names(routed.flow.instant) <- names(explicit.chan.inputs)                   # Assigne reach names to columsn in the isntant routed  flow
  reach.name <- reach.name[1:33]
  
  routed.flow.kinematic <- data.frame(matrix(nrow=nrow(explicit.chan.inputs),  # Initialize a dataframe of the same size of the channel inputs
                                             ncol=ncol(explicit.chan.inputs))) # Continued - for instant routed 
  names(routed.flow.kinematic) <- names(explicit.chan.inputs)      
  
  # FIRST ORDER REACHES - calculate for the first order reaches first
  for (reach.1 in reach.name[explicit.reach.table[,5]==1]) {      # For each of the reaches where in the explicit reach table the stream order is 1
    reach.char.1 <- as.character(reach.1)                         # Write the reach as a character (for referenceing purposes)
    routed.flow.instant[,reach.char.1] <-                         # Write to the routed.flow.instant data frame with the reach character as the column name
      explicit.chan.inputs[,as.character(reach.1)]                # Since it's a first order reach, the flow is just from the inputs
    routed.flow.kinematic[,reach.char.1] <- explicit.chan.inputs[,as.character(reach.1)]
  }
  
  # SECOND ORDER REACHES - calculate for the second order reaches
  second.order <- rev(reach.name[explicit.reach.table[,5]==2])    # Reverse the order because the lower 2nd order reaches might contain other 2nd order reaches
  for (reach.2 in second.order) {                                 # For each of the reaches where in the explicit reach table the stream order is 2 
    reach.char.2 <- as.character(reach.2)                         # Write the reach as a character for referencing
    routed.flow.instant[,reach.char.2] <-                         # Write to the routed.flow.instant data frame with reach
      explicit.chan.inputs[,as.character(reach.2)] +                                                             # Inflow from the uplands
      routed.flow.instant[,as.character(explicit.reach.table[,3][explicit.reach.table$link_no==reach.2])] + # Inflow from the first reach
      routed.flow.instant[,as.character(explicit.reach.table[,4][explicit.reach.table$link_no==reach.2])]   # Inflow from the second reach
    
    routed.flow.kinematic[,reach.char.2] <- kinematicWave(lateral.inflow=explicit.chan.inputs[,as.character(reach.2)],
                                                          inflow.reach.i=routed.flow.instant[,as.character(explicit.reach.table[,3][explicit.reach.table$link_no==reach.2])] + # Inflow from the first reach
                                                            routed.flow.instant[,as.character(explicit.reach.table[,4][explicit.reach.table$link_no==reach.2])],
                                                          width=explicit.reach.table$width_m2[explicit.reach.table$link_no==reach.2],
                                                          channel.length=explicit.reach.table$stream_length_ft[explicit.reach.table$link_no==reach.2]*0.3048,
                                                          slope=explicit.reach.table$stream_slope_ftpft[explicit.reach.table$link_no==reach.2],
                                                          manning.n=FR.Spatial$disc$groups$mann.n[1],
                                                          timestep=model.timestep*60*60)
    
  }
  
  #THIRD ORDER REACHES
  third.order <- rev(reach.name[explicit.reach.table[,5]==3])     # Reverse the order because the lower 3rd order reaches might contain other 3rd order reaches
  for (reach.3 in third.order) {                                  # For each of the reaches where in the explicit reach table the stream order is 3 
    reach.char.3 <- as.character(reach.3)                         # Write the reach as a character for referencing
    routed.flow.instant[,reach.char.3] <-                         # Write to the routed.flow.instant data frame with reach
      explicit.chan.inputs[,as.character(reach.3)] +                                                              # Inflow from the uplands
      routed.flow.instant[,as.character(explicit.reach.table[,3][explicit.reach.table$link_no==reach.3])] + # Inflow from the first reach
      routed.flow.instant[,as.character(explicit.reach.table[,4][explicit.reach.table$link_no==reach.3])]   # Inflow from the second reach
    
    routed.flow.kinematic[,reach.char.3] <- kinematicWave(lateral.inflow=explicit.chan.inputs[,as.character(reach.3)],
                                                          inflow.reach.i=routed.flow.instant[,as.character(explicit.reach.table[,3][explicit.reach.table$link_no==reach.3])] + # Inflow from the first reach
                                                            routed.flow.instant[,as.character(explicit.reach.table[,4][explicit.reach.table$link_no==reach.3])],
                                                          width=explicit.reach.table$width_m2[explicit.reach.table$link_no==reach.3],
                                                          channel.length=explicit.reach.table$stream_length_ft[explicit.reach.table$link_no==reach.3]*0.3048,
                                                          slope=explicit.reach.table$stream_slope_ftpft[explicit.reach.table$link_no==reach.3],
                                                          manning.n=FR.Spatial$disc$groups$mann.n[1],
                                                          timestep=model.timestep*60*60)
  }
  
  #FOURTH ORDER REACHES
  fourth.order <- rev(reach.name[explicit.reach.table[,5]==4])    # Reverse the order because the lower 4th order reaches might contain other 4th order reaches
  for (reach.4 in fourth.order) {                                 # For each of the reaches where in the explicit reach table the stream order is 4
    reach.char.4 <- as.character(reach.4)                         # Write the reach as a character for referencing
    routed.flow.instant[,reach.char.4] <-                         # Write to the routed.flow.instant data frame with reach
      explicit.chan.inputs[,as.character(reach.4)] +                                                             # Inflow from the uplands
      routed.flow.instant[,as.character(explicit.reach.table[,3][explicit.reach.table$link_no==reach.4])] + # Inflow from the first reach
      routed.flow.instant[,as.character(explicit.reach.table[,4][explicit.reach.table$link_no==reach.4])]   # Inflow from the second reach
    
    routed.flow.kinematic[,reach.char.4] <- kinematicWave(lateral.inflow=explicit.chan.inputs[,as.character(reach.4)],
                                                          inflow.reach.i=routed.flow.instant[,as.character(explicit.reach.table[,3][explicit.reach.table$link_no==reach.4])] + # Inflow from the first reach
                                                            routed.flow.instant[,as.character(explicit.reach.table[,4][explicit.reach.table$link_no==reach.4])],
                                                          width=explicit.reach.table$width_m2[explicit.reach.table$link_no==reach.4],
                                                          channel.length=explicit.reach.table$stream_length_ft[explicit.reach.table$link_no==reach.4]*0.3048,
                                                          slope=explicit.reach.table$stream_slope_ftpft[explicit.reach.table$link_no==reach.4],
                                                          manning.n=FR.Spatial$disc$groups$mann.n[1],
                                                          timestep=model.timestep*60*60)
  }
  
  if (length(run$Qsim)==12502) {
    Qout.lumped=as.vector(run$Qsim)/(3600)
    Qout.explicit.routed=routed.flow.instant[1:length(run$Qsim),1]/(3600)
  } else {
    Qout.lumped=as.vector(run$Qsim[1:nrow(run$fluxes$qbf)])/(3600)
    Qout.explicit.routed=routed.flow.instant[,1]/(3600)
  }
  
  
  if (length(run$Qsim) ==12502) {
    percent.diff=((run$Qsim-routed.flow.instant[1:length(run$Qsim),1])/run$Qsim)*100
    check.balances.new <- data.frame(Qout=run$Qsim,
                                     routed.flow.instant=routed.flow.instant[1:length(run$Qsim),1],
                                     percent.diff=percent.diff)
  } else {
    percent.diff=((run$Qsim[1:nrow(run$fluxes$qbf)]-routed.flow.instant[,1])/run$Qsim[1:nrow(run$fluxes$qbf)])*100
    check.balances.new <- data.frame(Qout=run$Qsim[1:nrow(run$fluxes$qbf)],
                                     routed.flow.instant=routed.flow.instant[,1],
                                     percent.diff=percent.diff)
  }
  
  check.balances.newest <- list(check.balances=check.balances.new,min.perc.diff=min(percent.diff),max.perc.diff=max(percent.diff),avg.perc.diff=mean(percent.diff))
  
  # Current units of routed.flow.instant are m^3/hr... convert to mm/s
  routed.flow.instant.mm_s <- sweep(routed.flow.instant,2,explicit.reach.table$US_area_m2,FUN='/')/3600*1000
  
  routing.out <- list(check.balances.new=check.balances.newest, routed.flow.instant.m3_hr=routed.flow.instant,
                      routed.flow.instant.mm_s=routed.flow.instant.mm_s,
                      routed.flow.kinematic.m3_hr=routed.flow.kinematic)
  
  
  #View(check.balances)
  #View(check.balances.new)
  #max(check.balances.new$percent.diff)
  # min(check.balances.new$percent.diff)
  return(routing.out)
  
  
}



### Function for spatial dynatopmodel components
DynatopSpatialFunctionExplicitReaches <- function(){
  
  # set WD
  setwd('C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/1 GIS ANALYSIS/')
  
  # Read in 1.5 m DEM, Flow length (unused), and soil Raster
  fr1mDEM <- raster('Rasters/FR1meterDEM.tif')                                                # Read in the 1 m FR DEM
  
  # NOTE AS OF 7/6/2021 - there has been an update to the raster package such and now the raster function appears to no longer be reading the crs (coordinate system) of the raster correctly. 
  # The Crs is: +proj=lcc +lat_0=36.3333333333333 +lon_0=-85.75 +lat_1=37.0833333333333 +lat_2=38.6666666666667 +x_0=1500000 +y_0=999999.9998984 +datum=NAD83 +units=us-ft +no_defs 
  
  crs(fr1mDEM) <- '+proj=lcc +lat_0=36.3333333333333 +lon_0=-85.75 +lat_1=37.0833333333333 +lat_2=38.6666666666667 +x_0=1500000 +y_0=999999.9998984 +datum=NAD83 +units=us-ft +no_defs'
  
  #fr1mFlowLen <- raster('Rasters/frFlowLen1m2.tif')                                           # Read in the flow length raster
  FRSoils1m <- raster('Rasters/FRSoils1mMASK.tif')                                            # Read in the soils raster
  crs(FRSoils1m) <- '+proj=lcc +lat_0=36.3333333333333 +lon_0=-85.75 +lat_1=37.0833333333333 +lat_2=38.6666666666667 +x_0=1500000 +y_0=999999.9998984 +datum=NAD83 +units=us-ft +no_defs'
  
  # Read in stream network
  FR_drn <- shapefile('Shapes/FR4000StreamNet.shp')
  bkf.width <- 12.16*(FR_drn$USAreaM2/1000^2*0.386102)^0.42*0.3048                            # estimate the width of the channel using the regional curves from Ashal Berry's thesis. Note this isn't perfect bc this equation is being extrapolated to smaller channels. 
  FR_drn_table <- as.data.frame(cbind('link_no'= FR_drn$LINKNO,'ds_link'= FR_drn$DSLINKNO,
                                      'us_link_1' = FR_drn$USLINKNO1,'us_link_2'=FR_drn$USLINKNO2,
                                      'stream_order'= FR_drn$strmOrder,'stream_length_ft'=FR_drn$Length,
                                      'stream_drop_ft'=FR_drn$strmDrop,'stream_slope_ftpft'=FR_drn$Slope,
                                      'stream_straightline_length_ft'=FR_drn$StraightL,'US_area_m2'=FR_drn$USAreaM2,
                                      'width_m2' = bkf.width))
  FR_drn_table$stream_length_ft[1] = 248.4
  
  # Create FR layers and TWI (labeled as atb in the model)
  layers <- build_layers(fr1mDEM,fill.sinks=TRUE)                                     # Create TWI and US contributing area raster layers
  
  # Build channels for FR
  chans.explicit <- build_chans2(dem=fr1mDEM, drn=FR_drn, buffer=2,                          # This automatically creates channels for the catchment
                                 chan.width=3.1,single.chan=F)                                    # Using this function, the spatially explicit reach information is retained. 
  chans.lumped <- build_chans(dem=fr1mDEM, drn=FR_drn, buffer=2,                            # This automatically creates channels for the catchment
                              chan.width=3.1)                                                  # using the drn from TauDEM. 
  
  
  # Discretize the layers for FR
  layers2 <- addLayer(fr1mDEM, layers$atb, FRSoils1m)                                 # Input layers to be used to discretize the watershed (different than above)
  FR.disc <- discretise(layers2, cuts=c(atb=10,FRSoils1mMASK=3),                      # Discretize the watershed into HRUs based on the TWI and soils data
                        chans=chans.lumped, area.thresh=0.05/100)                          # Note, it is imperative to specify an area.threshold >= 0.01, otherwise the model will not run
  FR.disc.explicit <- discretise(layers2, cuts=c(atb=10,FRSoils1mMASK=3),                    # Discretize the watershed into HRUs based on the TWI and soils data
                                 chans=chans.explicit, area.thresh=0.05/100)                           # Note, it is imperative to specify an area.threshold >= 0.01, otherwise the model will not run
  
  correct.ratio <- FR.disc$weights[3:nrow(FR.disc$weights),1]/
    rowSums(FR.disc.explicit$weights[,1:34])[35:nrow(FR.disc.explicit$weights)]
  # View(correct.ratio)
  if (rowSums(FR.disc.explicit$weights[,1:34])[35:nrow(FR.disc.explicit$weights)] == 0 &&
      FR.disc$weights[3:nrow(FR.disc$weights),1] > 0) {
    warning('There is zero flow in an explicit reach but positive flow in a lumped')
  }
  
  for (row in 1:(nrow(FR.disc$weights)-2)) {                                                 # here I am assuming that the ratios from HRUs into channels is correct for the lumped. So I am making sure that the ratios from HRUs in the explicit discretization matches that of the first
    FR.disc.explicit$weights[(34+row),1:34] <- FR.disc.explicit$weights[(34+row),1:34]*
      correct.ratio[row]
    if (is.na(FR.disc.explicit$weights[(34+row),1:34])) {
      FR.disc.explicit$weights[(34+row),1:34] <- 0
    }
  }
  
  FR.disc$groups$area <- FR.disc$groups$area/10.7639                                         # Convert from ft2 to m2
  FR.disc.explicit$groups$area <- FR.disc.explicit$groups$area/10.7639                       # convert from ft2 to m2
  
  #View(rowSums(FR.disc.explicit$weights[,1:34]))
  #View(FR.disc$weights[,1])
  #View(FR.disc$weights)
  #View(FR.disc.explicit$weights)
  #correct.ratio.new <- rowSums(FR.disc.explicit$weights[,1:34])[35:nrow(FR.disc.explicit$weights)]/FR.disc$weights[3:nrow(FR.disc$weights),1]
  #View(correct.ratio.new)                                                         # This should be one across the board, unless there is zero flow in from the explicit and the lumped
  #View(FR.disc$groups)
  #View(FR.disc.explicit$groups)
  #View(FR.disc$weights[,1])
  
  #setwd('C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/FR_Spatial_Function_23FEB2021')
  #write.csv(FR.disc$weights,'FRWeights_LumpedReach.csv')
  #write.csv(FR.disc.explicit$weights,'FRWeights_ExplicitReach.csv')
  #write.csv(FR.disc$groups, 'FRgroups_lumped.csv')
  #write.csv(FR.disc.explicit$groups,'FRgroups_explicit.csv')
  #writeRaster(chans.lumped,'chans.tif',overwrite=TRUE)
  #writeRaster(chans.explicit,'chans2.tif',overwrite=TRUE)
  #writeRaster(FR.disc$hru,'HRU.tif',overwrite=TRUE)
  #writeRaster(FR.disc.explicit$hru,'HRUs2.tif',overwrite=TRUE)
  
  # Build routing table for FR
  FR.RoutingTable <- build_routing_table(dem=fr1mDEM,breaks=30,chans=chans.lumped)   # The routing table for the time delay fucntion using network width assumptions
  FR.RoutingTable.explicit <- build_routing_table(dem=fr1mDEM,breaks=5,chans=chans.explicit)
  # write the spatial data to a list
  FRSpatial <- list("DEM"=fr1mDEM,"Soils"=FRSoils1m,"layers"=layers,"disc"=FR.disc,"RoutingTable"=FR.RoutingTable,"explicit.disc"=FR.disc.explicit,"explicit.RoutingTable"=FR.RoutingTable.explicit,"explicit.ChanTable"=FR_drn_table)
  return(FRSpatial)
}

build_chans2 <- function(dem,
                         drn,
                         chan.width=4,
                         atb=NULL,
                         buffer=10,
                         atb.thresh=0.8,
                         single.chan=TRUE)
{
  if(!is.null(drn))
  {
    # build the reach multi band raster and save
    message("Building raster for channel(s)...")
    reaches <- build.reach.raster(dem, drn, buffer=buffer,
                                  chan.width=chan.width)
    
    # single channel desired but multi-reach river network supplied
    #reaches[[1]] <- reaches[[1]]>0
    
    # ensure DEM and reaches have same extent ( processing can leave the former larger)
    reaches <- crop(reaches, dem)
    # Estimate proportion of river cells occupied by channel
    #  prop <- min(chan.width/xres(dem), 1)
    
  }
  else if(!is.null(atb))
  {
    # using the TWI to identify the channel
    reaches <- atb > atb.thresh*max(atb[], na.rm=TRUE)
    # Estimate proportion of river cells occupied by channel
    prop <- min(chan.width/xres(atb), 1)
    reaches[which(reaches[]==0)] <- NA
    cellprops <- reaches*prop
    reaches <- addLayer(reaches, cellprops)
  }
  else
  {
    stop("Provide vector channel data or raster layer to locate channels")
  }
  
  names(reaches) <- c("chans", "chanprops")
  return(reaches)
}

#############################################################
# raster of reach locations and cell proportion occupied
#############################################################
build.reach.raster <- function(dem, drn,
                               buffer=10,
                               
                               chan.width=4)
{
  # determnine the number of distinct reaches in the reach vector and compare with the chan.width vector
  # if greater then apply the width to all reaches
  drn <- as(drn, "SpatialLines")
  nchan <- length(drn@lines)
  if(length(chan.width) < nchan)
  {
    chan.width <- rep(chan.width, length.out=nchan)
    # have only supplied a single width across catchment
  }
  chan.width <- as.vector(chan.width)
  dem <- dem+0  # now in memory to prevent disk access
  
  # for reaches wider than cell size, expand the network so that adjacent cells are occupied
  buffer.width <- chan.width/2
  ichan <- which(chan.width>xres(dem))
  if(length(ichan)>0)
  {
    buffer.width[ichan] <- 0
  }
  # create buffer with specified width (on either side)
  drn <- rgeos::gBuffer(drn, width=buffer.width, byid=TRUE)
  
  reaches <- dem
  reaches[] <- NA
  # one row per channel
  message("Extracting reach cells...")
  rl <- raster::extract(dem, drn, cellnumbers=TRUE)
  
  # firt column the cell number occupied by the channel, second the reach index, third the width of this channel
  rlw <- lapply(1:length(rl),
                function(i)
                {
                  cbind(rl[[i]][,1], i, chan.width[i])
                  
                })
  
  rlw <- do.call(rbind, rlw)
  
  reaches[rlw[,1]] <- rlw[,2]
  
  # proportions of channel cells occupied by channel
  props <- reaches
  props[rlw[,1]] <-   pmin(rlw[,3]/xres(dem), 1)
  reaches <- addLayer(reaches,props)
  
  # calculate the cell proportions
  #	prop <- min(chan.width/xres(dem), 1)
  #	cellprops <- reaches*min(chan.width/xres(dem), 1)
  
  # add a layer with proportions of cell occuppied by channel. estimated by
  # proportion of cell size to channel width, probably close enough
  #	reaches <- addLayer(reaches, (reaches>0)*prop)
  
  names(reaches)=c("chan", "chanprop")
  
  return(reaches)
}

## Create a function to route water in a reach with the kinematic wave routing method
# Note: this assumes the kinematic approximation to the saint venant equations, no lateral inflow (e.g., flow comes in at the start of a reach), and a wide rectuangular channel 

kinematicWave <- function(lateral.inflow,inflow.reach.i,width,channel.length,slope,manning.n,timestep) {
  
  ## Explanation of arguments
  # lateral.inflow = vector of lateral inflow to reach j during timestep i (m^3/hr)
  # inflow.reach.i = vector of inflow from the previous reach (m^3/hr)
  # width = width of the reach (m)
  # length = legnth of the reach (m)
  # slope = slope of the channel (m/m)
  # manning.n = manning n of the channel (unitless) 
  # timestep = timestep of the channel (s)
  
  # calculate the alpha and beta components from Manning's equation 
  alpha.kinematic <- (manning.n*width^(2/3)/(1.49*(slope)^0.5)) ^ (3/5)
  beta.kinematic <- 0.6
  
  # Convert the lateral inflow (m^3/s) to flow per unit width 
  lateral.inflow <- lateral.inflow/channel.length
  
  # Convert the flow from m^3/hr to m^3/s
  inflow.reach.i <- inflow.reach.i/60/60
  lateral.inflow <- lateral.inflow/60/60
  
  # Initilize the outflow matrix
  Q.outflow <- matrix(nrow=length(inflow.reach.i),ncol=1)
  # Implment the numerical method 
  for (i in 1:length(inflow.reach.i)) {
    if (i==1) {
      Q.outflow[i] = inflow.reach.i[i]
    } else {
      Q.inflow.numerator <- timestep/channel.length*inflow.reach.i[i]+alpha.kinematic*beta.kinematic*((Q.outflow[i-1]+inflow.reach.i[i])/2)^(beta.kinematic-1)*Q.outflow[i-1]
      Q.lateral.numerator <- timestep*(lateral.inflow[i]+lateral.inflow[i-1])/2
      Q.inflow.denominator <- timestep/channel.length+alpha.kinematic*beta.kinematic*((Q.outflow[i-1]+inflow.reach.i[i])/2)^(beta.kinematic-1)
      
      Q.outflow[i] <- (Q.inflow.numerator+Q.lateral.numerator)/Q.inflow.denominator
    }
  }
  
  # Convert back to m^3/hr
  Q.outflow <- Q.outflow*60*60 
  return(Q.outflow)
}