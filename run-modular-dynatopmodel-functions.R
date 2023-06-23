# Title: run-modular-dynatopmodel-functions
# Description: These functions are used to run the dynamic TOPMODEL

# Function 0: Quiet the model
# Dynamic TOPMODEL autmatically outputs the result for each time step, which is annoying. I'm going to write a function to skip that...
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Function 1: Run dynatopmodel
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
  groups$m[groups$soils==4] <- m1.in          #m.in                                      # HEY! This needs to be dynamic! There could be more than 4 soil types!!! 
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
  run <- quiet(run.dtm(groups=groups,                             # Set the groups matrix with the parameter values 
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
  ))                                 # This function Runs Dynatopmodel (quietly) given the parameter values, weighting matrix, precipitation, PET, and observed stream flow for a given watershed.
  
  
  # Ensure Same Length of TS and Calculate GOF 
  
  q.sim <- zoo(window(run$qsim,start=dates.dfs$warmup.initial,    # Reduce the simulation to the calibration period - NOTE due to routing the simulation
                      end=calibration.final))                         # will last longer than the input period 
  
  
  names(q.sim) <- 'sim'                                           # Name the output
  
  out.model <- list(sim=q.sim,run=run)
  
  return(out.model)                                                     # Return the list from the function
}

# Function 2: Setup the calibration for dynamic TOPMODEL
Dynatophydromod <- function(param.values,                                        # Parameter vector containing values of dynatopmodel parameters
                            inner.timestep,                                      # Number of inner timesteps
                            rain,                                                # Timeseries of precipitation at dt to run the model (continuous data)
                            PET,                                                 # Timeseries of PET at dt needed to run the model (continuous data)
                            obs,                                                 # Timeseries of observed Q data to compare the simulated results 
                            disc,                                                # Spatial discretization (weighting matrix) and 'groups' parameter matrix 
                            RoutingTable=NULL,                                   # Routing table to discretize stream network
                            date,                                                # Optional
                            dates.df,                                            # Date frame of dates used for warm up and calibration periods 
                            model.drty,                                          # need the model directory to save everything to
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
  
  obs <- window(simLump$run$qobs,start=dates.df$calibration.initial,
                end = dates.df$calibration.final)
  qsim <- window(simLump$run$qsim, start = dates.df$calibration.initial,
                 end = dates.df$calibration.final)
  
  n <- length(qsim)                                                         # Removing the warming up period
  gof <- KGE(sim=as.numeric(log(qsim+1e-10)), obs=as.numeric(log(obs+1e-10)), method="2012")                             # Caclculate KGE of the log of the ts. added a small amount to avoid zeros. 
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
        write.csv(simLump$run$Qsim[1:length(simLump$run$fluxes$qbf[,1])],paste0(save.out,'Q_sim_big',iter,'.csv'))
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

# Function 3: Routing of Q from upper reaches to the watershed outlet
# Post processing and routing of Q after calibration...

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


# Function 4: Kinematic Wave Routing Method:
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