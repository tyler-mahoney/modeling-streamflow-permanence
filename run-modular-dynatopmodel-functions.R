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
  
  save.out <- paste0(model.drty, '/fluxes_stores/')
  fluxes <- simLump$run$fluxes
  names.fluxes <- names(fluxes)
  for (flux in 1:length(fluxes)) {
    flux.name <- names.fluxes[flux]
    for (iter in 1:600) {
      # note: 600 is arbitrary, If we run the model > 600 times there will be an issue. 
      out.flux.file <-
        paste0(save.out, 'flux', flux.name, iter, '.csv')
      if (!file.exists(out.flux.file)) {
        write.csv(fluxes[[flux]], out.flux.file)
        break
      }
    }
  }
  stores <- simLump$run$storages
  names.stores <- names(stores)
  for (store in 1:length(stores)) {
    store.name <- names.stores[store]
    for (iter in 1:600) {
      out.store.file <- paste0(save.out, 'store', store.name, iter, '.csv')
      if (!file.exists(out.store.file)) {
        write.csv(stores[[store]], out.store.file)
        break
      }
    }
  }
  
  for (iter in 1:600) {
    out.q.file <- paste0(save.out, 'qsim_specific', iter, '.csv')
    if (!file.exists(out.q.file)) {
      #write.csv(qsim,out.q.file)
      write.csv(simLump$run$Qsim[1:length(simLump$run$fluxes$qbf[, 1])],
                paste0(save.out, 'Q_sim_big', iter, '.csv'))
      write.csv(simLump$run$qsim[1:length(simLump$run$fluxes$qbf[, 1])],
                paste0(save.out, 'qsim_specific', iter, '.csv'))
      write.csv(gof, paste0(save.out, 'gof_log_kge', iter, '.csv'))
      
      break
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


explicit.routing.instant <- function(read.spatial,explicit.reach.table,run,model.timestep) {
  
  ############# READ IN REACH INPUTS AND CLEAN MATRICES #################
  num.reach <- nrow(explicit.reach.table)
  explicit.disc <- read.spatial$explicit.disc                          # Read in the explicit disc
  explicit.weights.full <- read.spatial$explicit.disc$weights
  explicit.groups <- explicit.disc$groups
  groups <- read.spatial$disc$groups
  
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
    
    
    flows.explicit.qbf[(num.reach+1):length(flows.explicit.qbf)] <- in.qbf[2:length(in.qbf)] 
    
    explicit.chan.inputs.qbf[iter,] <- as.vector((flows.explicit.qbf*explicit.groups$area) %*% explicit.weights.full)
    explicit.chan.inputs.qbf[iter,] <- explicit.chan.inputs.qbf[iter,] +
      (pex[1])*a.chan.explicit
    
    correction <- as.numeric(run$Qsim[iter])/sum(explicit.chan.inputs.qbf[iter,1:num.reach])
    
    explicit.chan.inputs.qbf[iter,] <- explicit.chan.inputs.qbf[iter,]*correction
    
    
  }
  
  # Check the water balance
  #  perc.diff.Qsim.explicit <- as.vector((run$Qsim-rowSums(explicit.chan.inputs.qbf[1:length(run$Qsim),1:num.reach]))/run$Qsim*100)
  
  #    check.balances <- data.frame(explicit=rowSums(explicit.chan.inputs.qbf[1:length(run$Qsim),1:num.reach]),
  #                                 Qsim=run$Qsim,
  #                                 diff.explicit.Qsim = rowSums(explicit.chan.inputs.qbf[1:length(run$Qsim),1:num.reach]) -
  #                                   run$Qsim,
  #                                 qin=lumped.qin[1:length(run$Qsim),1],
  #                                 diff.qin.Qsim = lumped.qin[1:length(run$Qsim),1]-
  #                                   run$Qsim,
  #                                 perc.diff.Qsim.explicit=perc.diff.Qsim.explicit)
  
  
  ############# ROUTE WATER FROM REACHES ############
  
  ####### Routing with no dispersion/lag #######
  reach.name <- paste0('R',explicit.reach.table$link_no)                      # Write the reach names/link numbers
  explicit.chan.inputs <- data.frame(explicit.chan.inputs.qbf)
  reach.name <- c(reach.name,paste0('HRU',groups$id[2:length(groups$id)]))      # Need to check if there are duplicate names
  duplicate.reaches.hru <- sum(duplicated(reach.name))
  # if duplcicated.reaches.hru > 0, then output a warning
  if (duplicate.reaches.hru > 0) {warning('reach names are duplicated. this will likely cause errors.')}
  
  names(explicit.chan.inputs) <- reach.name                       # Assign reach names to columns in the inflow matrix
  explicit.chan.inputs <- explicit.chan.inputs[,1:num.reach]
  
  routed.flow.instant <- data.frame(matrix(nrow=nrow(explicit.chan.inputs),  # Initialize a dataframe of the same size of the channel inputs
                                           ncol=ncol(explicit.chan.inputs))) # Continued - for instant routed  (m^3/hr)
  
  names(routed.flow.instant) <- names(explicit.chan.inputs)                   # Assigne reach names to columsn in the isntant routed  flow
  reach.name <- reach.name[1:num.reach]
  
  routed.flow.kinematic <- data.frame(matrix(nrow=nrow(explicit.chan.inputs),  # Initialize a dataframe of the same size of the channel inputs
                                             ncol=ncol(explicit.chan.inputs))) # Continued - for instant routed 
  names(routed.flow.kinematic) <- names(explicit.chan.inputs)      
  
  # We want to run the routing algorithm such that we start at the most distal reaches and work our way down stream. 
  
  # We need to calculate this for two types of first-order reaches: "true" first order reaches, which don't have any upstream reach and 
  # "downstream" first order reaches - those which have contribution from an upstream reach which is still a first order reach.
  
  # This account for breaking a reach up into multiple reaches. Note: if we were to do this for second- or third-order reaches, we'd need to do the samething for those reaches. 
  
  first.order.true <- rev(reach.name[which(explicit.reach.table$`stream_order`==1 & explicit.reach.table$`us_link_1`==-1 & explicit.reach.table$`us_link_2`==-1)]) 
  first.order.down.stream <- rev(reach.name[which(explicit.reach.table$`stream_order`==1 & explicit.reach.table$`us_link_1`!=-1 & explicit.reach.table$`us_link_2`==-1)]) 
  
  # TRUE FIRST ORDER REACHES - calculate for the first order reaches first at the highest point in the watershed
  for (reach.1 in first.order.true) {      # For each of the reaches where in the explicit reach table the stream order is 1
    reach.char.1 <- as.character(reach.1)                         # Write the reach as a character (for referenceing purposes)
    name.no.r <- as.integer(substr(reach.1,2,nchar(reach.1)))     # Convert the name to an integer because that's how it is in the explicit.reach.table
    
    routed.flow.instant[,reach.char.1] <-                         # Write to the routed.flow.instant data frame with the reach character as the column name
      explicit.chan.inputs[,as.character(reach.1)]                # Since it's a first order reach, the flow is just from the inputs
    routed.flow.kinematic[,reach.char.1] <- explicit.chan.inputs[,as.character(reach.1)]
  }
  
  # DOWNSTREAM FIRST ORDER REACHES - calculate for the first order reaches first DOWNSTREAM OF THOSE at the highest point in the watershed
  for (reach.1 in first.order.down.stream) {      # For each of the reaches where in the explicit reach table the stream order is 1 and there is a SINGLE contributing 1st order reach
    reach.char.1 <- as.character(reach.1)                         # Write the reach as a character (for referenceing purposes)
    name.no.r <- as.integer(substr(reach.1,2,nchar(reach.1)))     # Convert the name to an integer because that's how it is in the explicit.reach.table
    routed.flow.instant[,reach.char.1] <-                         # Write to the routed.flow.instant data frame with the reach character as the column name
      explicit.chan.inputs[,as.character(reach.1)] +                # Since it's a first order reach, the flow is just from the inputs
      routed.flow.instant[,paste0('R',as.character(explicit.reach.table$us_link_1[explicit.reach.table$link_no==name.no.r]))]
    
    routed.flow.kinematic[,reach.char.1] <- as.numeric(kinematicWave(lateral.inflow=explicit.chan.inputs[,as.character(reach.1)],                                                        inflow.reach.i=routed.flow.instant[,paste0('R',explicit.reach.table$us_link_1[explicit.reach.table$link_no==name.no.r])],
                                                                     width=explicit.reach.table$width_m2[explicit.reach.table$link_no==name.no.r],
                                                                     channel.length=explicit.reach.table$stream_length_ft[explicit.reach.table$link_no==name.no.r]*0.3048,
                                                                     slope=explicit.reach.table$stream_slope_ftpft[explicit.reach.table$link_no==name.no.r],
                                                                     manning.n=read.spatial$disc$groups$mann.n[1],
                                                                     timestep=model.timestep*60*60))
    
  }
  
  # SECOND ORDER REACHES - calculate for the second order reaches
  second.order <- rev(reach.name[explicit.reach.table$`stream_order`==2])    # Reverse the order because the lower 2nd order reaches might contain other 2nd order reaches
  for (reach.2 in second.order) {                                 # For each of the reaches where in the explicit reach table the stream order is 2 
    reach.char.2 <- as.character(reach.2)                         # Write the reach as a character for referencing
    name.no.r <- as.integer(substr(reach.2,2,nchar(reach.2)))     # Convert the name to an integer because that's how it is in the explicit.reach.table
    routed.flow.instant[,reach.char.2] <-                         # Write to the routed.flow.instant data frame with reach
      explicit.chan.inputs[,as.character(reach.2)] +                                                             # Inflow from the uplands
      routed.flow.instant[,paste0('R',as.character(explicit.reach.table$us_link_1[explicit.reach.table$link_no==name.no.r]))] + # Inflow from the first reach
      routed.flow.instant[,paste0('R',as.character(explicit.reach.table$us_link_2[explicit.reach.table$link_no==name.no.r]))]   # Inflow from the second reach
    
    routed.flow.kinematic[,reach.char.2] <- as.numeric(kinematicWave(lateral.inflow=explicit.chan.inputs[,as.character(reach.2)],                                                        inflow.reach.i=routed.flow.instant[,paste0('R',explicit.reach.table$us_link_1[explicit.reach.table$link_no==name.no.r])] + # Inflow from the first reach
                                                                       routed.flow.instant[,paste0('R',explicit.reach.table$us_link_2[explicit.reach.table$link_no==name.no.r])],
                                                                     width=explicit.reach.table$width_m2[explicit.reach.table$link_no==name.no.r],
                                                                     channel.length=explicit.reach.table$stream_length_ft[explicit.reach.table$link_no==name.no.r]*0.3048,
                                                                     slope=explicit.reach.table$stream_slope_ftpft[explicit.reach.table$link_no==name.no.r],
                                                                     manning.n=read.spatial$disc$groups$mann.n[1],
                                                                     timestep=model.timestep*60*60))
    
  }
  
  #THIRD ORDER REACHES
  third.order <- rev(reach.name[explicit.reach.table$`stream_order`==3])     # Reverse the order because the lower 3rd order reaches might contain other 3rd order reaches
  for (reach.3 in third.order) {                                  # For each of the reaches where in the explicit reach table the stream order is 3 
    reach.char.3 <- as.character(reach.3)                         # Write the reach as a character for referencing
    name.no.r <- as.integer(substr(reach.3,2,nchar(reach.3)))
    routed.flow.instant[,reach.char.3] <-                         # Write to the routed.flow.instant data frame with reach
      explicit.chan.inputs[,as.character(reach.3)] +                                                              # Inflow from the uplands
      routed.flow.instant[,paste0('R',as.character(explicit.reach.table$us_link_1[explicit.reach.table$link_no==name.no.r]))] + # Inflow from the first reach
      routed.flow.instant[,paste0('R',as.character(explicit.reach.table$us_link_2[explicit.reach.table$link_no==name.no.r]))]   # Inflow from the second reach
    
    routed.flow.kinematic[,reach.char.3] <- as.numeric(kinematicWave(lateral.inflow=explicit.chan.inputs[,as.character(reach.3)],                                                        inflow.reach.i=routed.flow.instant[,paste0('R',explicit.reach.table$us_link_1[explicit.reach.table$link_no==name.no.r])] + # Inflow from the first reach
                                                                       routed.flow.instant[,paste0('R',explicit.reach.table$us_link_2[explicit.reach.table$link_no==name.no.r])],
                                                                     width=explicit.reach.table$width_m2[explicit.reach.table$link_no==name.no.r],
                                                                     channel.length=explicit.reach.table$stream_length_ft[explicit.reach.table$link_no==name.no.r]*0.3048,
                                                                     slope=explicit.reach.table$stream_slope_ftpft[explicit.reach.table$link_no==name.no.r],
                                                                     manning.n=read.spatial$disc$groups$mann.n[1],
                                                                     timestep=model.timestep*60*60))
  }
  
  #FOURTH ORDER REACHES
  fourth.order <- rev(reach.name[explicit.reach.table$`stream_order`==4])    # Reverse the order because the lower 4th order reaches might contain other 4th order reaches
  for (reach.4 in fourth.order) {                                 # For each of the reaches where in the explicit reach table the stream order is 4
    reach.char.4 <- as.character(reach.4)                         # Write the reach as a character for referencing
    name.no.r <- as.integer(substr(reach.4,2,nchar(reach.4)))
    try(routed.flow.instant[,reach.char.4] <-                         # Write to the routed.flow.instant data frame with reach
          explicit.chan.inputs[,as.character(reach.4)] +                                                             # Inflow from the uplands
          routed.flow.instant[,paste0('R',as.character(explicit.reach.table$us_link_1[explicit.reach.table$link_no==name.no.r]))] + # Inflow from the first reach
          routed.flow.instant[,paste0('R',as.character(explicit.reach.table$us_link_2[explicit.reach.table$link_no==name.no.r]))] )  # Inflow from the second reach
    
    try(routed.flow.kinematic[,reach.char.4] <- as.numeric(kinematicWave(lateral.inflow=explicit.chan.inputs[,as.character(reach.4)],                                                        inflow.reach.i=routed.flow.instant[,paste0('R',explicit.reach.table$us_link_1[explicit.reach.table$link_no==name.no.r])] + # Inflow from the first reach
                                                                           routed.flow.instant[,paste0('R',explicit.reach.table$us_link_2[explicit.reach.table$link_no==name.no.r])],
                                                                         width=explicit.reach.table$width_m2[explicit.reach.table$link_no==name.no.r],
                                                                         channel.length=explicit.reach.table$stream_length_ft[explicit.reach.table$link_no==name.no.r]*0.3048,
                                                                         slope=explicit.reach.table$stream_slope_ftpft[explicit.reach.table$link_no==name.no.r],
                                                                         manning.n=read.spatial$disc$groups$mann.n[1],
                                                                         timestep=model.timestep*60*60)))
  }
  
  # Comparing m^3/s
  Qout.lumped=as.vector(run$Qsim[1:nrow(run$fluxes$qbf)])/(3600)
  Qout.explicit.routed=routed.flow.instant[,1]/(3600)
  
  percent.diff=((run$Qsim[1:nrow(run$fluxes$qbf)]-routed.flow.instant[,1])/run$Qsim[1:nrow(run$fluxes$qbf)])*100
  check.balances.new <- data.frame(Qout=run$Qsim[1:nrow(run$fluxes$qbf)],
                                   routed.flow.instant=routed.flow.instant[,1],
                                   percent.diff=percent.diff)
  
  check.balances.newest <- list(check.balances=check.balances.new,min.perc.diff=min(percent.diff),max.perc.diff=max(percent.diff),avg.perc.diff=mean(percent.diff))
  
  # Current units of routed.flow.instant are m^3/hr... convert to mm/s
  routed.flow.instant.mm_s <- sweep(routed.flow.instant,2,explicit.reach.table$US_area_m2,FUN='/')/3600*1000
  routed.flow.kinematic.mm_s <- sweep(routed.flow.kinematic,2,explicit.reach.table$US_area_m2,FUN='/')/3600*1000
  
  routing.out <- list(check.balances.new=check.balances.newest, routed.flow.instant.m3_hr=routed.flow.instant,
                      routed.flow.instant.mm_s=routed.flow.instant.mm_s,
                      routed.flow.kinematic.m3_hr=routed.flow.kinematic,
                      routed.flow.kinematic.mm_s=routed.flow.kinematic.mm_s)
  
  
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
  # inflow.reach.i = vector of inflow from the upstream reach (m^3/hr)
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

 
# Function 5: Fucntion to evaluate each behavioral parameterization and simulate streamflow permanence  
Headwater.evaluation.dynamic <- function(model.drty,qin.file, read.spatial, TS.logger.clean.all, TS.zoo.all, reach.identifiers, model.timestep, time.sim, ...) {
  # Read in explicit HRU results for fluxes and storages from the behavioral sets
  # NOTE WE MUST SET THE DIRECTORY TO THE PROPER CALIBRATION TEST FOLDER, AS OF 7/6/2021 IT IS CalibrationTest13
  fluxes.stores.dir <- paste0(model.drty,'/fluxes_stores/')
  kge.run <- read.csv(paste0(fluxes.stores.dir,'gof_log_kge',qin.file,'.csv'))
  
  
  
  
  
  if(kge.run[2]>=0.3){
    
    
    #setwd('C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/Calibration/CalibrationTest13/fluxes_stores')
    qin.run <- read.csv(paste0(fluxes.stores.dir,'fluxqin',qin.file,'.csv'))                           # Read qin 
    qin.run <- xts(qin.run[,-1],time.sim)                                            # Convert to xts
    qbf.run <- read.csv(paste0(fluxes.stores.dir,'fluxqbf',qin.file,'.csv'))                           # Read qbf
    qbf.run <- xts(qbf.run[,-1],time.sim)                                            # Convert to xts
    ae.run <- read.csv(paste0(fluxes.stores.dir,'fluxae',qin.file,'.csv'))                             # Read ae
    ae.run <- xts(ae.run[,-1],time.sim)                                              # Convert to xts
    rain.run <- read.csv(paste0(fluxes.stores.dir,'fluxrain',qin.file,'.csv'))                         # Read rain
    rain.run <- xts(rain.run[,-1],time.sim)                                          # Convert to xts
    qof.run <- read.csv(paste0(fluxes.stores.dir,'fluxqof',qin.file,'.csv'))                           # Read qof
    qof.run <- xts(qof.run[,-1],time.sim)                                            # Convert to xts
    Qsim.run <- read.csv(paste0(fluxes.stores.dir,'Q_sim_big',qin.file,'.csv'))                             # Read Qsim
    Qsim.run <- xts(Qsim.run[,-1],time.sim)                                          # Convert to xts
    compile.fluxes <- list(qin = qin.run, qbf = qbf.run, ae = ae.run,                # Compile fluxes 
                           rain = rain.run, qof = qof.run)                           # Into list
    compile.run <- list(fluxes=compile.fluxes, Qsim = Qsim.run)                      # Compile fluxes and Qsim into list
    # Run routing function -- note TOPMODEL puts out data in m/hr for reach-specific flow in runoff
    explicit.routing.run <- explicit.routing.instant(read.spatial,                     # Read in read.spatial and the compiled run
                                                     explicit.reach.table = explicit.reach.table,
                                                     run=compile.run,
                                                     model.timestep = model.timestep)# Run the explicit routing code for the simulation
    
    Q <- data.frame(explicit.routing.run$routed.flow.kinematic.mm_s)                   # Assign the routed flow instant (mm/s) to a variables
    r.names <- names(Q)                                                              # Get the reach names
    
    # Convert Q from mm/s to mm/day
    Q <- Q*60*60*24                                                                  # This is important for some of the empirical calculations we'll later need to do
    
    # Calculate the fdc of each reach 
    Flow.duration.reaches <- matrix(nrow = length(Q[,1]),ncol = length(r.names))   # initialize the flow duration matrix
    Q.in.rank <- matrix(nrow = length(Q[,1]),ncol = length(r.names))               # Initialize the Q in rank matrix
    colnames(Flow.duration.reaches) <- r.names                                       # Assign names to the flow duration matrix
    colnames(Q.in.rank) <- r.names                                                   # Assign names to the Q in rank matrix
    
    # Run the FDC code for the reaches
    for (r.name in 1:length(r.names)) {                                              # For each of the reaches
      
      Q.in <- Q[,r.name]                                                           # Assing the Q in for a reach
      order <- explicit.reach.table[r.name,6]                                      # Get the stream order
      Q.in.sort <- data.frame(flow.asc=sort(Q.in,decreasing =T))                   # Rank/sort the Q in
      rank <- 1:length(Q.in)                                                       # Get the rank
      Q.in.sort$rank <- rank                                                       # assign the rank for the sorted data dataframe
      Q.in.sort$Prob.Q.in <- 100*(Q.in.sort$rank/(length(Q.in)+1))                 # Calculate the probablility
      Flow.duration.reaches[,paste0('R',as.character(explicit.reach.table$link_no[r.name]))] <- Q.in.sort$Prob.Q.in                                # Record the probability 
      Q.in.rank[,paste0('R',as.character(explicit.reach.table$link_no[r.name]))] <- Q.in.sort$flow.asc                                             # Record the ascending flow
      
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
    max.threshold <- 20                                                             # This was the maximum transmissivity scaling exponent measured from the Prancevic and Kirchner 2018 Paper
    gamma.range <- seq(from=log(min.threshold), to = log(max.threshold),                       # Creates a sequency of 1000 thresholds for which we'll calculate the % correct
                       by = (log(max.threshold)-log(min.threshold))/1000)                 # Continued 
    
    # Format the data 
    names.loggers <- names(TS.logger.clean.all)
    
    # Initialize the lists for each of the loggers. 
    logger.TS <- list()                                                              # Initialize the logger Timeseries list
    logger.timeseries <- list()
    sim.obs.difference <- list()
    confusion.matrix.thresh <- list()
    out.thresh.df <- list()
    reach.slope <- data.frame(matrix(nrow=length(names.loggers)))
    reach.contributing.area <- data.frame(matrix(nrow=length(names.loggers)))
    
    # Calculate the respective error for each flow thresh for each sensor
    residual.on.off <- matrix(nrow=length(gamma.range),ncol=length(names.loggers))
    
    for (name.logger in 1:length(names.loggers)) {
      # Format the time series properly
      #time.sim <- time(result.optim.run$run$fluxes$qin)                                # Get the time stamp for the qsim
      time.sim <- data.frame(Date=time.sim)                                            # Set the col name to 'Date'
      logger.TS[[name.logger]] <- data.frame(state=TS.zoo.all[[name.logger]],
                                             Date=TS.logger.clean.all[[name.logger]]$Date)                  # Get the logger time series
      
      # Join the logger.TS to the time.sim 
      logger.timeseries[[name.logger]] <- left_join(time.sim,logger.TS[[name.logger]],by='Date')                                    # Join the logger data to the time series
      logger.timeseries[[name.logger]]$state <- ifelse(logger.timeseries[[name.logger]]$state==-1,NA,logger.timeseries[[name.logger]]$state)        # Set NA values
      logger.timeseries[[name.logger]] <- logger.timeseries[[name.logger]][complete.cases(logger.timeseries[[name.logger]]),]                      # Remove rows with NA values
      
      # initialize the matrix that records the difference in simulated and observed
      sim.obs.difference[[name.logger]] <- data.frame(
        matrix(nrow=length(logger.timeseries[[name.logger]]$Date),ncol=length(gamma.range)+1))
      
      # Get the contributing area and slope for the reach 
      reach.name.dynamic <- which(                                                 # Determine which row the reach name of the logger belongs to
        read.spatial$explicit.ChanTable$link_no==                                    # Compare the list of reach numbers to the current logger being evaluated
          as.numeric(reach.identifiers[name.logger]))                  # Define the reach name as a numeric
      reach.slope[name.logger,1] <-                                                  # Determine the slope of the reach (ft/ft) or (m/m)
        read.spatial$explicit.ChanTable$stream_slope_ftpft[reach.name.dynamic]       # Using the which function from above
      reach.contributing.area[name.logger,1] <-                                      # Determine the contributing area (m^2) using the which function from aboce
        read.spatial$explicit.ChanTable$US_area_m2[reach.name.dynamic]               # Using the function from above
      
      # Convert contributing area to km2
      reach.contributing.area[name.logger,1] <- 
        reach.contributing.area[name.logger,1]*1e-6
      
      
      
      # Run the analysis to simulate the 1/0 time series for the four sensors over the 1000 transmissivity scaling exponents
      for (gamma.in in 1:length(gamma.range)) {                                              
        # Write the simulated 1/0 time series for the average threshold
        
        # Determine the subsurface flow capacity of the valley below the reach. This equation is proposed in Godsey and Kirchner (2018) and Prancevic and Kirchner (2014)
        transmissivity <- exp(gamma.range[gamma.in])
        thresh <- transmissivity*reach.slope[name.logger,1]
        
        #percent.thresh <- Flow.duration.reaches[which(r.names==Reach.identifiers[name.logger])][nearest(Q.in.rank[which(r.names==Reach.identifiers[name.logger])],         # Look at the flow duration for reach 810 (associated with FC4) and 
        #                                                                                                thresh),1]                            # Look at the flow duration for reach 810 (associated with FC4) and 
        explicit.r.names <- names(explicit.routing.run$routed.flow.kinematic.m3_hr)
        
        
        reach.identifiers.second <- paste0('R',reach.identifiers)
        
        sim.on.off <- 
          ifelse(explicit.routing.run$routed.flow.kinematic.m3_hr[which(explicit.r.names==reach.identifiers.second[name.logger])]>thresh,           # USing the iterative threshold
                 1,0)                                                               # Set the values greater than the threshold to one, otherwise set to zero 
        sim.on.off <- data.frame(sim.on.off)                                           # Convert to a data frame
        #time.sim <- time(result.optim.run$run$fluxes$qin)                              # Get the time stamp for the qsim
        sim.on.off$Date <- time.sim$Date                                                    # Set the date of the data frame to the time.sim 
        
        # Compare the logger 1/0 to the simulated 1/0
        sim.obs.match <- left_join(logger.timeseries[[name.logger]],sim.on.off,by='Date')                           # Join the logger time series to the simulated time series by the date
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
        
        
      }
      
      
      residual.on.off[,name.logger] <- confusion.matrix.thresh[[name.logger]][,3]+
        confusion.matrix.thresh[[name.logger]][,4]
      
    }
    Sys.time()-x
    
    
    
    # Calculate the total residual error by summing the error for each sensor
    rowSums.residual.on.off <- data.frame(sum.error=rowSums(residual.on.off),transmissivity=exp(gamma.range))
    
    for (i in 1:length(names.loggers)) { 
      thresh.vals = data.frame(out.thresh.df[[i]]$flow.thresh)
      colnames(thresh.vals) <- paste0(reach.identifiers[i],' thresh')
      rowSums.residual.on.off <- cbind(rowSums.residual.on.off,thresh.vals)
    }
    
    
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
    Q.sub.c.best <- thresh.best*read.spatial$explicit.ChanTable$stream_slope_ftpft
    total.network.on.off <- data.frame(matrix(nrow=nrow(explicit.routing.run$routed.flow.kinematic.m3_hr),ncol=ncol(explicit.routing.run$routed.flow.kinematic.m3_hr)))
    colnames(total.network.on.off) <- r.names
    for (reach in 1:length(Q.sub.c.best)) {
      total.network.on.off[,reach] <- ifelse(explicit.routing.run$routed.flow.kinematic.m3_hr[,reach]>Q.sub.c.best[reach],1,0)
    }
    
    flow.network.thresh <- explicit.routing.run$routed.flow.kinematic.m3_hr
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
    
    percent.on.network <- data.frame(percent.on.network)
    rownames(percent.on.network) <- reach.names.list
    # Compile everything as a list
    out.headwater.evaluation <- list(                                       # Compile a list of everything we want to output
      flow.thresh.best.results=flow.thresh.best.results, # stats.on.off=stats.on.off,
      thresh.best=thresh.best, #best.percent.thresh=best.percent.thresh,
      total.network.on.off=total.network.on.off,
      flow.network.thresh=flow.network.thresh, percent.on.network=percent.on.network, 
      save.out.thresh.best= save.out.thresh.best,
      best.residual.on.off=best.residual.on.off,
      Q.sub.c.best=Q.sub.c.best,
      kge=kge.run,
      qin.file=qin.file,
      Q=Q
    )         
  } else {
    print('run does not meet minimum KGE threshold, we will not include this parameterization in the next stage')
    out.headwater.evaluation <- NA
  }
  
                                                                               # Continued
  
  return(out.headwater.evaluation)                                                 # return the list and exit the function
}
