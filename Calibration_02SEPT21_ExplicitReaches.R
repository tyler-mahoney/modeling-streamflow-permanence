## Calibration and Uncertainty Analysis for Falling Rock 
# Uses PSO calibration method to optimize the time series

# Insert Change 

# Read in the functions from this script ---- MAKE SURE TO UPDATE THIS AS NEEDED
source('DynamicTopmodel_FallingRock_Functions_02SEPT2021.R')

# Read in Libraries
library(dynatopmodel)                                           # Read in the dynatopmodel package
library(raster)                                                 # Read in the raster package to analyze raster files - NOTE it's necessary to load libraries in this order otherwise the code breaks.... 
library(hydroPSO)                                               # Read in hydroPSO algorithm (particle search optimizations)
library(hydroGOF)                                               # Read in hydroGOF package (goodness of fit)
library(hydroTSM)                                               # Read in hydroTSM package
library(boot)                                                   # Read in boot package
library(zoo)                                                    # Read in zoo package for time series analysis            
library(sensitivity)                                            # Read in the sensitivity analysis package
library(xts)                                                    # Read in the xts package also for time series analysis
library(ggplot2)                                                # Read in ggplot2 package for figure manipulation
library(dplyr)                                                  # Read in dplyr package
library(readr)                                                  # Read in read r package
library(lubridate)                                              # Read in lubridate for date manipulation
library(Evapotranspiration)                                     # Read in the Evapotranspiration packahge 
library(rasterVis)                                              # Read in the rasterVis package for raster visualizations
library(dismo)                                                  # Read in dismo package for k fold cross validation
library(rgl)                                                    # Read in the rgl package for 3D scatter plots
library(matrixStats)                                            # Read in the library for matrixStats
library(gridExtra)                                              # Read in the gridExtra package for GGPLOT2
library(tidyquant)

# REQUIRED FUNCTIONS: SEE functionS AT THE END OF THIS SCRIPT FOR NECESSARY FUNCTIONS FOR CALIBRATION

# Run the Spatial Equations and Load Discretization
FR.Spatial <- DynatopSpatialFunctionExplicitReaches()           # Runs the spatial discretization equations
FR.disc <- FR.Spatial$disc                                      # Load in the spatial discretization (weighting matrix) and 'groups' parameter matrix 
FR.RoutingTable <- FR.Spatial$RoutingTable                      # Load in the routing table
explicit.reach.table <- FR.Spatial$explicit.ChanTable           # Assign the explicit chan Table to a variable with reach bathymetry information

# Read the HRUs and assign unique soil groups NOTE: there has to be a more efficient way to do this, but right now it doesn't matter because this isn't being iterated. 
no.hru <- as.numeric(length(unique(FR.disc$hru)))               # Calculate the number of unique HRUs from the spatial function
j=2                                                             # Iterator for the soils vector    
soils <- array(no.hru+1)                                        # Initialize the soils vector
soils[1] <- 1                                                   # Set the initial row to 1, for the stream network, which is not an HRU
for (i in unique(FR.disc$hru)) {                                # The loop will output the soil type for each HRU
  soil.type <- FR.disc$layers$FRSoils1mMASK[FR.disc$hru==i]     # Identifies the soil type for a corresponding HRU
  soils[j] <- soil.type[1]                                      # Since there are multiple HRU units, the first soil type belonging to the HRU will be read
  j <- j+1                                                      # Increase the iterator
}                                                               # End the loop
FR.disc$groups <- cbind(FR.disc$groups,soils)                   # The soil type is bound to the groups layer

# Specify the Model Time Step (and internal time step)
model.timestep <- 2                                             # Time step of the model (in Hours)
in.timestep <- 1           #2                                     # Specify the number of internal time steps (unitless) - this is for the internal loop to calculate the baseflow

# Specify Calibration Period
calibration.initial <- "2003-11-23 20:00:00"                    # Calibration start       # calibration.initial <- "2010-11-24 12:00:00"
calibration.final <- "2006-09-30 20:00:00"                      # Calibration end         # calibration.final <- "2012-12-01 12:00:00"   

# Specify Warmup Period
warmup.initial <- "2003-05-23 20:00:00"                         # Warmup start            # warmup.initial <- "2010-07-24 12:00:00" 
warmup.final <- "2003-11-24 00:00:00"                           # Warmup end              # warmup.final <- "2010-11-24 12:00:00" 
#NOTE: FOR SOME REASON AT THIS TIME STEP THE WINDOW FUNCTION IS ADDING TWO TIME STEPS 

# Concatinate the dates
dates.df <- data.frame(calibration.initial, calibration.final,  # Concatinate the dates into a data frame 
                       warmup.initial, warmup.final)            # Continued
names.dates <- c('calibration.initial','calibration.final',     # Create names for the df
                 'warmup.initial','warmup.final')               # Continued
names(dates.df) <- names.dates                                  # Specify the date names 

# Specify the watershed area
watershed.area <- 970966.97                                     # Specify the watershed area (m^2) for Falling Rock

# Load in Met and Qobs time series from csv, manipulate to the desired time step, and write the function output to variables
# NOTE: we would also want to put validation initial/final dates.
# If the validation and calibration periods are not one after the other, need a second warm up period. 
input.timeseries <- Forcing.input(dt=model.timestep,            # Run the forcing input function
                                  dates.df=dates.df,            # Start and end dates and watershed area
                                  watershed.area=watershed.area,# Output is a list of three zoo objects for the rain, PET, and Qobs time series
                                  )                             # Input
rain.calib <- input.timeseries$P[1:(length(input.timeseries$P))]                     # Writing the P TS to a variable
PET.calib <- input.timeseries$PET[1:(length(input.timeseries$PET))]                  # Writing the PET TS to a variable
Q.obs.calib <- input.timeseries$Qobs[1:(length(input.timeseries$Qobs))]              # Writing the Qobs TS to a variable
dates.cal.warm <- time(rain.calib)                                                   # Extract dates for the timeseries

# Set up the model parameters; Read in min and max values (these values should be based on the literature).
v.of.min <- 10; v.of.max <- 150                                 # Overland flow velocity (m/hr) min and max
m.min <- -9.908; m.max <- -2.813                                # Form of exponential decline in conductivity (m) min and max
srz.max.min <-  .01; srz.max.max <- .75                         # Max root zone storage (m) min and max
srz.0.min <- 0.5; srz.0.max <- 1                                # Initial root zone storage (fraction) min and max
v.chan.min <- 500; v.chan.max <- 7000                           # Channel routing velocity (m/hr) min and max
natlog.T0.min <- 3; natlog.T0.max <- 16                         # Lat saturated transmissivity (m^2/hr) min and max
sd.max.min <- 0.2; sd.max.max <- 0.8                            # Max effective deficit of sat zone (m) min and max
td.min <- 0.01; td.max <- 100                                   # Unsat zone time delay (hr/m) min and max
mann.n.min <- 0.01; mann.n.max <- .15                           # Manning's (unitless) n min and max
S0.min <- 0.01; S0.max <- 0.3                                   # Nominal slope (fraction) min and max
CD.min <-  0.01; CD.max <- 0.5                                  # Capilary drive (unused) min and max
k0.min <-  0.1; k0.max <- 1000                                  # Initial saturated hydraulic conductivity (m/hr; unused) min and max
m1.min <- -9.908; m1.max <- -2.813                              # Assign range for the first m     
m2.min <- -9.908; m2.max <- -2.813                              # Assign range for the second m
m3.min <- -9.908; m3.max <- -2.813                              # Assign range for the third m 

# Concatinate parameter minmum, maximum, and names. Assign names to vectors
lower <- c(v.of.min, m.min, srz.max.min, srz.0.min, v.chan.min, # Vector of lower parametr values
           natlog.T0.min, sd.max.min, td.min,                   # Will be used to generate parameter sets 
           mann.n.min, S0.min, CD.min, k0.min,                  # Continued
           m1.min, m2.min, m3.min)                              # During calibration

upper <- c(v.of.max, m.max, srz.max.max, srz.0.max, v.chan.max, # Vector of upper parameter values
           natlog.T0.max, sd.max.max, td.max,                   # Will be used to generate parmeter sets 
           mann.n.max, S0.max, CD.max, k0.max,                  # Continued 
           m1.max, m2.max, m3.max)                              # During Calibration 

name.params <- c('v.of', 'm','srz.max','srz.0',                 # Create names of all parameters
                 'v.chan','natlog.T0','sd.max',                 # Included in Dynatopmodel
                 'td','mann.n','S0','CD','k0',                  # Continued
                 'm1','m2','m3')                                # To be implemented within the uncertainty and calibration
names(lower) <- name.params                                     # Specify names of columns for the lower part of the parameter space
names(upper) <- name.params                                     # Specify names of columns for the upper part of the parameter space

params <- lower+(upper-lower)/2                                 # Unused - but for a test run of the model

params.best.FC4 <- c(0.978,-6.75, 0.679,0.981,500,5.76,.606,1.58,0.0376,.181,.136,.0693,-9.91,-3.22,-2.81)

# Set the Working Directory (to run the calibration) - MAKE SURE TO RENAME THE FOLDER TO THE APPROPRIATE CALIBRATION TEST AS OF 6/27/2022 IT IS CalibrationOutput - a folder stored in GitHub
model.drty <- paste0(getwd(),'/CalibrationOutput',sep='')
#setwd(model.drty)                                                                  # Set the working directory again 

# Run a test of the model to make sure that things look good before running the PSO calibration 
result.test <- runPSOCalibrationTopmodel(param.values=params,                      # TEST - Run the model with the optimal parameter set
                                         inner.timesteps = in.timestep,            # input the inner.timesteps
                                         rains = rain.calib, PETs=PET.calib,       # All forcing components are the same
                                         obss=Q.obs.calib, discs=FR.disc,          # Spatial information comes from the dynatop spatial function
                                         RoutingTables=NULL,                       # See above
                                         dates.dfs = dates.df)                     # Calibration end date input
result.test.FC4 <- runPSOCalibrationTopmodel(param.values=params.best.FC4,                      # TEST - Run the model with the optimal parameter set
                                             inner.timesteps = in.timestep,            # input the inner.timesteps
                                             rains = rain.calib, PETs=PET.calib,       # All forcing components are the same
                                             obss=Q.obs.calib, discs=FR.disc,          # Spatial information comes from the dynatop spatial function
                                             RoutingTables=NULL,                       # See above
                                             dates.dfs = dates.df)                     # Calibration end date input


qsim.test <- window(result.test$sim, start = dates.df$calibration.initial,         # Window test sim
                    end = calibration.final)                                       # Continued
obs.test <- window(Q.obs.calib,start=dates.df$calibration.initial,                 # Window test obs
                   end = calibration.final)                                        # Continued

if (length(qsim.test)>length(obs.test)){
  qsim.test <- qsim.test[1:length(obs.test)]
} else {
  obs.test <- obs.test[1:length(qsim.test)]
}

gof.logKGE.test <- KGE(sim=as.numeric(log(qsim.test+1e-08)),                       # Test stats - log KGE
                       obs=as.numeric(log(obs.test+1e-08)), method="2012")         # Continued 
gof.KGE.test <- KGE(sim=as.numeric(qsim.test),                                     # Test stats - KGE
                    obs=as.numeric(obs.test), method="2012")                       # Continued 
gof.logNSE.test <- NSE(sim=as.numeric(log(qsim.test+1e-08)),                      # Test stats - log NSE
                       obs=as.numeric(log(obs.test+1e-08)))                       # Continued
gof.NSE.test <- NSE(sim=as.numeric(qsim.test),obs=as.numeric(obs.test))                                  # Test stats - NSE

gof.test <- data.frame(logKGE=gof.logKGE.test,KGE=gof.KGE.test,                    # Record all of the test stats
                       logNSE=gof.logNSE.test,NSE=gof.NSE.test)                    # Continued 
gof.test                                                                           # Print test stats 

result.test.2 <- Dynatophydromod(param.values=params,                              # SECOND TEST Parameter vector containing values of dynatopmodel parameters
                                 inner.timestep=in.timestep,                       # Number of inner timesteps
                                 rain=rain.calib,                                  # Timeseries of precipitation at dt to run the model (continuous data)
                                 PET=PET.calib,                                    # Timeseries of PET at dt needed to run the model (continuous data)
                                 obs=Q.obs.calib,                                  # Timeseries of observed Q data to compare the simulated results 
                                 disc=FR.disc,                                     # Spatial discretization (weighting matrix) and 'groups' parameter matrix 
                                 RoutingTable=NULL,                                # Routing table to discretize stream network
                                 date=dates.cal.warm,                              # Optional
                                 dates.df=dates.df,
                                 model.drty=model.drty                             # Date frame of dates used for warm up and calibration periods 
)


# Set up the general settings
Figures.drty.out <- paste0(model.drty, "/Figures")                                 # Create directroy to store the figures for the calibration
if (!file.exists(Figures.drty.out)) dir.create(Figures.drty.out,                   # Creates a figure directory
                                               recursive=TRUE)                     # If the output directory selected to store the figures does not exists, it is created:
RoutingTables <- NULL                                                              # Note: as of now since FR is so small, I am assuming that the routing table is null, meaning that there is no attenuation really. This assumption is okay for now given the small size of the catchment.
RoutingTable <- NULL                                                               # Note: as of now since FR is so small, I am assuming that the routing table is null, meaning that there is no attenuation really. This assumption is okay for now given the small size of the catchment.

# Run the PSO calibration procedure: Note we were having issues running the model with the model.FUN.args as the method to read the arguments. 
# NOTE: Now we run the hydroPSO algorithm with the model arguments as 'global' variables.
out <- hydroPSO(fn="hydromodInR",                               # This is an optional, but necessary specification that the calibration approach is a 'hydromod' with a model as an outside function
                lower=lower,                                    # Input the lower parameter range
                upper=upper,                                    # Input the upper parameter range
                method="spso2011",                              # Specify the algorithm to calibrate the model
                control=list(write2disk=TRUE, MinMax="max",     # These are the parameters for the HydroPSO calibration approach
                             npart=55, maxit=9, normalise=TRUE, # See the documentation for each of these items # Before i had at npart = 80 and maxit =25; the typical runs were 50 npart * 5 iter
                             REPORT=10, parallel="none",        # The model can either be run a certain number of iterations or it can reach a cut off threshold
                             reltol=1E-10),                     # There are also coefficients for how widely the search space is and another coefficient to control how quickly convergence occurs (Not shown, these are currently defaults)
                model.FUN="Dynatophydromod",                    # This is the function that calls dynatop
                model.FUN.args= list(obs = Q.obs.calib,                            # Timeseries of observed Q data to compare the simulated results #inner.timestep=in.timestep,            # Number of inner timesteps
                                     inner.timestep = in.timestep,                 # Inner timestep
                                     rain = rain.calib,                            # Timeseries of precipitation at dt to run the model 
                                     PET = PET.calib,                              # Timeseries of PET at dt needed to run the model 
                                     disc = FR.disc,                               # Spatial discretization (weighting matrix) and 'groups' parameter matrix 
                                     RoutingTable = NULL,                          # Routing table to discretize stream network
                                     date = dates.cal.warm,                        # Vector of dates for the calibration period
                                     dates.df = dates.df,                          # Date frame of dates used for warm up and calibration periods
                                     model.drty=model.drty),                       # Input the model directory                              
                obs = Q.obs.calib,                              # Note above we're the model.FUN.args wasn't working properly, so we're adding this back in as variables that will feed into the model
                inner.timestep = in.timestep,                   # Inner timestep   
                rain = rain.calib,                              # Timeseries of precipitation at dt to run the model 
                PET = PET.calib,                                # Timeseries of PET at dt needed to run the model
                disc = FR.disc,                                 # Spatial discretization (weighting matrix) and 'groups' parameter matrix 
                RoutingTable = NULL,                            # Routing table to discretize stream network
                date = dates.cal.warm,                          # Array of dates for each time step
                dates.df = dates.df,                            # Date frame of dates used for warm up and calibration periods
                model.drty = model.drty,                        # Model directory    
                # Date frame of dates used for warm up and calibration periods       
)                                                                                       
                                                                                         
# Read the optimal parameter values and rerun the model                                
optim.params <- out$par                                                            # Read the optimal parameter set
names(optim.params) <- names(params)                                               # Set names 
result.optim <- runPSOCalibrationTopmodel(param.values=optim.params,               # Run the model with the optimal parameter set
                                          inner.timesteps = in.timestep,           # input the inner.timesteps
                                          rains = rain.calib, PETs=PET.calib,      # All forcing components are the same
                                          obss=Q.obs.calib, discs=FR.disc,         # Spatial information comes from the dynatop spatial function
                                          RoutingTables=NULL,                      # See above
                                          dates.dfs = dates.df)                    # Calibration end date input
sim.zoo <- (window(result.optim$sim, start = dates.df$calibration.initial,         # window the simulation
                         end = dates.df$calibration.final))                        # Will convert to zoo obj later
obs.zoo <- (window(Q.obs.calib, start = dates.df$calibration.initial,              # Window the observed
                  end = dates.df$calibration.final))                               # Will convert to zoo obj later
                                                                
if (length(obs.zoo)<length(sim.zoo)) {                          # Ensure that the lengths of the obs and sim are the same
  sim.zoo <- sim.zoo[1:length(obs.zoo)]                         # Just to makes sure they are the same length
} else {                                                        # Sometimes they can get off if there's problems with the aggregate function in dynatop
  obs.zoo <- obs.zoo[1:length(sim.zoo)]                         # Unsure why this occurs, but for example of the dt is 2 then there is like a 4 time step offset
}                                                               # TS should be the same length now

sim.zoo <- zoo(sim.zoo)                                         # Convert to zoo object
obs.zoo <- zoo(obs.zoo)                                         # Convert to zoo object

# Plot the optimal simulation and the observed q
#x11()                                                           # Open a new window
#plot(sim.zoo)                                                   # Plot the simulation
#lines(obs.zoo, col = "red")                                     # Change the color of the line
#plot.zoo(cbind(sim.zoo, obs.zoo))                               # Plot the observed
#plot.zoo(cbind(sim.zoo, obs.zoo),                               # Combine the sim and obs
#         plot.type = "single",                                  # Plot type specified
#         col = c("red", "blue"),                                # Colors of the lines
#         ylab='Q (m/hr)',                                       # Specify y axis
#         xlab='date')                                           # Specify x axis

# Process the results of the calibration for the ggof plot
#sprintf('The KGE is: %2f',out$value)                            # Print the kge score
#main.title <- paste0('Falling Rock Dynatopmodel log(3 x m)')    # Specify title for the simulation
#sim.numeric <- as.numeric(sim.zoo)                              # Write sim as a numeric 
#obs.numeric <- as.numeric(obs.zoo)                              # Write obs as a numeric           # Write obs a zoo          
#x11()                                                           # Create a new window
#ggof(sim=sim.zoo,obs=obs.zoo, main=main.title,                  # Run the ggof which also creates a plot
#     ylab="Q, [m/hr]", cex=0.3, lty=c(1,1) )                    # Specifications of the ggof plot

# Plot the results of the PSO calibration
#residuals <- sim.zoo.new - obs.zoo.new                          # Used to plot the difference in the simulated and observed Q (currently unused)
x11()                                                           # Open a new window
main.title <- paste0('Falling Rock Dynatopmodel log(m)')        # Specify title for the simulation
plot_results(do.png=FALSE,                                      # plots do not go to PNG figures instead of to the screen
             MinMax="max",                                      # 'best' parameter set is the one with the highest GoF
             beh.thr=0.3,                                       # parameter sets with GoF >= 0.3 will be considered as behavioural
             do.pairs=TRUE,                                     # a basic correlation plot is produced among parameters and GoF
             gof.name="KGE2012",                                # name used for axis titles in plots related to parameters
             legend.pos="right",                                # position of the legend in the convergence plot
             main=main.title)                                   # user-defined title for most of the output figures

# Read the results of the calibration and process the results into quantiles 
# NOTE: the code was having trouble reading in the obs file produced by hydroPSO. I just had to delete the title column from the text file and it worked properly.
# Read the results of the calibration and process the results into quantiles 
# NOTE: the code was having trouble reading in the obs file produced by hydroPSO. I just had to delete the title column from the text file and it worked properly.
PSO.drty <- paste0("/PSO.out/")                                 # Specify the directory where all of the PSO files are written to
res <- read_results(drty.out=PSO.drty, MinMax="max",            # Read the results of the calibration; MAX is noting that we are looking for maximum KGE vals
                    beh.thr=0.3,modelout.cols=NULL)             # 0.3 is the threshold for specifying behavioral model realizations
params <- res[["params"]]                                       # Record the optimal parameters from the results
gofs <- res[["gofs"]]                                           # Read the goodness of fit values
Qsims <- res[["model.values"]]                                  # Model values for the simulations
model.best <- res[["model.best"]]                               # Read the best model results
model.obs <- res[["model.obs"]]                                 # Outputthe observed results 
optim.params <- data.frame(t(res$best.param))                   # Transpose the data and write as dataframe to be read in by the model
obs.zoo <- (window(Q.obs.calib,                                 # Windowing the obs.zoo
                   start = dates.df$calibration.initial,        # start
                   end = dates.df$calibration.final))           # End
obs.zoo <- zoo(obs.zoo[1:length(model.best)])                   # Convert to zoo and set length to the length of the best model. 
dates.cal <- time(obs.zoo)                                      # Extract the dates from the obs.zoo

result.optim.run <- runPSOCalibrationTopmodel(param.values=optim.params[[1]],      # Run the model with the optimal parameter set
                                              inner.timesteps = in.timestep,           # input the inner.timesteps
                                              rains = rain.calib, PETs=PET.calib,      # All forcing components are the same
                                              obss=Q.obs.calib, discs=FR.disc,         # Spatial information comes from the dynatop spatial function
                                              RoutingTables=NULL,                      # See above
                                              dates.dfs = dates.df)                    # Calibration end date input


result.optim <- zoo(model.best,dates.cal)                       # Create zoo of the best model
obs.zoo.new <- zoo(model.obs,dates.cal)                         # Create new zoo of the model obs

# Calculate statistics of model run
pbias(sim=result.optim,obs=obs.zoo.new)                                            # Calculate pbias of sim and obs
gof.logKGE <- KGE(sim=as.numeric(log(result.optim+1e-10)),                         # Calculate log KGE
                  obs=as.numeric(log(model.obs+1e-10)), method="2012")             # Continued
gof.logNSE <- dynatopmodel::NSE(qsim=as.numeric(log(result.optim+1e-10)),          # Calculate log NSE
                                qobs=as.numeric(log(model.obs+1e-10)))                           # Continued
gof.KGE <- KGE(sim=as.numeric(result.optim+1e-10),obs=as.numeric(model.obs+1e-10), # Calculate KGE
               method='2012')                                                      # Continued 
gof.NSE <- dynatopmodel::NSE(qsim=as.numeric(result.optim+1e-10),                  # Calculate NSE 
                             qobs=as.numeric(model.obs+1e-10))                                   # Continued 
gof.metrics <- data.frame(logKGE=gof.logKGE,KGE=gof.KGE,                           # Record all of the test stats
                          logNSE=gof.logNSE,NSE=gof.NSE)                           # Continued 
gof.metrics                                                                        # print stats


## Process Quantiles of the model results 
params.025.50.975 <- wquantile(params, weights=gofs,            # Determine the weighted quantiles based on the gofs
                               byrow=FALSE,                     # Specifying wquanitle inputs
                               probs=c(0.025, 0.5, 0.975),      # probabilities to be calculated
                               normwt=TRUE, verbose=FALSE)      # Other inputs
params.025.50.best.975 <- cbind(params.025.50.975[, c(1,2)],    # Cbind to determine how far past the median
                                Best=as.numeric(optim.params[[1]]),                # This is realated to the beset parameters
                                params.025.50.975[, 3] )        # Which column is needed 
colnames(params.025.50.best.975)[4] <- "97.5%"                  # Column names
round( params.025.50.best.975, 2)                               # Rounding the parameters - for viewing
n <- ncol(Qsims)                                                # number of time steps
Qsim.025.q50.q975 <- wquantile(Qsims, weights=gofs, byrow=FALSE,# This is based on the Qsims this time 
                               probs=c(0.025, 0.5, 0.975),      # This is the quantiles that will be calculated
                               normwt=TRUE, verbose=FALSE)      # It's normalized
round( head(Qsim.025.q50.q975), 3)                              # Prints the quantiles Qsim
q025 <- zoo(Qsim.025.q50.q975[,1], dates.cal)                   # Transform the 2.5 quantiles to zoo
q975 <- zoo(Qsim.025.q50.q975[,3], dates.cal)                   # Transform the 97.5 quantile to zoo
q025 <- window(q025,start=dates.df$calibration.initial,         # Window the q025
               end=dates.df$calibration.final)                  # Continued
q975 <- window(q975,start=dates.df$calibration.initial,         # Window the q975
               end=dates.df$calibration.final)                  # Continued
pf <- pfactor(x=obs.zoo.new, lband=q025, uband=q975, na.rm=TRUE)# Calculate the P factor
pf                                                              # Print the P factor
rf <- rfactor(x=obs.zoo.new, lband=q025, uband=q975, na.rm=TRUE)# Calculate the R factor
rf                                                              # Print the R factor 

precip.out <- window(rain.calib,                                # Window precip for writing model
                     start=dates.df$calibration.initial,        # Start
                     end=dates.df$calibration.final)            # End
obs.out <- window(Q.obs.calib,                                  # Window obs for writing model
                  start=dates.df$calibration.initial,           # Start
                  end=dates.df$calibration.final)               # End
optim.out <- window(result.optim.run$run$qsim,                       # Window optim result for writing model
                    start=dates.df$calibration.initial,         # Start
                    end=dates.df$calibration.final)             # End
precip.out <- precip.out[1:12502]                               # Ensure length is 12502 - NOTE THIS IS BAD PRACTICE AND WILL NEED TO CHANGE IF dt CHANGES
obs.out <- obs.out[1:12502]                                     # Ensure length is 12502 - NOTE THIS IS BAD PRACTICE AND WILL NEED TO CHANGE IF dt CHANGES

model.save <- data.frame(precip=precip.out,obs=obs.out,         # Compile data frame for writnig out
                         optim=optim.out,q025=q025,q975=q975)   # Compile continued 
#write.csv(model.save,file='modelsave.csv')                      # Save it


# Plot the 95PPU
x11()                                                           # Open a new window
main.uncert <- paste0("Uncertainty bounds Dynatopmodel:",       # Title of the uncertainty
                      " Falling Rock")                          # Continued
Qobs.col <- "black"                                             # Color of the q obs
best.sim.col <- "blue"                                          # Color of the best simulation
bands.col <- "lightblue"                                        # Color of the uncertainty bands
Qsim.best.zoo <- zoo(result.optim, time(obs.zoo.new) )          # Zoo of the best Qsim
plot(obs.zoo.new, xaxt="n", xlab="", type="n", main=main.uncert,# Plot the best
     ylab="Q, [m/hr]",ylim=c(0,.003))                           # empty plotting area
drawTimeAxis(obs.zoo.new)                                       # The time axis
plotbandsonly(lband=q025, uband=q975)                           # Plot the 95 PPU
lines(obs.zoo.new, col=Qobs.col, lwd=0.3)                       # Plot Qobs
lines(Qsim.best.zoo, col=best.sim.col, lwd=0.3)                 # Plot best simulation
grid()                                                          # Turn on the grid
legend("topright", legend=c("Qobs", "best.sim", "95PPU"),       # Add in the legend
       lty=c(1, 1, NA),                                         # Specifcs of the legend
       pch=c(NA, NA, 15), col=c(Qobs.col, best.sim.col,         # Continued
                                bands.col),                     # Continued 
       bty="n", cex = 1.2)                                      # Continued 
legend("topleft", legend=c(paste0("P-factor: ", round(pf, 2)),  # Legend for P and R factors 
                           paste0("R-factor: ", round(rf, 2)) ),# R factor
       bty="n", cex = 1.2)                                      # Specifics of the legend

## DYNATOP POSTPROCESSING CODE

## Run explicit Routing function for the optimal run
# Note - to calculate the Q in mm/s I divided by the area upstream of the outlet of each reach...
#explicit.routing <- explicit.routing.instant(FR.Spatial,                           # Run the explicit routing function
#                                             run=result.optim.run$run)             # Returns Q in m/hr and m^3/hr for each reach
time.sim <- time(result.optim.run$sim)                                             # Get the time of the sim out
time.qbf <- time(result.optim.run$run$fluxes$qbf)                                  # Get the time of the qbf sim

## Read in the 1/0 timeseries
# Read data from the loggers stored in the following directory
logger.dir <- paste0(getwd(),'/LoggerData',sep='')
files.logger <- list.files(logger.dir,pattern='*clean.csv')                        # List files in the above directory
files.logger <- files.logger[1:4]                                                  # Right now only care about the FC sensors (1-4) in the directory
Reach.identifiers <- c('X3722','X3242','X2538','X810') # is FC3 2538?? 
Reach.identifiers.second <- c('3722','3242','2538','810')
TS.zoo.all <- list()
TS.logger.clean.all <- list()
for (file.name in 1:length(files.logger)) {                                                  # Iterate over each file
  # Read in data and format
  short.name <- substr(files.logger[file.name],start=1,
                       stop=(nchar(files.logger[file.name])-10))               # Create name for the iteration
  TS.logger <- read.csv(paste0(logger.dir,'/',files.logger[file.name]))                          # Read the CSV for the logger
  TS.logger.clean <- data.frame(Date=TS.logger$RoundTime,                          # Get only the relevant data Date
                                State=TS.logger$ApproxState)                       # Read approx state into data frame
  TS.logger.clean$State[which(TS.logger.clean$State=='Missing data')] <- 3         # Set missing data to 3
  TS.logger.clean$Date <- (ymd_hms(TS.logger.clean$Date))                          # Convert naming convention of the date
  TS.logger.clean$Date2 <- as.Date(TS.logger.clean$Date)                           # Convert to base date 
  TS.logger.clean$State <- as.numeric(TS.logger.clean$State)                       # Get the state
  TS.zoo <- zoo(TS.logger.clean$State,TS.logger.clean$Date)                        # Convert to zoo
  
  # Calculate the statistics of the TS
  on <- sum(TS.logger.clean$State==1, na.rm=T)                                     # Sum the 'wet' instances
  off <- sum(TS.logger.clean$State==0, na.rm =T)                                   # Sum the 'dry' instances
  Total.logger <- on+off                                                           # Get the total number of readings
  percent.on <- on/Total.logger*100                                                # Convert wet to percent

  if (file.name == 1) {                                              # write out the percentage on
    name = paste0(short.name,'percent')                                            # First file
    logger.out.percent <- data.frame(percent.on)                                   # write percent on
  } else {                                                                         # other files 
    name = paste0(short.name,'percent')                                            # set thier name
    logger.out.percent <- cbind(logger.out.percent,percent.on)                     # bind the reading to the previous
  }
  TS.zoo.all[[file.name]] <- TS.zoo
  TS.logger.clean.all[[file.name]] <- TS.logger.clean
}
colnames(logger.out.percent)=files.logger                                          # Set the names of the columns of the percent out 
names(TS.zoo.all) <- files.logger
names(TS.logger.clean.all) <- files.logger 

# Read in the files for behavorial simulations
# Set the directory NOTE: USER MUST CHANGE THE FOLDER FOR WHATEVER CALIBRATION TEST WE ARE ON, AS OF 6/27/2022 THIS IS CalibrationOutput from Github
qin.files.dir <- paste0(model.drty,'/fluxes_stores',sep='') #' C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/Calibration/CalibrationTest13/fluxes_stores'
qin.files <- list.files(qin.files.dir,pattern="fluxqin")        # Read in all the files with the 'fluxqin' name

headwater.postprocess <- c()                                    # Initialize the postprocessing vector 

# CHANGE THE 400 TO 1 LATER
# Run the post processing code (Headwater.evaluation) for the behavioral parameter sets
for (qin.file in 1:length(qin.files)) {                         # Iterate over every behavioral parmeter set
  print(qin.file)                                               # Pring the qin.file number 
  print(Sys.time())
  headwater.postprocess[[qin.file]] <-                          # Run the Headwater.evaluation function for each behavioral parameter set
    Headwater.evaluation.dynamic(qin.file,FR.Spatial,                   # Inputs: qin.file and FR.Spatial configuration
                         TS.logger.clean.all,TS.zoo.all,                # Also the logger data
                         logger.out.percent,
                         Reach.identifiers,
                         Reach.identifiers.second,
                         model.timestep,
                         result.optim.run,
                         time.qbf)                    # Also the percent that the logger data is 'wet'
}

# Post process the list generated for the behavioral sets to read out the correct.state.percent for both years, for 2003, for 2005, and the flow threshold
# For both study years, this is using the average best percent threshold, deriving a flow and calculating the percent correct: 
kge.in.files <- list.files(path=qin.files.dir,pattern='gof_log_kge')
# CHANGE THE 400 TO 1 LATER
kge.in.files <- kge.in.files[1:length(qin.files)]
qin.files <- qin.files[1:length(qin.files)]
kge.in <- matrix(nrow=length(kge.in.files))
for (kge.in.file in 1:length(kge.in.files)) {
  raw.kge <- read.csv(paste0(qin.files.dir,'/gof_log_kge',kge.in.file,'.csv'))
  kge.in[kge.in.file] <- raw.kge$x
}

# Ouput a vector of the optimal flow theshold for each run 
optimal.flow.thresh <- data.frame(matrix(nrow=length(qin.files),ncol=1))

for (iter.postprocess in 1:length(qin.files)) {
  optimal.flow.thresh[iter.postprocess,1] <- headwater.postprocess[[iter.postprocess]]$thresh.best       # This is the optimal transmissivity that was calibrated. m^3/hr
}
optimal.flow.thresh.mm.day <- optimal.flow.thresh
colnames(optimal.flow.thresh.mm.day) <- 'flow.thresh'

# Output a matrix of the percent correct for each optimal threshold for each FC sensor will be size  4 x number of behavioural parameter sets - #Run 105 is the optimal run for FC4
percent.correct.state.FC <- data.frame(matrix(nrow=length(qin.files),ncol=length(Reach.identifiers)))

for (iter.postprocess in 1:length(qin.files)) { 
  for (Reach.id in 1:length(Reach.identifiers)) {
    percent.correct.state.FC[iter.postprocess,Reach.id] <- headwater.postprocess[[iter.postprocess]]$flow.thresh.best.results$out.thresh.data[[Reach.id]]$percent.correct[which(headwater.postprocess[[iter.postprocess]]$flow.thresh.best.results$out.thresh.data[[Reach.id]]$transmissivity==optimal.flow.thresh[iter.postprocess,1])]
     }
  }

colnames(percent.correct.state.FC) <- Reach.identifiers

# Output a vector of the total error for each optimal threshold
total.error.state.FC <- data.frame(matrix(nrow=length(qin.files),ncol=length(Reach.identifiers)))

for (iter.postprocess in 1:length(qin.files)) { 
  for (Reach.id in 1:length(Reach.identifiers)) {
    total.error.state.FC[iter.postprocess,Reach.id] <- headwater.postprocess[[iter.postprocess]]$flow.thresh.best.results$out.thresh.data[[Reach.id]]$sum.incorrect[which(headwater.postprocess[[iter.postprocess]]$flow.thresh.best.results$out.thresh.data[[Reach.id]]$transmissivity==optimal.flow.thresh[iter.postprocess,1])]
  }
}

# output a vector of the total correct for each optimal threshold
total.correct.state.FC <- data.frame(matrix(nrow=length(qin.files),ncol=length(Reach.identifiers)))

for (iter.postprocess in 1:length(qin.files)) {
  for (Reach.id in 1:length(Reach.identifiers)) {
    total.correct.state.FC[iter.postprocess,Reach.id] <- headwater.postprocess[[iter.postprocess]]$flow.thresh.best.results$out.thresh.data[[Reach.id]]$sum.correct[which(headwater.postprocess[[iter.postprocess]]$flow.thresh.best.results$out.thresh.data[[Reach.id]]$transmissivity==optimal.flow.thresh[iter.postprocess,1])]
  }
}

Q.sub.c.reaches.all <- (matrix(ncol=length(qin.files),nrow=length(FR.Spatial$explicit.ChanTable[,1])))
for (iter.postprocess in 1:length(qin.files)) {
  Q.sub.c.reaches.all[,iter.postprocess] <- headwater.postprocess[[iter.postprocess]]$Q.sub.c.best
}
Q.sub.c.mean <- rowMeans(Q.sub.c.reaches.all)
Q.sub.c.stdev <- rowSds(Q.sub.c.reaches.all)

Optim.transmissivity

# Calculate the total percent correct for the entire study period over all of the FC sensors
total.percent.correct.FC <- data.frame(sum.correct=rowSums(total.correct.state.FC),sum.incorrect=rowSums(total.error.state.FC))
total.percent.correct.FC$percent.correct <- total.percent.correct.FC$sum.correct/(total.percent.correct.FC$sum.correct+total.percent.correct.FC$sum.incorrect)*100

# Output for each behavioural parameter set the percentage that each reach is on 
network.percent.on <- (matrix(nrow=length(FR.Spatial$explicit.ChanTable$link_no),ncol=length(qin.files)))

for (iter.postprocess in 1:length(qin.files)) {
  network.percent.on[,iter.postprocess] <- headwater.postprocess[[iter.postprocess]]$percent.on.network
}

network.percent.on.df <- data.frame(network.percent.on)
rownames(network.percent.on.df) <- paste0('X',FR.Spatial$explicit.ChanTable$link_no)

# Calculate the mean and standard deviation of the percent on 
network.percent.on.stats <- data.frame(mean=rowMeans(network.percent.on))
network.percent.on.stats$stdev <-rowSds((network.percent.on))
network.percent.on.stats$var <- rowVars(network.percent.on)

# Output for each time series with which the flow is on and off as a list 
simulated.network.flow <- list()

for (iter.postprocess in 1:length(qin.files)) {
  simulated.network.flow[[iter.postprocess]] <- headwater.postprocess[[iter.postprocess]]$flow.network.thresh # tHIS IS IN M^3/HR
}

# determine the total network length for each time step by: multiplying the on/off state of the reach by the reach length, then doing a row sum.
# Record the network length at each timestep
reach.length.m <- FR.Spatial$explicit.ChanTable$stream_length_ft*0.3048
network.length.postprocess <- matrix(nrow=nrow(result.optim.run$run$fluxes$qbf),ncol=length(qin.files))
for (iter.postprocess in 1:length(qin.files)) { 
  network.length.postprocess[,iter.postprocess] <- rowSums(sweep(headwater.postprocess[[iter.postprocess]]$total.network.on.off,MARGIN=2,reach.length.m,`*`))
}
network.length.total.fig <- (headwater.postprocess[[qin.file]]$total.network.on.off[1:nrow(network.length.postprocess),])
row.names(network.length.total.fig) <- q_specific.run$time
#rownames(network.length.total) <- seq(1,12504)

# Convert from m/hr (the output of dynatopmmodel) to mm/day: 1 m / 1 hr * 1000 mm / 1 m * 24 hr / 1 day
sim.plot <- Qsim.best.zoo*24000
obs.plot <- obs.zoo.new*24000
lower.plot <- q025*24000
upper.plot <- q975*24000
FC4.plot <- result.test.FC4$sim[2210:length(result.test.FC4$sim)]*24000
upper.plot <- ifelse(FC4.plot<upper.plot,upper.plot,FC4.plot)+.2
obs.minus.sim.plot <- obs.plot-sim.plot
outlet.sim.obs <- data.frame(time=as.Date(time(Qsim.best.zoo)),sim=sim.plot,obs=obs.plot,
                             lower=lower.plot,upper=upper.plot,FC4.plot=FC4.plot,obs.minus.sim=obs.minus.sim.plot) # This code converts everything from m/hr to mm/day
#write.csv(outlet.sim.obs,'C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/6 WRITE UP/SciHub/dischargeTimeSeries.csv')
#upper.lower.band <- data.frame(lower=q025,upper=q975)
#require(scales)
sim.obs.outlet.plot <- ggplot(data=outlet.sim.obs,aes(x=time,y=sim))+geom_ribbon(data=outlet.sim.obs,aes(ymin=lower,ymax=upper),fill='grey70')+geom_line(data=outlet.sim.obs,aes(x=time,y=obs),color='red',linetype='solid') +
  geom_line(data=outlet.sim.obs,aes(x=time,y=FC4.plot),color='aquamarine4',linetype='solid')+theme(legend.position = "none")
sim.obs.outlet.plot <- sim.obs.outlet.plot+geom_line(color='black',size=.55)+ theme_bw() + xlab('Date') +
  scale_x_date(date_breaks = '6 months')+ylab('Discharge (mm/d)*')+ylim(0,20.0)
sim.obs.outlet.plot <-  sim.obs.outlet.plot + theme(axis.title.x=element_blank())
sim.obs.outlet.plot

sim.obs.outlet.plot.2 <- ggplot(data=outlet.sim.obs,aes(x=time,y=sim))+geom_ribbon(data=outlet.sim.obs,aes(ymin=lower,ymax=upper),fill='grey70')+geom_line(color='black',size=.55)+theme(legend.position = "none")
sim.obs.outlet.plot.2 <- sim.obs.outlet.plot.2+geom_ma(ma_fun=SMA,n=48,data=outlet.sim.obs,aes(x=time,y=obs),color='red',linetype='solid') + theme_bw() + xlab('Date') +
  scale_x_date(date_breaks = '6 months')+ylab('Discharge (mm/d)*')+ scale_y_continuous(trans='log10',labels=label_number(accuracy =0.01),limits=c(0.01,100))#+ylim(0,20.0)
sim.obs.outlet.plot.2 <-  sim.obs.outlet.plot.2 + theme(axis.title.x=element_blank())
sim.obs.outlet.plot.2

mean.network.length <- rowMeans(network.length.postprocess)
mean.network.length.df <- data.frame(length=mean.network.length,time=time.qbf)
median.network.length <- rowMedians(network.length.postprocess)
sd.network.length <- rowSds(network.length.postprocess)
sd.network.length.min <- mean.network.length-sd.network.length
sd.network.length.max <- mean.network.length+sd.network.length
plot.network.length <- data.frame(length=mean.network.length,time=time.qbf,sd.min=sd.network.length.min,sd.max=sd.network.length.max)
network.length.start <- which(plot.network.length$time==time(Qsim.best.zoo[1]))
plot.network.length <- plot.network.length[network.length.start:length(plot.network.length$length),]

length.plot <- ggplot(data=plot.network.length,aes(x=time,y=length))+xlab('Date')+ylab('Mean network length (m)')
#length.plot <- length.plot + geom_ma(ma_fun=SMA,n=36, linetype='solid') 
length.plot <- length.plot + geom_ribbon(data=plot.network.length, aes(ymin=sd.min,ymax=sd.max),fill='grey')+theme_bw()+
  geom_ma(ma_fun=SMA,n=12, linetype='solid') 

plot.network.percent <- data.frame(percent=100*plot.network.length$length/sum(FR.Spatial$explicit.ChanTable$stream_length_ft*0.3048),
                                   time=as.Date(plot.network.length$time),
                                   sd.min.percent=100*plot.network.length$sd.min/sum(FR.Spatial$explicit.ChanTable$stream_length_ft*0.3048),
                                   sd.max.percent=100*plot.network.length$sd.max/sum(FR.Spatial$explicit.ChanTable$stream_length_ft*0.3048))
percent.length.plot=ggplot(data=plot.network.percent,aes(x=time,y=percent))+xlab('Date') + ylab('Mean Percent Flowing (%)') +
  geom_ribbon(data=plot.network.percent,aes(ymin=sd.min.percent,ymax=sd.max.percent),fill='grey') +
  geom_ma(ma_fun=SMA,n=12,linetype='solid',size=1.05,color='black') + ylim(40,99) + theme_bw() +
  scale_x_date(date_breaks = '6 months') + theme(axis.title.x=element_blank())

ggplot(data=plot.network.percent, aes(x=percent)) + geom_histogram(binwidth=.9,colour = 'black', fill = 'white') + theme_classic()

percent.length.discharge <- data.frame(date=plot.network.percent$time[1:12502],percent=plot.network.percent$percent[1:12502],discharge=outlet.sim.obs$sim)

#roll.mean.network <- rollmean(zoo(plot.network.length$length),k=36,align='center')
#length.plot.2 <- length.plot+geom_line(data=roll.mean.network,aes(y=length))
gsim.obs.outlet.plot <- ggplotGrob(sim.obs.outlet.plot)
gpercent.length.plot <- ggplotGrob(percent.length.plot)
maxWidth = grid::unit.pmax(gsim.obs.outlet.plot$widths[2:5],gpercent.length.plot$widths[2:5])
gsim.obs.outlet.plot$widths[2:5] <- as.list(maxWidth)
gpercent.length.plot$widths[2:5] <- as.list(maxWidth)
options(scipen=10000)
x11()
grid.arrange(gsim.obs.outlet.plot,gpercent.length.plot,ncol=1)



sim.obs.outlet.plot.1 <- ggplot(data=outlet.sim.obs,aes(x=time,y=obs.minus.sim.plot))+#geom_line(data=outlet.sim.obs,aes(x=time,y=obs),color='red',linetype='solid') +
  theme(legend.position = "none")
sim.obs.outlet.plot.1 <- sim.obs.outlet.plot.1+geom_line(color='black',size=.55)+ theme_bw() + xlab('Date') +
  scale_x_date(date_breaks = '6 months')+ylab('Discharge Deficit (mm/d)*')+ylim(-20,20.0)
sim.obs.outlet.plot.1 <-  sim.obs.outlet.plot.1 + theme(axis.title.x=element_blank())
sim.obs.outlet.plot.1

sim.obs.outlet.plot.2 <- ggplot(data=outlet.sim.obs,aes(x=time,y=sim))+geom_ribbon(data=outlet.sim.obs,aes(ymin=lower,ymax=upper),fill='grey70')+geom_line(data=outlet.sim.obs,aes(x=time,y=obs),color='red',linetype='solid') +
  geom_line(data=outlet.sim.obs,aes(x=time,y=FC4.plot),color='aquamarine4',linetype='solid')+theme(legend.position = "none")
sim.obs.outlet.plot.2 <- sim.obs.outlet.plot.2+geom_line(color='black',size=.55)+ theme_bw() + xlab('Date') +
  scale_x_date(date_breaks = '6 months')+ylab('Discharge (mm/d)*')+ylim(0,20.0)
sim.obs.outlet.plot.2 <-  sim.obs.outlet.plot.2 + theme(axis.title.x=element_blank())
sim.obs.outlet.plot.2


gsim.obs.outlet.plot.1 <- ggplotGrob(sim.obs.outlet.plot.1)
gsim.obs.outlet.plot.2 <- ggplotGrob(sim.obs.outlet.plot.2)
gpercent.length.plot <- ggplotGrob(percent.length.plot)
maxWidth = grid::unit.pmax(gsim.obs.outlet.plot.1$widths[2:5],gsim.obs.outlet.plot.2$widths[2:5],gpercent.length.plot$widths[2:5])
gsim.obs.outlet.plot.1$widths[2:5] <- as.list(maxWidth)
gsim.obs.outlet.plot.2$widths[2:5] <- as.list(maxWidth)
gpercent.length.plot$widths[2:5] <- as.list(maxWidth)
options(scipen=10000)
x11()
grid.arrange(gsim.obs.outlet.plot.1,gsim.obs.outlet.plot.2,gpercent.length.plot,ncol=1)

q.vs.length <- data.frame(Q.sim=sim.plot,length = plot.network.length$length[1:length(sim.plot)],time=plot.network.length$time[1:length(sim.plot)])
View(q.vs.length)
q.vs.length.plot <- ggplot(data=q.vs.length,aes(x=Q.sim,y=length/1000))+geom_point()+scale_x_log10()+scale_y_log10()+ylab('Flowing Stream Length (km)')+
  xlab('Discharge (mm/day)')+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

length(q.vs.length$Q.sim[q.vs.length$Q.sim>2])
q.vs.length.2 <- q.vs.length[q.vs.length$Q.sim<1,]
q.vs.length.3 <- q.vs.length[q.vs.length$Q.sim>=1,]
q.vs.length$class <-ifelse(q.vs.length$Q.sim>1,'High Flow','Low Flow')

h.l.plot <- ggplot() + geom_point(data=q.vs.length.2,aes(x=Q.sim,y=length/1000),color='black') +
  geom_point(data=q.vs.length.3,aes(x=Q.sim,y=length/1000,color='grey')) + scale_x_log10()+
  scale_y_log10()+geom_smooth()
h.l.plot

h.l.plot.low <- ggplot(data=q.vs.length,aes(x=Q.sim,y=length/1000)) +geom_point()+ geom_smooth(aes(color=class),method=lm,fullrange=TRUE)+scale_x_log10()+
  scale_y_log10(limits=c(2,6)) + scale_color_brewer(palette="Paired")+theme_bw() + ylab('Length (km')

h.l.plot.low <- h.l.plot.low + labs(x='Discharge (mm/day)',y='Stream Length (km)')+theme(legend.position ='none')
x11()
h.l.plot.low
x11()
q.vs.length.plot


#write.csv(q.vs.length,'qvslength.csv')
#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('green','red'))

# Create a data frame for the flow threshold, total percent correct, and the kge...
kge.totalpercentcorrect.flowthresh <- data.frame(kge=kge.in,total.percent.correct=total.percent.correct.FC$percent.correct,flow.thresh=optimal.flow.thresh.mm.day)

# Assign colors to the data frame
# based on the y values
kge.totalpercentcorrect.flowthresh$color <- rbPal(10)[as.numeric(cut(kge.totalpercentcorrect.flowthresh$total.percent.correct,breaks = 10))]
kge.totalpercentcorrect.flowthresh$color.thresh <- rbPal(10)[as.numeric(cut(kge.totalpercentcorrect.flowthresh$flow.thresh,breaks = 10))]

# Plot the percent correct versus the kge 

plot.1 <- ggplot(data=kge.totalpercentcorrect.flowthresh,
                 aes(x=kge,y=total.percent.correct,color=(flow.thresh))) + geom_point() + xlab('KGE') + ylab('Total Percent Correct (%)') +
  scale_colour_viridis_c() + labs(color='flow.thresh')
plot.2 <- ggplot(data=kge.totalpercentcorrect.flowthresh,aes(x=kge,y=((flow.thresh)),color=(total.percent.correct))) + geom_point() + xlab('KGE') + ylab('Log (flow thresh [m^3/hr])') +
  labs(color = 'Total Percent Correct (%)') + scale_color_viridis_c()
x11()
grid.arrange(plot.1,plot.2)

plot(x=kge.totalpercentcorrect.flowthresh$kge,y=kge.totalpercentcorrect.flowthresh$total.percent.correct, col=kge.totalpercentcorrect.flowthresh$color.thresh, pch=16,
     xlab='kge',ylab='total percent (%)')
plot(y=kge.totalpercentcorrect.flowthresh$flow.thresh,x=kge.totalpercentcorrect.flowthresh$kge,col=kge.totalpercentcorrect.flowthresh$color, pch=16,
     xlab='flow mm/day',ylab='kge')


# Create a data frame for the flow threshold, FC1 percent correct, and the kge and plot
kge.FC1percentcorrect.flowthresh <- data.frame(kge=kge.in,FC1.percent.correct=percent.correct.state.FC$X3722,flow.thresh=(optimal.flow.thresh.mm.day))
kge.FC1percentcorrect.flowthresh$color <- rbPal(10)[as.numeric(cut(kge.FC1percentcorrect.flowthresh$FC1.percent.correct,breaks = 10))]
kge.FC1percentcorrect.flowthresh$color.thresh <- rbPal(10)[as.numeric(cut(kge.FC1percentcorrect.flowthresh$flow.thresh,breaks = 10))]
plot(x=kge.FC1percentcorrect.flowthresh$kge,y=kge.FC1percentcorrect.flowthresh$FC1.percent.correct, col=kge.totalpercentcorrect.flowthresh$color)
plot(x=kge.FC1percentcorrect.flowthresh$kge,y=kge.FC1percentcorrect.flowthresh$flow.thresh, col=kge.FC1percentcorrect.flowthresh$color)

# Create a data frame for the flow threshold, FC2 percent correct, and the kge and plot
kge.FC2percentcorrect.flowthresh <- data.frame(kge=kge.in,FC2.percent.correct=percent.correct.state.FC$X3242,flow.thresh=(optimal.flow.thresh.mm.day))
kge.FC2percentcorrect.flowthresh$color <- rbPal(10)[as.numeric(cut(kge.FC2percentcorrect.flowthresh$FC2.percent.correct,breaks = 10))]
kge.FC2percentcorrect.flowthresh$color.thresh <- rbPal(10)[as.numeric(cut(kge.FC2percentcorrect.flowthresh$flow.thresh,breaks = 10))]
plot(x=kge.FC2percentcorrect.flowthresh$kge,y=kge.FC2percentcorrect.flowthresh$FC2.percent.correct, col=kge.totalpercentcorrect.flowthresh$color, pch=16)
plot(x=kge.FC2percentcorrect.flowthresh$kge,y=kge.FC2percentcorrect.flowthresh$flow.thresh, col=kge.FC2percentcorrect.flowthresh$color, pch=16)

# Create a data frame for the flow threshold, FC3 percent correct, and the kge and plot
kge.FC3percentcorrect.flowthresh <- data.frame(kge=kge.in,FC3.percent.correct=percent.correct.state.FC$X2538,flow.thresh=(optimal.flow.thresh.mm.day))
kge.FC3percentcorrect.flowthresh$color <- rbPal(10)[as.numeric(cut(kge.FC3percentcorrect.flowthresh$FC3.percent.correct,breaks = 10))]
kge.FC3percentcorrect.flowthresh$color.thresh <- rbPal(10)[as.numeric(cut(kge.FC3percentcorrect.flowthresh$flow.thresh,breaks = 10))]
plot(x=kge.FC3percentcorrect.flowthresh$kge,y=kge.FC3percentcorrect.flowthresh$FC3.percent.correct, col=kge.totalpercentcorrect.flowthresh$color, pch=16)
plot(x=kge.FC3percentcorrect.flowthresh$kge,y=kge.FC3percentcorrect.flowthresh$flow.thresh, col=kge.FC3percentcorrect.flowthresh$color, pch=16)

# Create a data frame for the flow threshold, FC4 percent correct, and the kge and plot
kge.FC4percentcorrect.flowthresh <- data.frame(kge=kge.in,FC4.percent.correct=percent.correct.state.FC$X810,flow.thresh=(optimal.flow.thresh.mm.day))
kge.FC4percentcorrect.flowthresh$color <- rbPal(10)[as.numeric(cut(kge.FC4percentcorrect.flowthresh$FC4.percent.correct,breaks = 10))]
kge.FC4percentcorrect.flowthresh$color.thresh <- rbPal(10)[as.numeric(cut(kge.FC4percentcorrect.flowthresh$flow.thresh,breaks = 10))]
plot(x=kge.FC4percentcorrect.flowthresh$kge,y=kge.FC4percentcorrect.flowthresh$FC4.percent.correct, col=kge.totalpercentcorrect.flowthresh$color, pch=16)
plot(x=kge.FC4percentcorrect.flowthresh$kge,y=kge.FC4percentcorrect.flowthresh$flow.thresh, col=kge.FC4percentcorrect.flowthresh$color, pch=16)

# plot all four on a two by two plot
x11()
par(mfrow=c(2,2))
plot(x=kge.FC1percentcorrect.flowthresh$kge,y=kge.FC1percentcorrect.flowthresh$FC1.percent.correct, col=kge.totalpercentcorrect.flowthresh$color, pch=16, main= 'FC1',
     xlab='KGE Outlet Streamflow', ylab = 'Percent Correct (%)')
plot(x=kge.FC2percentcorrect.flowthresh$kge,y=kge.FC2percentcorrect.flowthresh$FC2.percent.correct, col=kge.totalpercentcorrect.flowthresh$color, pch=16, main= 'FC2',
     xlab='KGE Outlet Streamflow', ylab = 'Percent Correct (%)')
plot(x=kge.FC3percentcorrect.flowthresh$kge,y=kge.FC3percentcorrect.flowthresh$FC3.percent.correct, col=kge.totalpercentcorrect.flowthresh$color, pch=16, main= 'FC3',
     xlab='KGE Outlet Streamflow', ylab = 'Percent Correct (%)')
plot(x=kge.FC4percentcorrect.flowthresh$kge,y=kge.FC4percentcorrect.flowthresh$FC4.percent.correct, col=kge.totalpercentcorrect.flowthresh$color, pch=16, main= 'FC4',
     xlab='KGE Outlet Streamflow', ylab = 'Percent Correct (%)')

plot.fc1 <- ggplot(data=kge.FC1percentcorrect.flowthresh,aes(x=kge,y=FC1.percent.correct,color=flow.thresh))+ geom_point() + xlab('KGE') + ylab('Percent Flow State Correct (%)') +scale_colour_viridis_c() + labs(color='Transmissivity (m^3/hr)') + ggtitle('(a) FC1') + ylim(60,100)+theme_bw()+theme(legend.position='none')
plot.fc2 <- ggplot(data=kge.FC2percentcorrect.flowthresh,aes(x=kge,y=FC2.percent.correct,color=flow.thresh))+ geom_point() + xlab('KGE') + ylab('Percent Flow State Correct (%)') +scale_colour_viridis_c() + labs(color='Transmissivity (m^3/hr)') + ggtitle('(b) FC2') + ylim(60,100)+theme_bw()+theme(legend.position='none')
plot.fc3 <- ggplot(data=kge.FC3percentcorrect.flowthresh,aes(x=kge,y=FC3.percent.correct,color=flow.thresh))+ geom_point() + xlab('KGE') + ylab('Percent Flow State Correct (%)') +scale_colour_viridis_c() + labs(color='Transmissivity (m^3/hr)') + ggtitle('(c) FC3') + ylim(60,100)+theme_bw()+theme(legend.position='none')
plot.fc4 <- ggplot(data=kge.FC4percentcorrect.flowthresh,aes(x=kge,y=FC4.percent.correct,color=flow.thresh))+ geom_point() + xlab('KGE') + ylab('Percent Flow State Correct (%)') +scale_colour_viridis_c() + labs(color='Transmissivity (m^3/hr)') + ggtitle('(d) FC4') + ylim(60,100)+theme_bw()+theme(legend.position='none')

x11()
grid.arrange(plot.fc1,plot.fc2,plot.fc3,plot.fc4)
plot.legend <- ggplot(data=kge.FC1percentcorrect.flowthresh,aes(x=kge,y=FC1.percent.correct,color=flow.thresh))+ geom_point() + xlab('KGE') + ylab('Total Percent Correct (%)') +scale_colour_viridis_c() + labs(color=bquote('Transmissivity ('~m^3~'/hr)')) + ggtitle('FC1') + ylim(95,100)
plot.legend

hist(log(kge.FC1percentcorrect.flowthresh$flow.thresh))
Norm.Trans <- (transmissivit.norm= kge.FC1percentcorrect.flowthresh$flow.thresh-mean(kge.FC1percentcorrect.flowthresh$flow.thresh)/sd(kge.FC1percentcorrect.flowthresh$flow.thresh))
Norm.trans.in <- data.frame(transmissivity=kge.FC1percentcorrect.flowthresh$flow.thresh,norm.transmissivity=Norm.Trans)

hist.transmissivity <- ggplot(kge.FC1percentcorrect.flowthresh, aes(x=flow.thresh, fill = flow.thresh))+geom_histogram(aes(y=..density..), binwidth=0.4,colour='black',fill='white')+ # c('#440154FF','#404788FF','#39568CFF','#287D8EFF','#1F968BFF','#29AF7FFF','#55C667FF','#73D055FF','#95D840FF','#B8DE29FF','#FDE725FF')
  geom_density(size=1.25)+xlab('Transmissivity') + ylab('Density') + theme_bw()+theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) 
x11()
hist.transmissivity + scale_fill_viridis_c(lim=c(min(kge.FC1percentcorrect.flowthresh$flow.thresh),4.67),begin= 0, end =1) + scale_size(range=c(1,100)) +theme(legend.key.height = unit(1,'cm'), legend.key.width= unit(4,'cm'), legend.position='bottom')

hist.flow.thresh <- hist(kge.FC1percentcorrect.flowthresh$flow.thresh, breaks=13)

#ggplot(kge.FC1percentcorrect.flowthresh, aes(flow.thresh)) + geom_histogram(binwidth=0.4,fill=c('#440154FF','#404788FF','#39568CFF','#287D8EFF','#1F968BFF','#29AF7FFF','#55C667FF','#73D055FF','#95D840FF','#B8DE29FF','#FDE725FF')) + scale_fill_viridis_c()

# plot the FDC

#setwd('C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/Calibration/CalibrationTest13/fluxes_stores')
qin.run <- read.csv(paste0(model.drty,'fluxes_stores', 'fluxqin',qin.file,'.csv'))                           # Read qin 
qin.run <- xts(qin.run[,-1],time.qbf)                                            # Convert to xts
qbf.run <- read.csv(paste0(model.drty,'fluxes_stores', 'fluxqbf',qin.file,'.csv'))                           # Read qbf
qbf.run <- xts(qbf.run[,-1],time.qbf)                                            # Convert to xts
ae.run <- read.csv(paste0(model.drty,'fluxes_stores','fluxae',qin.file,'.csv'))                             # Read ae
ae.run <- xts(ae.run[,-1],time.qbf)                                              # Convert to xts
rain.run <- read.csv(paste0(model.drty,'fluxes_stores','fluxrain',qin.file,'.csv'))                         # Read rain
rain.run <- xts(rain.run[,-1],time.qbf)                                          # Convert to xts
qof.run <- read.csv(paste0(model.drty,'fluxes_stores','fluxqof',qin.file,'.csv'))                           # Read qof
qof.run <- xts(qof.run[,-1],time.qbf)                                            # Convert to xts
Qsim.run <- read.csv(paste0(model.drty,'fluxes_stores','qsim',qin.file,'.csv'))                             # Read Qsim
Qsim.run <- xts(Qsim.run[,-1],time.qbf)                                          # Convert to xts
q_specific.run <- read.csv(paste0(model.drty,'fluxes_stores','qsim_specific',qin.file,'.csv'))
q_specific.run <- data.frame(Q=q_specific.run[,-1]*24000,time=time.qbf)
compile.fluxes <- list(qin = qin.run, qbf = qbf.run, ae = ae.run,                # Compile fluxes 
                       rain = rain.run, qof = qof.run)                           # Into list
compile.run <- list(fluxes=compile.fluxes, Qsim = Qsim.run)                      # Compile fluxes and Qsim into list

# Run routing function -- note TOPMODEL puts out data in m/hr for reach-specific flow in runoff
explicit.routing.run <- explicit.routing.instant(FR.Spatial,                     # Read in FR.Spatial and the compiled run
                                                 run=compile.run,
                                                 model.timestep = model.timestep)# Run the explicit routing code for the simulation
Q <- data.frame(explicit.routing.run$routed.flow.instant.mm_s)                   # Assign the routed flow instant (mm/s) to a variables
r.names <- names(Q)                                                              # Get the reach names

# Convert Q from mm/s to mm/day
Q <- Q*24*60*60                                                                # This is important for some of the empirical calculations we'll later need to do

Q.compare.test <- data.frame(sim_routing=Q$X3802, sim_dtop=q_specific.run,diff = (Q$X3802-q_specific.run)/q_specific.run*100)


View(explicit.routing.run$routed.flow.instant.m3_hr)
volume.reaches <- explicit.routing.run$routed.flow.instant.m3_hr*model.timestep
total.outflow.volume <- sum(volume.reaches$`3802`)
total.firstorder.volume <- sum(volume.reaches[,17:33])
total.firstorder.volume/total.outflow.volume

temporal.volume <- matrix(nrow = length(volume.reaches$`3802`))

for (ts in 1:length(volume.reaches$`3802`)) {
  temporal.volume <- sum(volume.reaches[ts,17:33])/volume.reaches$`3802`[ts]
}

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

Ylabel_parse=parse(text='Discharge (mm d^-1*)')
flow.duration.outlet <- data.frame(flow.duration = Flow.duration.reaches$X3802,Q.in=Q.in.rank$X3802)
options(sciepen=10)
x11()
par(mar=c(5,5,5,5))
plot(x=Flow.duration.reaches$X3802,y=Q.in.rank$X3802,ylab=expression(paste('Discharge ', '(m d'^'-1',')')),yaxt='n',xlab = 'Exceedance Percentile', type ='l', log ='y')
axis(2,at=c(0.005,0.01,0.1,1,25,200,1000), labels=c(0.005,0.01,0.1,1,25,200,1000))
View(flow.duration.outlet)

#write.csv(flow.duration.outlet,'C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/6 WRITE UP/SciHub/outFlowDuration.csv')

# Read in stream network
FR_drn <- shapefile('SpatialInputData/FR4000StreamNet.shp')
FR_boundary <- shapefile('SpatialInputData/FRSubcatchments.shp')
FR_DEM <- raster('SpatialInputData/fr1meterDEM.tif')

# Add the average percent connected for the reaches 
FR_drn$percent.connected <- network.percent.on.stats$mean
test.df.percent.on <- data.frame(reach.no=FR_drn$LINKNO,percent.on=network.percent.on.stats$mean)
#setwd(paste0(getwd(),'/Shapefiles'))
#shapefile(FR_drn,'FR_drn_mean.shp')


#x11()
#plot(FR_DEM)
#lines(FR_boundary)
#lines(FR_drn)

# Plot the 95PPU
x11()                                                           # Open a new window
main.uncert <- paste0("Uncertainty bounds Dynatopmodel:",       # Title of the uncertainty
                      " Falling Rock")                          # Continued
Qobs.col <- "black"                                             # Color of the q obs
best.sim.col <- "blue"                                          # Color of the best simulation
bands.col <- "lightblue"                                        # Color of the uncertainty bands
Qsim.best.zoo <- zoo(result.optim[3000:6000], time(obs.zoo.new[3000:6000]) )          # Zoo of the best Qsim
Qsim.best.zoo.plot <- Qsim.best.zoo
plot(Qsim.best.zoo.plot, xaxt="n", xlab="", type="n", main=main.uncert,# Plot the best
     ylab="Q, [m/hr]",ylim=c(0,.0015), 
     )                           # empty plotting area
drawTimeAxis(obs.zoo.new[3000:6000])                                       # The time axis
plotbandsonly(lband=q025[3000:6000], uband=q975[3000:6000])                           # Plot the 95 PPU
                     # Plot Qobs
lines(Qsim.best.zoo.plot, col=best.sim.col, lwd=0.3)                 # Plot best simulation
grid()                                                          # Turn on the grid
legend("topright", legend=c( "best.sim", "95PPU"),       # Add in the legend
       lty=c(1, NA),                                         # Specifcs of the legend
       pch=c(NA, 15), col=c(best.sim.col,         # Continued
                                bands.col),                     # Continued 
       bty="n", cex = 1.2)                                      # Continued 
legend("topleft", legend=c(paste0("P-factor: ", round(pf, 2)),  # Legend for P and R factors 
                           paste0("R-factor: ", round(rf, 2)) ),# R factor
       bty="n", cex = 1.2)                                      # Specifics of the legend


# Calculate statistics for the post processing data
postprocess.stats <- data.frame(mean=mean(correct.state.percent$correct.percent),  # Calculate the mean correct percent for all behavioral sets
                                max=max(correct.state.percent$correct.percent),    # Calculate max correct percent 
                                min=min(correct.state.percent$correct.percent),    # Calculate min correct percent
                                sd=sd(correct.state.percent$correct.percent))      # Calculate the standard deviation
postprocess.stats.2003 <-                                                          # For 2003
  data.frame(mean=mean(correct.state.percent.2003$correct.percent),                # Calculate the mean correct percent for all behavioral sets 2003
             max=max(correct.state.percent.2003$correct.percent),                  # Calculate max correct percent 2003 
             min=min(correct.state.percent.2003$correct.percent),                  # Calculate min correct percent 2003
             sd=sd(correct.state.percent.2003$correct.percent))                    # Calculate the standard deviation 2003
postprocess.stats.2005 <-                                                          # For 2005
  data.frame(mean=mean(correct.state.percent.2005$correct.percent),                # Calculate the mean correct percent for all behavioral sets 2005
             max=max(correct.state.percent.2005$correct.percent),                  # Calculate max correct percent 2005
             min=min(correct.state.percent.2005$correct.percent),                  # Calculate min correct percent 2005
             sd=sd(correct.state.percent.2005$correct.percent))                    # Calculate the standard devaition for 2005
flow.thresh.stats <-                                                               # For the flow threshold
  data.frame(mean=mean(flow.thresh.postproccessing$flow.threshold),                # Calculate the mean flow threhsold for all behavioral sets 
             max=max(flow.thresh.postproccessing$flow.threshold),                  # Calculate max flow threshold
             min=min(flow.thresh.postproccessing$flow.threshold),                  # Calculate min flow threshold
             sd=sd(flow.thresh.postproccessing$flow.threshold))                    # Calculate the standard deviation for flow threshold
percent.thresh.stats <-                                                            # For the flow threshold
  data.frame(mean=mean(percent.thresh.postproccessing$percent.threshold),          # Calculate the mean flow threhsold for all behavioral sets 
             max=max(percent.thresh.postproccessing$percent.threshold),            # Calculate max flow threshold
             min=min(percent.thresh.postproccessing$percent.threshold),            # Calculate min flow threshold
             sd=sd(percent.thresh.postproccessing$percent.threshold))              # Calculate the standard deviation for flow threshold

# Plot histograms of the results 
plot(hist(correct.state.percent$correct.percent),                                  # Plot the histogram for the correct percent for the study years
     main='Histogram correct percent',xlab='correct percent')                      # Title and labels
plot(hist(correct.state.percent.2003$correct.percent),                             # Plot the histogram for the correct percent 2003 
     main='Histogram correct percent 2003',xlab='correct percent')                 # Title and labels
plot(hist(correct.state.percent.2005$correct.percent),                             # Plot the histogram for the correct percent 2005 
     main='Histogram correct percent 2005',xlab='correct percent')                 # Title and labels 
plot(hist(log(flow.thresh.postproccessing$flow.threshold)),                        # Plot the histogram for the log of the flow threshold - log because flow threshold is so small
     main='Flow threshold',xlab='log(Runoff [mm/s])')                              # Title and labels
plot(hist((flow.thresh.postproccessing$flow.threshold)),                           # Plot the histogram for the log of the flow threshold - log because flow threshold is so small
     main='Flow threshold',xlab='(Runoff [mm/s])')                                 # Title and labels
plot(hist(percent.thresh.postproccessing$percent.threshold),                        # Plot the histogram for the log of the flow threshold - log because flow threshold is so small
     main='Percent threshold',xlab='[%]')                                          # Title and labels


#instant.kinematic <- explicit.routing.run$routed.flow.instant.m3_hr-explicit.routing.run$routed.flow.kinematic.m3_hr
#x11()

# Function to read in the results of a run, run explicit routing code, output the fdc of reaches, estimate the 'threshold flow',
# estimate the time that the threshold flow is exceeded, determine the timing of the 'on' and 'off' of the sensor...
# Inputs: qbf/Qin for each HRU for a behavioral simulation, FR.Spatial configuration, sensor inputs
Headwater.evaluation.dynamic <- function(qin.file, FR.Spatial, TS.logger.clean.all, TS.zoo.all, logger.out.percent, Reach.identifiers, Reach.identifiers.second,model.timestep, result.optim.run, time.qbf, ...) {
  # Read in explicit HRU results for fluxes and storages from the behavioral sets
  # NOTE WE MUST SET THE DIRECTORY TO THE PROPER CALIBRATION TEST FOLDER, AS OF 7/6/2021 IT IS CalibrationTest13
  fluxes.stores.dir <- paste0(getwd(),'/CalibrationOutput/fluxes_stores/')
  #setwd('C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/2 DYNAMIC TOPMODEL ANALYSIS/1 FR Test/Calibration/CalibrationTest13/fluxes_stores')
  qin.run <- read.csv(paste0(fluxes.stores.dir,'fluxqin',qin.file,'.csv'))                           # Read qin 
  qin.run <- xts(qin.run[,-1],time.qbf)                                            # Convert to xts
  qbf.run <- read.csv(paste0(fluxes.stores.dir,'fluxqbf',qin.file,'.csv'))                           # Read qbf
  qbf.run <- xts(qbf.run[,-1],time.qbf)                                            # Convert to xts
  ae.run <- read.csv(paste0(fluxes.stores.dir,'fluxae',qin.file,'.csv'))                             # Read ae
  ae.run <- xts(ae.run[,-1],time.qbf)                                              # Convert to xts
  rain.run <- read.csv(paste0(fluxes.stores.dir,'fluxrain',qin.file,'.csv'))                         # Read rain
  rain.run <- xts(rain.run[,-1],time.qbf)                                          # Convert to xts
  qof.run <- read.csv(paste0(fluxes.stores.dir,'fluxqof',qin.file,'.csv'))                           # Read qof
  qof.run <- xts(qof.run[,-1],time.qbf)                                            # Convert to xts
  Qsim.run <- read.csv(paste0(fluxes.stores.dir,'qsim',qin.file,'.csv'))                             # Read Qsim
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
  
  # Convert Q from mm/s to mm/day
  Q <- Q*60*60*24                                                                  # This is important for some of the empirical calculations we'll later need to do
  
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
      
      # Determine the subsurface flow capacity of the valley below the reach. This equation is proposed in Godsey and Kirchner (2018) and Prancevic and Kirchner (2014)
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
  #data.directory <- paste('C:/Users/david/OneDrive/Desktop/EPA/',        # Define the directory string
  #                        'EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/',                        # I'm breaking up the paste for readability
  #                        '3 DATA/2-PROJECT1-ROBINSON-FOREST/',         # This isn't the working directory where the PSO items will be stored
  #                        '3-UK-ROBINSONFOREST-DATA/TIME_SERIES/',      # This is just the WD for reading in input files
  #                        sep='')  
  #setwd(data.directory)                                                 # Set the working directory
  
  Met.data.in <- read_csv('HydroInputData/CWS_Rain_clean_hourly.csv')                                    # Read in the data with forcing info
  Q.data.in <- read_csv('HydroInputData/Qobs2000_2015.csv')                              # Read in the observed Q data for falling rock (from Chris/Science base)
  PET.data.in <- read_csv('HydroInputData/WATER_pub_2013May2_2.csv')                   # Read in the PETdata
  Temp.data.in <- read_csv('HydroInputData/Air_Temp.csv')                              # Read in the temp data
  Day.data.in <- read_csv('HydroInputData/FR_precip.csv')                              # Read in the daylight hours
  
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
    
    
    pex <- in.rain-in.ae
    pex[pex<0] <- 0
    
    flows.explicit.qbf <- rep(0,nrow(explicit.weights.full))  
                            
    
    flows.explicit.qbf[(num.reach+1):length(flows.explicit.qbf)] <- in.qbf[2:length(in.qbf)] 
    
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
  
  ## set WD
  #setwd('C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/4 ANALYSIS/1 GIS ANALYSIS/')
  
  # Read in 1.5 m DEM, Flow length (unused), and soil Raster
  fr1mDEM <- raster('SpatialInputData/FR1meterDEM.tif')                                                # Read in the 1 m FR DEM
  
  # NOTE AS OF 7/6/2021 - there has been an update to the raster package such and now the raster function appears to no longer be reading the crs (coordinate system) of the raster correctly. 
  # The Crs is: +proj=lcc +lat_0=36.3333333333333 +lon_0=-85.75 +lat_1=37.0833333333333 +lat_2=38.6666666666667 +x_0=1500000 +y_0=999999.9998984 +datum=NAD83 +units=us-ft +no_defs 
  
  crs(fr1mDEM) <- '+proj=lcc +lat_0=36.3333333333333 +lon_0=-85.75 +lat_1=37.0833333333333 +lat_2=38.6666666666667 +x_0=1500000 +y_0=999999.9998984 +datum=NAD83 +units=us-ft +no_defs'
  
  #fr1mFlowLen <- raster('Rasters/frFlowLen1m2.tif')                                           # Read in the flow length raster
  FRSoils1m <- raster('SpatialInputData/FRSoils1mMASK.tif')                                            # Read in the soils raster
  crs(FRSoils1m) <- '+proj=lcc +lat_0=36.3333333333333 +lon_0=-85.75 +lat_1=37.0833333333333 +lat_2=38.6666666666667 +x_0=1500000 +y_0=999999.9998984 +datum=NAD83 +units=us-ft +no_defs'
  
  # Read in stream network
  FR_drn <- shapefile('SpatialInputData/FR4000StreamNet.shp')
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
  #if (rowSums(FR.disc.explicit$weights[,1:34])[35:nrow(FR.disc.explicit$weights)] == 0 &&
  #    FR.disc$weights[3:nrow(FR.disc$weights),1] > 0) {
   # warning('There is zero flow in an explicit reach but positive flow in a lumped')
  #}
  
  ##for (row in 1:(nrow(FR.disc$weights)-2)) {                                                 # here I am assuming that the ratios from HRUs into channels is correct for the lumped. So I am making sure that the ratios from HRUs in the explicit discretization matches that of the first
  #  FR.disc.explicit$weights[(34+row),1:34] <- FR.disc.explicit$weights[(34+row),1:34]*
  #    correct.ratio[row]
  #  if (is.na(FR.disc.explicit$weights[(34+row),1:34])) {
  #    FR.disc.explicit$weights[(34+row),1:34] <- 0
  #  }
  #}
  
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
