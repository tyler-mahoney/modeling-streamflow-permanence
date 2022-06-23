# Script to read in 15-minute time series data from Chris Barton and manipulate the timeseries, as needed 
# Result of the script is to have serially complete time series data for 15-min time intervals
# from 2003-2004 and 2011-2012

# read in raw csv of 2011-2012 data. 
# Note that the dynatopmodel brompton$pe dataset skips the 2:00-3:00 hour on 3/11/2012 for daylight savings. 
# The Code below will treat the change of from EST to EDT as follows: the hour on 3/13/2011 from 2:00-3:00 that would
# otherwise be deleted here will just be rerun
# Same with 2012. 
library(zoo)
library(xts)
setwd('C:/Users/DMAHONEY/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/Tyler/EPA/3 PROJECT 1 KENTUCKY HEADWATER STREAMS/3 DATA/2-PROJECT1-ROBINSON-FOREST/3-UK-ROBINSONFOREST-DATA/TIME_SERIES')
# Read in 2010 data and convert to timeseries. Note that the date in the csv must be the form: yyyy/mm/dd hr:mm:ss
FR2010 <- read.csv('FR2010.csv')
FR2010$Date.and.Time <- as.POSIXct(FR2010$Date.and.Time)
FR2010.discharge <- data.frame('Q_cms'=FR2010$cms)
rownames(FR2010.discharge) <- FR2010$Date.and.Time
FR2010.discharge <- as.xts(FR2010.discharge)

# Read in 2011 data
FR2011 <- read.csv('FR2011.csv')
FR2011$Date.and.Time <- as.POSIXct(FR2011$Date.and.Time, tz = "EST")
correctTS2011 <- seq.POSIXt(from=ISOdatetime(2011,01,01,00,0,0,tz="EST"),to=ISOdate(2011,12,31,23,45,0,tz="EST"),by="15 min")
Q.cms <- vector(length=length(correctTS2011))
FR2011.discharge <- data.frame(Q.cms)

idxmatch <- match(correctTS2011,FR2011$Date.and.Time)

for (i in 1:length(correctTS2011))
{
  FR2011.discharge$Q.cms[i] <- as.numeric(FR2011$cms[idxmatch[i]])
  if (is.na(FR2011.discharge$Q.cms[i]))
  {
    FR2011.discharge$Q.cms[i] = as.numeric(FR2011.discharge$Q.cms[i-1])
    # Here I am simlply stating that if htere is no discharge for a time step, then I will set the discharge
    # to that of the previous time step
  }
}

FR2011.discharge <- xts(x=FR2011.discharge, order.by=correctTS2011)

#check if the spacing for the timeseries is equal
dt.ser <- diff(as.numeric(index(FR2011.discharge)))
if(!all(dt.ser[]==dt.ser[1]))
{
  warning('theres a problem with the time series')
}

which(dt.ser[]!=900)
View(FR2011.discharge)


# Read in 2012 data
FR2012 <- read.csv('FR2012.csv')
FR2012$Date.and.Time <- as.POSIXct(FR2012$Date.and.Time, tz = 'EST')
correctTS2012 <- seq.POSIXt(from=ISOdatetime(2012,01,01,00,0,0,tz="EST"),to=ISOdate(2012,12,06,11,00,0,tz="EST"),by="15 min")
Q.cms <- vector(length=length(correctTS2012))
FR2012.discharge <- data.frame(Q.cms)

idxmatch <- match(correctTS2012,FR2012$Date.and.Time)

for (i in 1:length(correctTS2012))
{
  FR2012.discharge$Q.cms[i] <- FR2012$cms[idxmatch[i]]
  if (is.na(FR2012.discharge$Q.cms[i]))
  {
    FR2012.discharge$Q.cms[i] = FR2012.discharge$Q.cms[i-1] 
    # Here I am simlply stating that if htere is no discharge for a time step, then I will set the discharge
    # to that of the previous time step
  }
}

FR2012.discharge <- xts(x=FR2012.discharge, order.by=correctTS2012)

#check if the spacing for the timeseries is equal
dt.ser <- diff(as.numeric(index(FR2012.discharge)))
if(!all(dt.ser[]==dt.ser[1]))
{
  warning('theres a problem with the time series')
}

which(dt.ser[]!=900)
View(FR2012.discharge)


# Read in 2002 data
FR2002 <- read.csv('FR2002.csv') # Note, the input data format should be "%YYYY-%m-%d $H:$mm$ss
# Read in the timeseries and convert the string to a POSIX
FR2002$Time <- as.POSIXct(FR2002$Time, tz = "EST")
# create a uniformly spaced time series
correctTS2002 <- seq.POSIXt(from=ISOdatetime(2002,11,01,00,0,0,tz="EST"),to=ISOdate(2002,12,31,23,45,0,tz="EST"),by="15 min")
Q.cms <- vector(length=length(correctTS2002))
FR2002.discharge <- data.frame(Q.cms)

idxmatch <- match(correctTS2002,FR2002$Time)

for (i in 1:length(correctTS2002))
{
  FR2002.discharge$Q.cms[i] <- FR2002$cms[idxmatch[i]]
  if (is.na(FR2002.discharge$Q.cms[i]))
  {
    FR2002.discharge$Q.cms[i] = FR2002.discharge$Q.cms[i-1] 
    # Here I am simlply stating that if htere is no discharge for a time step, then I will set the discharge
    # to that of the previous time step
  }
}

FR2002.discharge <- xts(x=FR2002.discharge, order.by=correctTS2002)

#check if the spacing for the timeseries is equal
dt.ser <- diff(as.numeric(index(FR2002.discharge)))
if(!all(dt.ser[]==dt.ser[1]))
{
  warning('theres a problem with the time series')
}

which(dt.ser[]!=900)
View(FR2002.discharge)



# Read in 2003 data
FR2003 <- read.csv('FR2003.csv') # Note, the input data format should be "%YYYY-%m-%d $H:$mm$ss
# Read in the timeseries and convert the string to a POSIX
FR2003$Time <- as.POSIXct(FR2003$Time, tz = "EST")
# create a uniformly spaced time series
correctTS2003 <- seq.POSIXt(from=ISOdatetime(2003,01,01,00,0,0,tz="EST"),to=ISOdate(2003,12,31,23,45,0,tz="EST"),by="15 min")
Q.cms <- vector(length=length(correctTS2003))
FR2003.discharge <- data.frame(Q.cms)

idxmatch <- match(correctTS2003,FR2003$Time)

for (i in 1:length(correctTS2003))
{
  FR2003.discharge$Q.cms[i] <- FR2003$cms[idxmatch[i]]
  if (is.na(FR2003.discharge$Q.cms[i]))
  {
    FR2003.discharge$Q.cms[i] = FR2003.discharge$Q.cms[i-1] 
    # Here I am simlply stating that if htere is no discharge for a time step, then I will set the discharge
    # to that of the previous time step
  }
}

FR2003.discharge <- xts(x=FR2003.discharge, order.by=correctTS2003)

#check if the spacing for the timeseries is equal
dt.ser <- diff(as.numeric(index(FR2003.discharge)))
if(!all(dt.ser[]==dt.ser[1]))
{
  warning('theres a problem with the time series')
}

which(dt.ser[]!=900)
View(FR2003.discharge)




# Read in 2004 data
FR2004 <- read.csv('FR2004.csv') # Note, the input data format should be "%YYYY-%m-%d $H:$mm$ss
# Read in the timeseries and convert the string to a POSIX
FR2004$Time <- as.POSIXct(FR2004$Time, tz = "EST")
# create a uniformly spaced time series
correctTS2004 <- seq.POSIXt(from=ISOdatetime(2004,01,01,00,0,0,tz="EST"),to=ISOdate(2004,11,01,23,45,0,tz="EST"),by="15 min")
Q.cms <- vector(length=length(correctTS2004))
FR2004.discharge <- data.frame(Q.cms)

idxmatch <- match(correctTS2004,FR2004$Time)

for (i in 1:length(correctTS2004))
{
  FR2004.discharge$Q.cms[i] <- FR2004$cms[idxmatch[i]]
  if (is.na(FR2004.discharge$Q.cms[i]))
  {
    FR2004.discharge$Q.cms[i] = FR2004.discharge$Q.cms[i-1] 
    # Here I am simlply stating that if htere is no discharge for a time step, then I will set the discharge
    # to that of the previous time step
  }
}

FR2004.discharge <- xts(x=FR2004.discharge, order.by=correctTS2004)

#check if the spacing for the timeseries is equal
dt.ser <- diff(as.numeric(index(FR2004.discharge)))
if(!all(dt.ser[]==dt.ser[1]))
{
  warning('theres a problem with the time series')
}

which(dt.ser[]!=900)
View(FR2004.discharge)

FR2002.2004.discharge = rbind(FR2002.discharge,FR2003.discharge,FR2004.discharge)
FR2010.2012.discharge = rbind(FR2010.discharge,FR2011.discharge,FR2012.discharge)
typeof(FR2010.discharge)
typeof(FR2011.discharge)
typeof(FR2012.discharge)
