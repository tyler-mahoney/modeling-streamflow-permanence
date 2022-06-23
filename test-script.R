# Test script to upload to Github 

library(ggplot2)
rand.df <- data.frame('rand.one'=runif(10000,0,1),'rand.two'=runif(10000,0,1))  # Create df of random numbers

ggplot(data=rand.df,aes(x=rand.one,y=rand.two)) + geom_point()                  # Plot the x vs y random numbers
