###### Powell Center: Phenological patterns of mosquitoes #######

# Travis McDevitt-Galles
# 08/24/2021
# title: 02_Sim_Model

# Simulating data and parameter values to validate current model structure
# before adding complexity or additional levels to the model

library(dplyr)
library(ggplot2)
library(patchwork)
library(rstan)
library(rstanarm)
library(matrixStats)

#Set working directory
setwd("C:/Users/tmcdevitt-galles/Documents/Population_dynamics")

# simulating population data for a single plot overtime. 
# 
# Data frame should represent 5 years of sampling with ~ 25 sampling events
# per year that occur during the summer portion. Count data will have a fixed
# offset just for consistence. There are 6 parameters we will specify and will
# want the model to extract, 
#   1) Rho: how much is the current count dependent on the previous count
#   2) Beta[1]: the intercet
#   3) Beta[2]: the effect of temperature on population
#   4) Beta[3]: the effect of precipitation on population
#   5) Beta[4]: the interactive effects of temp and precip on population
#   6) Sigma: error in population model 
#   
#   We will not simulate covariates (Weather data) but instead will use real 
#   data 


### Setting the values:
###   Rho = .9
###   Beta[1] = -3
###   Beta[2] = 2
###   Beta[3] = .5
###   Beta[4] = -.5
###   Sigma = 1.5


## Lets make it simple by using the total number of sampling events 

real.df <- read.csv("Data/simple_weather.csv")

# adding julian dates to dataframe
# 
real.df$Julian <- real.df$DOY+(real.df$Year-min(real.df$Year))*365

covar.df <- real.df

## Plotting the data
covar.df %>%  ggplot(aes(x=Julian, y=TMIN) ) + geom_point()

covar.df %>%  ggplot(aes(x=Julian, y=PPT) ) + geom_point()

# creating an empty vector to store simulated data

y.df <- data.frame( Julian = covar.df$Julian )

y.df$Count <- NA

y.df$Offset <- rbeta(nrow(y.df),9,3)   

X <- model.matrix(y.df$Offset ~ (covar.df$TMIN) *
                  (covar.df$PPT) )

## establising known parameters

betas <- c(-3,.2,0,0)

rho <- .8

sigma <- 1.5


y.df <- y.df %>% arrange( Julian)

y.df$Count[1] <-0

r <- rep(NA, nrow(y.df)-1 )

# simulating data ?

for(i in 2:nrow(y.df)){
  
r[i] <- (sum(X[i-1,]*betas))*(y.df$Count[i-1])+((y.df$Count[i-1])*rho)
  
  y.df$Count[i] <- rlnorm( 1,(r[i]) ,1.5)
}

y.df %>% ggplot(aes(x=Julian, y= exp((Count) ))) + geom_point()


y.df %>% ggplot(aes(x=Julian, y= log10((Count)+1 ))) + geom_point()


Plot_Samples <- nrow(y.df)
## testing to see if i can get these values back


stan_d <- list( N = nrow(y.df), P = ncol(X),
                X= X, G = length(Plot_Samples),
                group_samp = c(Plot_Samples),
                Y= y.df$Count, offset = y.df$Offset )

ar_output3 <- stan( 'Scripts/Stan/Ar_take2.stan', 
                    data=stan_d, iter = 5000,
                    control = list(max_treedepth = 20))