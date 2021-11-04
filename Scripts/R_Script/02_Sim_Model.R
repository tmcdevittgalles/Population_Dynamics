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
library(gamm4)

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


weather.df <- read.csv("Data/simple_weather.csv")



weather.df$Julian <- weather.df$DOY+(weather.df$Year-min(weather.df$Year))*365


### Temperature
gam.temp <- gam( TMIN ~ s(Julian,k=25, bs= "cr"),
                 data=weather.df, family="gaussian")


summary(gam.temp)
## lets try and plot this with new data

pred.df <- dplyr::select( weather.df, c("Julian"))

dum.df <-  unique(pred.df)

weather.df$STemp<- predict(gam.temp, newdata = dum.df)


weather.df  %>% 
  ggplot(  aes(x=Julian, y=STemp) ) + 
  geom_point( aes(x= Julian,y=TMIN),size=1)+
  geom_line(size=2,alpha=.75)+ theme_classic()

weather.df$mPPT <- caTools::runmean(weather.df$PPT,30)

plot(weather.df$Julian,weather.df$mPPT, type = "l")
plot(weather.df$Julian,weather.df$TMIN, type = "l")
plot(weather.df$Julian,weather.df$STemp, type = "l")

X <- model.matrix(~scale(weather.df$STemp)*scale(weather.df$mPPT))

y.df <- data.frame( Julian = weather.df$Julian )

y.df$Count <- NA

y.df$Offset <- rbeta(nrow(y.df),9,3)   

## establising known parameters

betas <- c(-5, 2, .5,-1)

rho <- .8

sigma <- 4


y.df <- y.df %>% arrange( Julian)

y.df$Count[1] <-0

r <- rep(NA, nrow(y.df)-1 )

# simulating data ?

for(i in 2:nrow(y.df)){
  
r[i] <- (sum(X[i-1,]*betas))+((y.df$Count[i-1])*rho)
  
  y.df$Count[i] <-  ( rnorm( 1,(r[i])*y.df$Offset[i] ,sigma) )
}

y.df %>% ggplot(aes(x=Julian, y= ((Count) ))) + geom_point()


y.df %>% ggplot(aes(x=Julian, y= log10(exp(Count)+1 ))) + geom_point()

y.df$Julian <- 1:nrow(y.df)

y.df$rCount <- as.integer( round( exp(y.df$Count), 0 ))

y.df %>% ggplot(aes(x=Julian, y= log10(rCount+1 ))) + geom_point()+
  theme_classic()+ylab("Simulated mosquito count")+ xlab("Julian date")+ 
  
  theme( legend.key.size = unit(1.5, "cm"),
         legend.title =element_text(size=14,margin = margin(r =10,unit = "pt")),
         legend.text=element_text(size=14,margin = margin(r =10, unit = "pt")), 
         legend.position = "top",
         axis.line.x = element_line(color="black") ,
         axis.ticks.y = element_line(color="black"),
         axis.ticks.x = element_line(color="black"),
         axis.title.x = element_text(size = rel(1.8)),
         axis.text.x  = element_text(vjust=0.5, color = "blac`k",size=14),
         axis.text.y  = element_text(vjust=0.5,color = "black",size=14),
         axis.title.y = element_text(size = rel(1.8), angle = 90) ,
         strip.text.x = element_text(size=20) )

## testing to see if i can get these values back

stan_d <- list( N = nrow(y.df), P = ncol(X),
                Time = nrow(X),
                Julian = y.df$Julian,
                X= X,
                Y= y.df$rCount, offset = y.df$Offset )

ar_output3a <- stan( 'Scripts/Stan/Ar_take_3.stan', 
                    data=stan_d, iter = 2000,
                    control = list(max_treedepth = 10))



print(ar_output3a, pars = c("rho", 'beta', "sigma"))
traceplot(ar_output3)


saveRDS(ar_output3a ,"bayes_sim.rds")

saveRDS(y.df,"sim_data.rds")

ar_output3a <- readRDS("bayes_sim.rds")
 
 
post <- extract(ar_output3a)

## plotting rho values

rho.df <- as.data.frame(post$rho)

colnames(rho.df) <- "Rho"

sig.df <- as.data.frame(post$sigma)

colnames(sig.df) <- "Sigma"
## plotting rho values

beta.df <- as.data.frame(post$beta)


colnames(beta.df) <- c("Intercept", "TMin", "PPT", "TMin:PPT")

beta.df <- cbind.data.frame(beta.df, rho.df)

beta.df <- cbind.data.frame(beta.df, sig.df)



beta.df <- tidyr::pivot_longer(beta.df, cols= 1:6, names_to="Parameters",
                               values_to = "Estimate")

beta.df$Parameters <- factor( beta.df$Parameters, level=c("Rho", "Intercept", "TMin", "PPT",
                                                          "TMin:PPT", "Sigma"))

known.df <- data.frame( Pars = as.factor(c("Rho", "Intercept", "TMin", "PPT",
                                           "TMin:PPT", "Sigma")),
                        Value = as.numeric(c(.8,-5,2,.5,-1,4)) )

ggplot(beta.df, aes(x=Parameters, y=Estimate))+  geom_violin(alpha=.5,fill="gray")+
  geom_boxplot(color="black", width=.05, outlier.size = 0) +
  geom_hline(yintercept = 0, size=1 ) +theme_classic()+
  ylab("Estimate")+ xlab("Parameters")+
  geom_point(data=known.df, aes(x=Pars, y = Value),size=3, color="red")+
  theme( legend.key.size = unit(1.5, "cm"),
         legend.title =element_text(size=14,margin = margin(r =10,unit = "pt")),
         legend.text=element_text(size=14,margin = margin(r =10, unit = "pt")), 
         legend.position = "top",
         axis.line.x = element_line(color="black") ,
         axis.ticks.y = element_line(color="black"),
         axis.ticks.x = element_line(color="black"),
         axis.title.x = element_text(size = rel(1.8)),
         axis.text.x  = element_text(vjust=0.5, color = "black",size=14),
         axis.text.y  = element_text(vjust=0.5,color = "black",size=14),
         axis.title.y = element_text(size = rel(1.8), angle = 90) ,
         strip.text.x = element_text(size=20) )


### Using less data points ( subsetting to what we survey in the real data)


full.df <- read.csv("Data/simple_df.csv")

dim(full.df) ## 1923 X 24

names(full.df)

## Data structure
## 
##  Species: Aedes vexans
##  Domain : D05 - great lakes
##  Site: UNDE
##  Plots: 11
##  Years: 5 , 2014,2016,2017, 2018, 2019
##  Unique sampling events : 1923

## adding julian date to data frame

full.df$Julian <- (full.df$DOY+(full.df$Year-min(full.df$Year))*365) +1

## subsetting the simulated data set 

subY.df <- y.df %>% filter( Julian %in% full.df$Julian)


subY.df %>% ggplot(aes(x=Julian, y= (rCount ))) + geom_point()


stan_d <- list( N = nrow(subY.df), P = ncol(X),
                Time = nrow(X),
                Julian = subY.df$Julian,
                X= X,
                Y= subY.df$rCount, offset = subY.df$Offset )

ar_output3s <- stan( 'Scripts/Stan/Ar_take_3.stan', 
                     data=stan_d, iter = 5000,
                     control = list(max_treedepth = 20))


print(ar_output3s, pars = c("rho", 'beta', "sigma"))
traceplot(ar_output3s)


post <- extract(ar_output3s)


#saveRDS(ar_output3 ,"bayes_sim.rds")


## plotting rho values

rho.df <- as.data.frame(post$rho)

colnames(rho.df) <- "Rho"

sig.df <- as.data.frame(post$sigma)

colnames(sig.df) <- "Sigma"
## plotting rho values

beta.df <- as.data.frame(post$beta)


colnames(beta.df) <- c("Intercept", "TMin", "PPT", "TMin:PPT")

beta.df <- cbind.data.frame(beta.df, rho.df)

beta.df <- cbind.data.frame(beta.df, sig.df)



beta.df <- tidyr::pivot_longer(beta.df, cols= 1:6, names_to="Parameters",
                               values_to = "Estimate")

beta.df$Parameters <- factor( beta.df$Parameters, level=c("Rho", "Intercept", "TMin", "PPT",
                                                          "TMin:PPT", "Sigma"))

known.df <- data.frame( Pars = as.factor(c("Rho", "Intercept", "TMin", "PPT",
                                           "TMin:PPT", "Sigma")),
                        Value = as.numeric(c(.8,-1,1,-.5,.15,1.5)) )

ggplot(beta.df, aes(x=Parameters, y=Estimate))+  geom_violin(alpha=.5,fill="gray")+
  geom_boxplot(color="black", width=.05, outlier.size = 0) +
  geom_hline(yintercept = 0, size=1 ) +theme_classic()+
  ylab("Estimate")+ xlab("Parameters")+
  geom_point(data=known.df, aes(x=Pars, y = Value),size=3, color="red")+
  theme( legend.key.size = unit(1.5, "cm"),
         legend.title =element_text(size=14,margin = margin(r =10,unit = "pt")),
         legend.text=element_text(size=14,margin = margin(r =10, unit = "pt")), 
         legend.position = "top",
         axis.line.x = element_line(color="black") ,
         axis.ticks.y = element_line(color="black"),
         axis.ticks.x = element_line(color="black"),
         axis.title.x = element_text(size = rel(1.8)),
         axis.text.x  = element_text(vjust=0.5, color = "black",size=14),
         axis.text.y  = element_text(vjust=0.5,color = "black",size=14),
         axis.title.y = element_text(size = rel(1.8), angle = 90) ,
         strip.text.x = element_text(size=20) )




plot.df <- as.data.frame(colMedians(as.matrix(post$m)))

test <- colQuantiles(as.matrix(post$m), probs=c(0.025,0.975))

plot.df$Q_2.5  <- test[,1]

plot.df$Q_97.5  <- test[,2]
plot.df$Julian <- 1:nrow(plot.df)

colnames(plot.df)[1] <- "Count"


plot.df %>% 
  ggplot( aes(x=Julian, y=log(exp(Count)+1))) +
  geom_ribbon(aes(x=Julian, ymin= log(exp(Q_2.5)+1), 
                  ymax= log(exp(Q_97.5)+1) ),alpha=.2 )+
  geom_line(alpha=.85, size=1) +
  geom_point(data=y.df, aes(x=Julian, y=log(( exp(Count) *Offset) +1)),
                              alpha=.25)+
  ylab("log(Predicted mosquito count)") + 
  xlab("Julian date")+ theme_classic()+
  theme( legend.key.size = unit(1.5, "cm"),
         legend.title =element_text(size=14,margin = margin(r =10,unit = "pt")),
         legend.text=element_text(size=14,margin = margin(r =10, unit = "pt")), 
         legend.position = "top",
         axis.line.x = element_line(color="black") ,
         axis.ticks.y = element_line(color="black"),
         axis.ticks.x = element_line(color="black"),
         axis.title.x = element_text(size = rel(1.8)),
         axis.text.x  = element_text(vjust=0.5, color = "black",size=14),
         axis.text.y  = element_text(vjust=0.5,color = "black",size=14),
         axis.title.y = element_text(size = rel(1.8), angle = 90) ,
         strip.text.x = element_text(size=20) )

