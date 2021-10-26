###### Powell Center: Phenological patterns of mosquitoes #######

# Travis McDevitt-Galles
# 08/24/2021
# title: 04_Pop_Model_Novel_Species

# Recreating Mevin's modeling framework to see if i can retain

library(dplyr)

library(ggplot2)
library(patchwork)
library(rstan)
library(rstanarm)
library(matrixStats)
library(gamm4)
library(here)

#Set working directory
setwd("C:/Users/tmcdevitt-galles/Documents/Population_dynamics")


### Loading in full dataset
load("Data/Mosquito_Data_Clean.Rda" ) 


full.df <- complete.df %>%  filter(Site=="UNDE")

full.df <- full.df %>% filter(  SciName == "Coquillettidia perturbans" |
                                SciName == "Aedes vexans")


dim(full.df) ## 1923 X 24

names(full.df)

## Data structure
## 
##  Species: Aedes vexans
##  Domain : D01 - Northeast
##  Site: UNDE
##  Plots: 11
##  Years: 5 , 2014,2016,2017, 2018, 2019
##  Unique sampling events : 1649
## adding julian date to data frame

full.df$Year <- as.integer(full.df$Year)

full.df$Julian <- julian( full.df$Date , 
        origin = as.Date("2014-01-01"))

full.df$Julian <- full.df$Julian +1

stan.df <- full.df %>% filter(SubsetWeight >0)

stan.df <- stan.df %>% group_by(SciName, Julian, DOY, Plot, Domain, Year) %>% 
  summarise(
    Count = sum(Count),
    TrapHours = sum(TrapHours),
    SubsetWeight = sum( SubsetWeight),
    TotalWeight = sum( TotalWeight )
  )


stan.df$Offset <- (stan.df$TrapHours/24) * 
  (stan.df$SubsetWeight/stan.df$TotalWeight)



stan.df %>% 
  ggplot(aes(x=DOY,y= Count*Offset  ,color=as.factor(Year)) )+ geom_point()

stan.df %>% 
  ggplot(aes(x=Julian,y= (Count*Offset), color=SciName ))+
  geom_point()+     scale_color_brewer(palette="Set1", name="SciName")  +
  theme_classic()+ylab("Mosquito Count")+ xlab("Julian date")+ 
 
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

#### Add in the weather data


load("Data/DailyPrismMod.Rda" ) 

site.df <- unique(select(ungroup(complete.df), c("Site","Plot")))


contigus.df <- right_join(contigus.df,site.df, by="Plot")


weather.df <- contigus.df %>% filter(Site=="UNDE")
weather.df$Year <- as.integer(weather.df$Year)

weather.df$Julian <- julian( weather.df$Date , 
                                    origin = as.Date("2014-01-01"))

weather.df <- weather.df %>% group_by(Site, Year, Julian, Date) %>% 
  summarize(
    TMIN = mean(TMIN),
    PPT = mean(PPT)
  )



### Temperature
gam.temp <- gam( TMIN ~ s(Julian,k=25, bs= "cr"),
                 data=weather.df, family="gaussian")


summary(gam.temp)
## lets try and plot this with new data

pred.df <- dplyr::select( weather.df, c("Julian"))

dum.df <-  unique(pred.df)



weather.df$STemp<- as.numeric(predict(gam.temp, newdata = dum.df))


weather.df  %>% 
  ggplot(  aes(x=Julian, y=STemp) ) + 
  geom_point( aes(x= Julian,y=TMIN),size=1)+
  geom_line(size=2,alpha=.75)+ theme_classic()

weather.df$mPPT <- caTools::runmean(weather.df$PPT,30)

plot(weather.df$Julian,weather.df$mPPT, type = "l")
plot(weather.df$Julian,weather.df$TMIN, type = "l")
plot(weather.df$Julian,weather.df$STemp, type = "l")

weather.df <- weather.df %>% filter( Year > 2013)

weather.df$Julian <- julian( weather.df$Date , 
                             origin = as.Date("2014-01-01"))

weather.df$Julian <- weather.df$Julian +1

weather.df$Track <- 1:nrow(weather.df)

weather.df1 <- weather.df

weather.df1$SciName <- "Aedes vexans"

weather.df2 <- weather.df

weather.df2$SciName <- "Coquillettidia perturbans"

weather.df <- rbind.data.frame(weather.df1, weather.df2)


stan.df %>% 
ggplot(aes(x=Julian,y= log10(Count*Offset +1) ))+
  geom_point(aes(color=SciName))+     scale_color_brewer(palette="Set1", name="SciName")  +
  theme_classic()+ylab("Mosquito Count")+ xlab("Julian date")+ 
  geom_line(data=weather.df, aes(x= Julian, y= STemp/5),size=2)+
  
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

Weather <- model.matrix(~scale(weather.df$STemp)*scale(weather.df$mPPT)*weather.df$SciName)


spp <- model.matrix(~0 + stan.df$SciName )

stan_d <- list( N = nrow(stan.df), P = ncol(Weather),
                Time = max(weather.df$Julian), nSpp = ncol(spp),
                Julian = stan.df$Julian,
                X= Weather,
                Spp = spp,
                Y= stan.df$Count, offset = stan.df$Offset )

ar_output5 <- stan( 'Scripts/Stan/Two_Spp_SSM2.stan', 
                    data=stan_d, iter = 2000,
                    control = list(max_treedepth = 15))

## Start time 5:06

 print(ar_output5, pars = c("beta", 'rho'))
traceplot(ar_output5, pars =c("beta"))

post <- extract(ar_output5)

plot(x=1:nrow(weather.df),y =exp(colMedians(as.matrix(post$m[,,2]))), type= "l")
points(x=stan.df$Julian, y= log10((stan.df$Count*stan.df$Offset+1)))

saveRDS(ar_output5 ,"two_spp2.rds")

ar_output5 <- readRDS("bayes_back.rds")

## plotting rho values

rho.df <- as.data.frame(post$rho)

colnames(rho.df) <- "Rho"
## plotting rho values

beta.df <- as.data.frame(post$beta)


colnames(beta.df) <- c("Intercept",  "Temp", "PPT","C.Per",
                       "Temp:PPT", "Temp:C.Per", 
                       "PPT:C.Per", "ThreeWay")

beta.df <- cbind.data.frame(beta.df, rho.df)

beta.df <- tidyr::pivot_longer(beta.df, cols= 1:8, names_to="Parameters",
                               values_to = "Estimate")



ggplot(beta.df, aes(x=Parameters, y=Estimate))+  geom_violin(alpha=.5,fill="gray")+
  geom_boxplot(color="black", width=.05, outlier.size = 0) +
  geom_hline(yintercept = 0, size=1 ) +theme_classic()+
  ylab("Estimate")+ xlab("Parameters")+
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



plot1.df <- beta.df %>% filter( Parameters =="PPT" |
                                  Parameters == "TMin"|
                                  Parameters == "TMin:PPT")



ggplot(plot1.df, aes(x=Parameters, y=Estimate))+  geom_violin(alpha=.5,fill="gray")+
  geom_boxplot(color="black", width=.05, outlier.size = 0) +
  geom_hline(yintercept = 0, size=1 ) +theme_classic()+
  ylab("Estimate")+ xlab("Parameters")+
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



plot.df <- as.data.frame(colMedians(as.matrix(post$m[,,2])))

test <- colQuantiles(as.matrix(post$m[,,2]), probs=c(0.025,0.975))

plot.df$Q_2.5  <- test[,1]

plot.df$Q_97.5  <- test[,2]
plot.df$Julian <- 1:nrow(plot.df)

colnames(plot.df)[1] <- "Count"

plot.df$SciName <- "Coquillettidia perturbans"


plot1.df <- as.data.frame(colMedians(as.matrix(post$m[,,1])))

test <- colQuantiles(as.matrix(post$m[,,1]), probs=c(0.025,0.975))

plot1.df$Q_2.5  <- test[,1]

plot1.df$Q_97.5  <- test[,2]
plot1.df$Julian <- 1:nrow(plot1.df)

colnames(plot1.df)[1] <- "Count"

plot1.df$SciName <- "Aedes vexans"

plot.df <- rbind.data.frame(plot.df, plot1.df)

plot.df %>% 
  ggplot( aes(x=Julian, y=log(exp(Count)+1))) +
  geom_ribbon(aes(x=Julian, ymin= log(exp(Q_2.5)+1), 
               ymax= log(exp(Q_97.5)+1),fill=SciName ),alpha=.2 )+
  geom_line(aes(color=SciName),alpha=.85, size=2)+
  geom_point(data=stan.df, aes(x=Julian, y=log(( Count *Offset) +1),
                               color=as.factor(SciName)), alpha=.25)+
  scale_color_brewer(palette="Set1", name="Species")+
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

#### Plotting growth rate #####
grow.df <- as.data.frame(colMedians(as.matrix(post$r[,,2])))

test <- colQuantiles(as.matrix(post$r[,,2]), probs=c(0.025,0.975))

grow.df$Q_2.5  <- test[,1]

grow.df$Q_97.5  <- test[,2]
grow.df$Julian <- 1:nrow(grow.df)

colnames(grow.df)[1] <- "Count"

grow.df$SciName <- "Coquillettidia perturbans"


grow1.df <- as.data.frame(colMedians(as.matrix(post$r[,,1])))

test <- colQuantiles(as.matrix(post$r[,,1]), probs=c(0.025,0.975))

grow1.df$Q_2.5  <- test[,1]

grow1.df$Q_97.5  <- test[,2]
grow1.df$Julian <- 1:nrow(grow1.df)

colnames(grow1.df)[1] <- "Count"

grow1.df$SciName <- "Aedes vexans"

grow.df <- rbind.data.frame(grow.df, grow1.df)

grow.df %>% 
  ggplot( aes(x=Julian, y=exp(Count))) +
  geom_hline(yintercept = 1,size=.5, linetype = "dashed")+
  geom_ribbon(aes(x=Julian, ymin= exp(Q_2.5), 
                  ymax= (exp(Q_97.5)),fill=SciName ),alpha=.2 )+
  geom_line(aes(color=SciName),alpha=.85, size=1)+
  scale_fill_brewer(palette="Set1", name="Species")+
  scale_color_brewer(palette="Set1", name="Species")+

  ylab("log(Predicted mosquito count)") + 
  xlab("Julian date")+ theme_classic()+
  theme( legend.key.size = unit(1, "cm"),
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

