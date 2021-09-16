###### Powell Center: Phenological patterns of mosquitoes #######

# Travis McDevitt-Galles
# 08/24/2021
# title: 01_Pop_Model

# Recreating Mevin's modeling framework to see if i can retain

library(dplyr)
library(ggplot2)
library(patchwork)
library(rstan)
library(rstanarm)
library(matrixStats)
library(gamm4)

#Set working directory
setwd("C:/Users/tmcdevitt-galles/Documents/Population_dynamics")


### Loading in simple dataset

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

full.df$Julian <- full.df$DOY+(full.df$Year-min(full.df$Year))*365

## modifying the dataset 

stan.df <- full.df %>% filter(SubsetWeight >0 & TotalWeight >0 )

stan.df <- stan.df %>% group_by(SciName, Julian, DOY, Plot, Domain, Year) %>% 
  summarise(
    Count = sum(Count),
    TrapHours = sum(TrapHours),
    SubsetWeight = sum( SubsetWeight),
    TotalWeight = sum( TotalWeight ),
    Tmin7 = mean(Tmin7),
    PPT14 = mean(PPT14)
  )


stan.df$Offset <- (stan.df$TrapHours/24) * 
                   (stan.df$SubsetWeight/stan.df$TotalWeight)



stan.df %>% 
  ggplot(aes(x=DOY,y= Count*Offset  ,color=as.factor(Year)) )+ geom_point()

stan.df %>% 
  ggplot(aes(x=Julian,y= log10(Count*Offset+1), color=as.factor(Year) ))+
        geom_point()+     scale_color_brewer(palette="Set1", name="Year")  +
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



###### Getting data ready for stan model

Plot_Samples <- stan.df %>%  group_by(Plot) %>% summarize(Samples = n())
         

stan.df <- stan.df %>% arrange(Plot, Julian)


co.matrix <- model.matrix(stan.df$Count~ scale(stan.df$Tmin7) *
                            scale(stan.df$PPT14))

stan_d <- list( N = nrow(stan.df), P = ncol(co.matrix),
                X= co.matrix, G = nrow(Plot_Samples),
                group_samp = Plot_Samples$Samples,
                Y= stan.df$Count, offset = stan.df$Offset )

ar_output4 <- stan( 'Scripts/Stan/Ar_take2.stan', 
                      data=stan_d, iter = 5000,
                    control = list(max_treedepth = 20))


saveRDS(ar_output4 ,"bayes.rds")

## print estimated coefficients
readRDS("bayes.rds")

print(ar_output4, pars = c("rho", 'beta'))
traceplot(ar_output3)

post <- extract(ar_output4)

plot(ar_output4)


## plotting rho values

rho.df <- as.data.frame(post$rho)


colnames(rho.df) <- c((Plot_Samples$Plot))


rho.df <- tidyr::pivot_longer(rho.df, cols= starts_with("UNDE"),
                              names_to="Plot", values_to = "estimate")



ggplot(rho.df, aes(x=Plot, y=estimate))+ geom_violin(alpha=.5,fill="gray")+
  geom_boxplot(color="black", width=.05, outlier.size = 0) +theme_classic()+
  ylab("Rho estimate")+
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



## plotting rho values

beta.df <- as.data.frame(post$beta)


colnames(beta.df) <- c("Intercept", "TMin", "PPT", "TMin:PPT")


beta.df <- tidyr::pivot_longer(beta.df, cols= 1:4, names_to="Beta",
                               values_to = "Estimate")



ggplot(beta.df, aes(x=Beta, y=Estimate))+  geom_violin(alpha=.5,fill="gray")+
  geom_boxplot(color="black", width=.05, outlier.size = 0) +
  geom_hline(yintercept = 0, size=1 ) +theme_classic()+
  ylab("Beta estimate")+ xlab("Parameters")+
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





plot.df <- stan.df

  plot.df$Pred <- c(colMedians(as.matrix(post$m)))

test <- colQuantiles(as.matrix(post$m), probs=c(0.025,0.975))

plot.df$Q_2.5  <- test[,1]

plot.df$Q_97.5  <- test[,2]

plot.df %>% 
  ggplot( aes(x=Julian, y=log10(Count*Offset+1))) + 
  geom_ribbon(aes(x=Julian, ymin=Q_2.5*.52, ymax=(Q_97.5*.52)))+
  geom_point(color="red", size=.7)+ facet_wrap(~Plot, scales ="free")

plot.df %>% 
  ggplot( aes(x=Julian, y=log10((Count*Offset)+1))) + 
  geom_point(aes(x= Julian, y= Pred*.52),size=1)+
  geom_point(color="red", size=.7)+ facet_wrap(~Plot, scales ="free")


plot.df %>% 
  ggplot( aes(x=Julian, y=log10((Count*Offset)+1), color = as.factor(Year))) +
  geom_point(aes(x= Julian, y= Pred*Offset),size=2)+
  geom_point(color="black", size=1.5,alpha=.5)+ylab("Mosquito Count") + 
  xlab("Julian date")+  facet_wrap(~Plot)+
  scale_color_brewer(palette = "Set1",
                     name="Predicted count year")+ theme_classic()+
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


plot.df %>% 
  ggplot( aes(x=Julian, y=log10(Count*Offset+1))) + 
  geom_ribbon(aes(x=Julian, ymin=0, ymax=(Q_97.5*Offset)))+
  geom_point(color="red", size=1.5)+ facet_wrap(~Year, scales ="free")

plot.df %>% 
  ggplot( aes(x=exp(Q_97.5), y=exp(Pred))) +geom_point() +
  geom_abline(slope=1, intercept=0)





#### New attempt at doing some back casting i guess?

weather.df <- read.csv("Data/simple_weather.csv")



weather.df$Julian <- weather.df$DOY+(weather.df$Year-min(weather.df$Year))*365

Weather <- model.matrix(~ scale(weather.df$PPT)*scale(weather.df$TMIN))






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

Weather <- model.matrix(~scale(weather.df$STemp)*scale(weather.df$mPPT))



stan_d <- list( N = nrow(stan.df), P = ncol(Weather),
                Time = nrow(Weather),
                Julian = stan.df$Julian,
                X= Weather,
                Y= stan.df$Count, offset = stan.df$Offset )

ar_output5 <- stan( 'Scripts/Stan/Ar_take_3.stan', 
                    data=stan_d, iter = 2000,
                    control = list(max_treedepth = 10))


print(ar_output5, pars = c("rho", 'beta'))
traceplot(ar_output5)

post <- extract(ar_output5)

plot(x=1:nrow(weather.df),y =exp(colMedians(as.matrix(post$m))), type= "l")
points(x=stan.df$Julian, y= log10((stan.df$Count*stan.df$Offset+1)))

saveRDS(ar_output5 ,"bayes_back.rds")

ar_output5 <- readRDS("bayes_back.rds")

## plotting rho values

rho.df <- as.data.frame(post$rho)

colnames(rho.df) <- "Rho"
## plotting rho values

beta.df <- as.data.frame(post$beta)


colnames(beta.df) <- c("Intercept", "TMin", "PPT", "TMin:PPT")

beta.df <- cbind.data.frame(beta.df, rho.df)

beta.df <- tidyr::pivot_longer(beta.df, cols= 1:5, names_to="Parameters",
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



plot.df <- as.data.frame(colMedians(as.matrix(post$m)))

test <- colQuantiles(as.matrix(post$m), probs=c(0.025,0.975))

plot.df$Q_2.5  <- test[,1]

plot.df$Q_97.5  <- test[,2]
plot.df$Julian <- 1:nrow(plot.df)

colnames(plot.df)[1] <- "Count"

plot.df %>% 
  ggplot( aes(x=Julian, y=log10(exp(Count)+1))) +
  geom_ribbon(aes(x=Julian, ymin= log10(exp(Q_2.5)+1), 
                  ymax= log10(exp(Q_97.5)+1) ), fill ="grey",alpha=.7 )+
  geom_line(color="black",alpha=.75, size=2)+
  geom_point(data=stan.df, aes(x=Julian, y=log10(( Count /TrapHours) +1),
                               color=as.factor(Year)) )+
  scale_color_brewer(palette="Set1", name="Year")+
  ylab("log10(Predicted mosquito count)") + 
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


### 

test <- colQuantiles(as.matrix(post$m), probs=c(0.025,0.975))

plot.df$Q_2.5  <- test[,1]

plot.df$Q_97.5  <- test[,2]
plot.df$Julian <- 1:nrow(plot.df)

colnames(plot.df)[1] <- "Count"

at1 <- plot.df %>% filter(Julian >800 & Julian < 1100)

stan1.df <- stan.df %>% filter(Julian >800 & Julian < 1100)

clor <- RColorBrewer::brewer.pal(1,"Set1")

at1 %>% 
  ggplot( aes(x=Julian, y=log10(exp(Count)+1))) +
  geom_ribbon(aes(x=Julian, ymin= log10(exp(Q_2.5)+1), 
                  ymax= log10(exp(Q_97.5)+1) ), fill ="grey",alpha=.7 )+
  geom_line(color="black", alpha=.75, size=2) +
  geom_point(data=stan1.df, aes(x=Julian, y=log10(( Count /TrapHours) +1),
                              color=as.factor(Year)) )+
  scale_color_manual( values ="#377EB8", name="Year")+
  ylab("log10(Predicted mosquito count)") + ylim(0,2.5)+
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





plot.df %>% 
  ggplot( aes(x=Julian, y=log10((Count*Offset)+1))) + 
  geom_point(aes(x= Julian, y= Pred*.52),size=1)+
  geom_point(color="red", size=.7)+ facet_wrap(~Plot, scales ="free")


plot.df %>% 
  ggplot( aes(x=Julian, y=log10((Count*Offset)+1), color = as.factor(Year))) +
  geom_point(aes(x= Julian, y= Pred*Offset),size=2)+
  geom_point(color="black", size=1.5,alpha=.5)+ylab("Mosquito Count") + 
  xlab("Julian date")
  scale_color_brewer(palette = "Set1",
                     name="Predicted count year")+ theme_classic()+
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


plot.df %>% 
  ggplot( aes(x=Julian, y=log10(Count*Offset+1))) + 
  geom_ribbon(aes(x=Julian, ymin=0, ymax=(Q_97.5*Offset)))+
  geom_point(color="red", size=1.5)+ facet_wrap(~Year, scales ="free")

plot.df %>% 
  ggplot( aes(x=exp(Q_97.5), y=exp(Pred))) +geom_point() +
  geom_abline(slope=1, intercept=0)



### Removing a year


nRm.df <- stan.df %>% filter(Year != 2016)

stan_d <- list( N = nrow(nRm.df ), P = ncol(Weather),
                Time = nrow(Weather),
                Julian = nRm.df$Julian,
                X= Weather,
                Y= nRm.df$Count, offset = nRm.df$Offset )

ar_output5 <- stan( 'Scripts/Stan/Ar_take_3.stan', 
                    data=stan_d, iter = 2000,
                    control = list(max_treedepth = 10))


print(ar_output5, pars = c("rho", 'beta'))
traceplot(ar_output5)

post <- extract(ar_output5)

plot(x=1:nrow(weather.df),y =exp(colMedians(as.matrix(post$m))), type= "l")
points(x=stan.df$Julian, y= log10((stan.df$Count*stan.df$Offset+1)))

saveRDS(ar_output5 ,"bayes_back.rds")


## plotting rho values

rho.df <- as.data.frame(post$rho)

colnames(rho.df) <- "Rho"
## plotting rho values

beta.df <- as.data.frame(post$beta)


colnames(beta.df) <- c("Intercept", "TMin", "PPT", "TMin:PPT")

beta.df <- cbind.data.frame(beta.df, rho.df)

beta.df <- tidyr::pivot_longer(beta.df, cols= 1:5, names_to="Parameters",
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




plot.df <- as.data.frame(colMedians(as.matrix(post$m)))

test <- colQuantiles(as.matrix(post$m), probs=c(0.025,0.975))

plot.df$Q_2.5  <- test[,1]

plot.df$Q_97.5  <- test[,2]
plot.df$Julian <- 1:nrow(plot.df)

colnames(plot.df)[1] <- "Count"

stan2.df <- stan.df %>% filter(Year != 2016)

plot.df %>% 
  ggplot( aes(x=Julian, y=log10(exp(Count)+1))) +
  geom_ribbon(aes(x=Julian, ymin= log10(exp(Q_2.5)+1), 
                  ymax= log10(exp(Q_97.5)+1) ), fill ="grey",alpha=.7 )+
  geom_line(color="black",alpha=.75, size=2)+
  geom_point(data=stan2.df, aes(x=Julian, y=log10(( Count /TrapHours) +1),
                               color=as.factor(Year)) )+
  scale_color_brewer(palette="Set1", name="Year")+
  ylab("log10(Predicted mosquito count)") + 
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


### 

test <- colQuantiles(as.matrix(post$m), probs=c(0.025,0.975))

plot.df$Q_2.5  <- test[,1]

plot.df$Q_97.5  <- test[,2]
plot.df$Julian <- 1:nrow(plot.df)

colnames(plot.df)[1] <- "Count"

at1 <- plot.df %>% filter(Julian >800 & Julian < 1100)

stan1.df <- stan.df %>% filter(Julian >800 & Julian < 1100)

clor <- RColorBrewer::brewer.pal(1,"Set1")

at1 %>% 
  ggplot( aes(x=Julian, y=log10(exp(Count)+1))) +
  geom_ribbon(aes(x=Julian, ymin= log10(exp(Q_2.5)+1), 
                  ymax= log10(exp(Q_97.5)+1) ), fill ="grey",alpha=.7 )+
  geom_line(color="black", alpha=.75, size=2) +
  geom_point(data=stan1.df, aes(x=Julian, y=log10(( Count /TrapHours) +1),
                                color=as.factor(Year)) )+
  scale_color_manual( values ="#377EB8", name="Year")+
  ylab("log10(Predicted mosquito count)") +# ylim(0,2.5)+
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



