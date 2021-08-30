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
setwd("C:/Users/tmcdevitt-galles/Documents/Population_model")

#
