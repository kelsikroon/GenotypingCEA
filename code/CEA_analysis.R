# CEA strategies - November 2024
library(foreach)
library(doParallel)
library(dplyr)
library(haven)
library(lubridate)
library(reshape2)
library(xlsx)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggtext)

file.end <- format(Sys.Date(), "%Y%m%d")


#### Analysis preparation #### 
source('CEA_functions.R')
load("Bootstrap/bootstrap_non16"); params_non16 <- pars[! pars[,3] > 0.95,]
load("Bootstrap/bootstrap_16"); params_16 <- pars[! pars[,3] > 0.95,]
rm(pars)

set.seed(130601)

bmd.strategies <- c(rep("referAll", 5), rep("refer7", 4), rep("refer1618", 3), rep("repeat", 2))
normal.strategies <- c("referAll", "refer7", "refer1618", "repeat", "repeat7", "refer7", "refer1618", 
                       "repeat", "repeat7", "refer1618", "repeat", "repeat7","repeat", "repeat7" )
strategies <- data.frame(bmd.strategies, normal.strategies)

myCluster <- makeCluster(4) # type of cluster
registerDoParallel(myCluster)

onsetTimes <- seq(0, 10, 1)
n.reps <- 1000

#### 1. No bootstrapping #### 
resultsCEA_onerun <- foreach(i = 1:length(bmd.strategies)) %dopar% {
  load('~/Desktop/PhD/Projects/Colposcopy CEA/Analysis/ref.data.cens.RData') # data already organised/prepared for analysis
  source("CEA_functions.R") # functions for CEA and to get number of referrals
  source("CEA_parameters.R") # file for setting parameters
  library(dplyr)
  temp <- list() # loop through all combinations of BMD and normal strategies and calculate the NHB
  for (t in 1:length(onsetTimes)){ # loop through the 7 time points and store results for each strategy combination of normal & BMD at that time since onset
    temp[[t]] <- npBS_NHB(ref.data.cens, bmd.strategies[i], normal.strategies[i], n.reps, onsetTimes[t], trichot=F, sample.vink=F, bootstrap=F)
  }
  temp
}

## organise the data into a nice data frame and take difference from most aggressive strategy (refer all), and save into a excel file for manuscript
resultsCEA_onerun <- organiseAndSave(resultsCEA_onerun, "onerun")

# plot: 1100x600 size
plotCEA(resultsCEA_onerun, bootstrap = F, sample.vink = F) + geom_vline(xintercept=4, col='grey', lwd=0.5, lty=2)

#### 2. Fixed estimates from Vink #### 
resultsCEA_temp <- foreach(i = 1:length(bmd.strategies)) %dopar% {
  load('~/Desktop/PhD/Projects/Colposcopy CEA/Analysis/ref.data.cens.RData') # data already organised/prepared for analysis
  source("CEA_functions.R") # functions for CEA and to get number of referrals
  source("CEA_parameters.R") # file for setting parameters
  library(dplyr)
  temp <- list() # loop through all combinations of BMD and normal strategies and calculate the NHB
  for (t in 1:length(onsetTimes)){ # loop through the 7 time points and store results for each strategy combination of normal & BMD at that time since onset
    temp[[t]] <- npBS_NHB(ref.data.cens, bmd.strategies[i], normal.strategies[i], n.reps, onsetTimes[t], trichot=F, sample.vink=F, bootstrap=T)
  }
  temp
}

## organise the data into a nice data frame and take difference from most aggressive strategy (refer all), and save into a excel file for manuscript
resultsCEA <- organiseAndSave(resultsCEA_temp, "fixed") # load("resultsCEA_fixed_20241120.RData")

# plot: 1100x600 size
plotCEA(resultsCEA, bootstrap = T, sample.vink = F)


#### 3. Sampled parameters from Vink #### 
resultsCEA_bootstrap_temp <- foreach::foreach(i = 1:length(bmd.strategies)) %dopar% {
  load('~/Desktop/PhD/Projects/Colposcopy CEA/Analysis/ref.data.cens.RData') # data already organised/prepared for analysis
  source("CEA_functions.R") # functions for CEA and to get number of referrals
  source("CEA_parameters.R") # file for setting parameters
  library(dplyr)
  temp <- list() # loop through all combinations of BMD and normal strategies and calculate the NHB
  for (t in 1:length(onsetTimes)){ # loop through the 7 time points and store results for each strategy combination of normal & BMD at that time since onset
    temp[[t]] <- npBS_NHB(ref.data.cens, bmd.strategies[i], normal.strategies[i], n.reps, onsetTimes[t], trichot=F, sample.vink=T, bootstrap=T)
  }
  temp
}

## organise the data into a nice data frame and take difference from most aggressive strategy (refer all), and save into a excel file for manuscript
resultsCEA_bs <- organiseAndSave(resultsCEA_bootstrap_temp, "sampled")

# plot: 950x500 size
plotCEA(resultsCEA_bs, bootstrap = T, sample.vink = T)
