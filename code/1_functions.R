options(dplyr.summarise.inform = FALSE) # keep dplyr messages quiet 


#------------------------------------------
# Function to show k decimal places of x even if there are zeros
spec_dec <- function(x, k) trimws(format(round(x, k), nsmall=k)) 

#------------------------------------------
# Function to perform one rep of a bootstrap sample for all strategies on the same resampled data set - used for the probabilistic analysis 
one_rep <- function(dat, timesinceonset, sample.vink=F){
  if (sample.vink){ 
    shape16 <- sample(params_16[,1], 1, replace=T) ; scale16 <- sample(params_16[,2], 1, replace=T) ; weight16 <- sample(params_16[,3], 1, replace=T)
    shape <- sample(params_non16[,1], 1, replace=T) ; scale <- sample(params_non16[,2], 1, replace=T) ; weight <-  sample(params_non16[,3], 1, replace=T)
  }else{
    shape16 <- 3.33 ; scale16 <- 9.67 ; weight16 <- 0.44
    shape <- 9.14 ; scale <- 2.49 ; weight <- 0.31
  }
  Probcancer5years_1618 <-  pgamma(5+timesinceonset,shape=shape16,scale=scale16) - pgamma(timesinceonset,shape=shape16,scale=scale16) 
  Probcancer5years_1618 <- Probcancer5years_1618 / (1 - pgamma(timesinceonset,shape=shape16,scale=scale16)) * cumulativeprogcin3for1618_lifetime
  Probcancer5years_no1618 <- pgamma(5+timesinceonset,shape=shape,scale=scale) - pgamma(timesinceonset,shape=shape,scale=scale)
  Probcancer5years_no1618 <- Probcancer5years_no1618 / (1 - pgamma(timesinceonset,shape=shape,scale=scale) ) * cumulativeprogcin3for1618_lifetime *(weight/weight16) 
  
  probcancerParams <- c(Probcancer5years_1618, Probcancer5years_no1618)
  pob_params <- pobInput(dat) # calculate POBASCAM parameters from the resampled data set (proportion of regressive CIN)
  
  nmbs <- rep(NA, length(ascus.strategies)) # make empty vector to store bootstrap estimates 
  referred <- rep(NA, length(ascus.strategies))
  cxca_2rounds <- rep(NA, length(ascus.strategies))
  cin2_2rounds <- rep(NA, length(ascus.strategies))
  cin3_2rounds <- rep(NA, length(ascus.strategies))
  cin2_1rounds <- rep(NA, length(ascus.strategies))
  for (i in 1:length(ascus.strategies)){
    # calculate the costs and QALY's for the resampled data using only the intervention group 
    temp_res <- computeCostsandQALYloss(dat[dat$grp =='i',], hpv.strategies[i], ascus.strategies[i], normal.strategies[i], pob_params, probcancerParams) 
    
    #if(i %in% c(9, 11, 17)) print(c(i, temp_res[1:3], temp_res[6], cin23_2rounds = unname(temp_res[4]) + unname(temp_res[5])))
    
    nmbs[i] <- (-temp_res[2]*ICERthreshold - temp_res[1])/1057 # change to NMB
    # table for supplementary material:
    referred[i] <- temp_res[3]*100000/1057
    cxca_2rounds[i] <- temp_res[6]*100000/1057
    cin2_2rounds[i] <- temp_res[4]*100000/1057
    cin3_2rounds[i] <- temp_res[5]*100000/1057
    cin2_1rounds[i] <- temp_res[7]*100000/1057
  }
  NMB_d <- nmbs - nmbs[1]
  result <- data.frame(strategy = 1:19, hpv.strategies, ascus.strategies, normal.strategies, time = timesinceonset, nmbs, NMB_d, referred, cxca_2rounds, cin23_2rounds = cin2_2rounds + cin3_2rounds)
  return(result)
}


#------------------------------------------
# Function that computes costs and QALYS for a specific strategy (edited from Hans' code)
# --> Input parameters from "2_parameters.R" file 
# * Costrepeattest,CostnoCIN2,CostCIN2,CostCIN3,CostCancer,Costdeath, (Various sources)
# * LossCIN23,Lossdiagnosiscancer,Lossmonitor,Losspreterminal,Lossterminal (McDonald et al., 2017, CCC)
# * disce,discc (Dutch guidelines for economic)
# * Maxlifeyearslost,Notdeathfromothercauses5yr
# * Probcancersdetected_1618,Probcancersdetected_no1618 (national registry data)

# --> Input parameters that are obtained from POBASCAM data (recalculated for each resampled dataset)
# Numsubjects,ProbregressiveCIN2_1618,ProbregressiveCIN2_no1618,

# --> Input parameters that are obtained by running the most aggressive strategy (refer everyone)
# MaxnumberCIN2round1_1618,MaxnumberCIN2round1_no1618,MaxnumberCIN3round1_1618,MaxnumberCIN3round1_no1618,

# --> Input parameters that are obtained for the specific strategy under consideration by running the inputNumbers() function on the resampled POBASCAM data
# CIN2round1_1618,CIN2round1_no1618,CIN3round1_1618,CIN3round1_no1618,NoCIN2round1, 
# Firsttest,Repeattest,Probcancer5years,Survivalyr,Survivalyr2,

computeCostsandQALYloss <- function(dat, hpvStrategy, ascusStrategy, normalStrategy, pobParams, probCancerParams)   {
  ProbregressiveCIN2_no1618 <- pobParams$ProbregressiveCIN2_no1618
  ProbregressiveCIN2_1618 <- pobParams$ProbregressiveCIN2_1618
  
  Probcancer5years_1618 <- probCancerParams[1]
  Probcancer5years_no1618 <- probCancerParams[2]
  
  Numsubjects <- dim(dat)[1]
  referAll_res <- unname(inputNumbers_HPVfirst(dat, 'none', "referAll", "referAll")) # get input numbers for most aggressive strategy for the data (changes when resampled) (used for maximum number of cases)
  MaxnumberCIN2round1_1618 <- referAll_res[1] 
  MaxnumberCIN2round1_no1618 <- referAll_res[2] 
  MaxnumberCIN3round1_1618 <- referAll_res[3] 
  MaxnumberCIN3round1_no1618 <- referAll_res[4] 
  
  input_data <- unname(inputNumbers_HPVfirst(dat, hpvStrategy, ascusStrategy, normalStrategy)) # get input numbers for current strategy for the data  
  
  CIN2round1_1618 <- input_data[1] # assign values based on output from inputNumbers_HPVfirst function
  CIN2round1_no1618 <- input_data[2]
  CIN3round1_1618 <- input_data[3]
  CIN3round1_no1618 <- input_data[4]
  NoCIN2round1 <- input_data[5]
  Firsttest <- input_data[6]
  Repeattest <- input_data[7]
  
  MissednumberCIN2round1_1618 <-   MaxnumberCIN2round1_1618 - CIN2round1_1618
  MissednumberCIN3round1_1618 <- (1 - Probcancersdetected_1618)*(MaxnumberCIN3round1_1618 - CIN3round1_1618)
  MissednumberCIN2round1_no1618 <- MaxnumberCIN2round1_no1618 - CIN2round1_no1618
  MissednumberCIN3round1_no1618 <- (1 - Probcancersdetected_no1618)*(MaxnumberCIN3round1_no1618 - CIN3round1_no1618)
  Missedcancersround1 <- Probcancersdetected_1618*(MaxnumberCIN3round1_1618 - CIN3round1_1618) + Probcancersdetected_no1618*(MaxnumberCIN3round1_no1618 - CIN3round1_no1618)
  
  PersistentnumberCIN2round2_1618 <- MissednumberCIN2round1_1618*(1 - ProbregressiveCIN2_1618)*(1 - Probcancer5years_1618) 
  PersistentnumberCIN3round2_1618 <- MissednumberCIN3round1_1618*(1 - Probcancer5years_1618) 
  PersistentnumberCIN2round2_no1618 <- MissednumberCIN2round1_no1618*(1 - ProbregressiveCIN2_no1618)*(1 - Probcancer5years_no1618) 
  PersistentnumberCIN3round2_no1618 <- MissednumberCIN3round1_no1618*(1 - Probcancer5years_no1618) 
  
  CIN2round2 <- PersistentnumberCIN2round2_1618 + PersistentnumberCIN2round2_no1618
  CIN3round2 <- PersistentnumberCIN3round2_1618 + PersistentnumberCIN3round2_no1618
  Cancerround2new  <- MissednumberCIN2round1_1618*(1 - ProbregressiveCIN2_1618)*Probcancer5years_1618 +
    MissednumberCIN3round1_1618*Probcancer5years_1618 +
    MissednumberCIN2round1_no1618*(1 - ProbregressiveCIN2_no1618)*Probcancer5years_no1618 +
    MissednumberCIN3round1_no1618*Probcancer5years_no1618
  Cancerround2old <- Missedcancersround1
  
  Probsurvival <- Survivalyr[10]
  Probsurvival2 <- Survivalyr2[10]
  
  CIN2round1 <- CIN2round1_1618 + CIN2round1_no1618
  CIN3round1 <- (1 - Probcancersdetected_1618)*CIN3round1_1618 + (1 - Probcancersdetected_no1618)*CIN3round1_no1618 
  Cancerround1 <- Probcancersdetected_1618*CIN3round1_1618 + Probcancersdetected_no1618*CIN3round1_no1618
  
  Survivalyr <- c(1, Survivalyr)
  Survivalyr2 <- c(1, Survivalyr2)

  Cost <- NoCIN2round1*CostnoCIN2 +
    CIN2round1*CostCIN2 +
    CIN3round1*CostCIN3 +
    Firsttest*Costfirsttest +
    (Repeattest*Costrepeattest)/discc^1 +   # +1 year discounting
    (Cancerround1*Costcancer)/discc^0.25 +  # +0.25 year discounting
    (CIN2round2*CostCIN2 + CIN3round2*CostCIN3)/discc^5 +
    (Notdeathfromothercauses5yr*Costcancer*(Cancerround2old))/discc^(5 + 0.25) +  # +0.25 year discounting
    (Notdeathfromothercauses5yr*Costcancer*(Cancerround2new))/discc^(5 + 0.25)    # +0.25 year discounting
  
  for (i in 0:9){ # add costs for cancer non-survivors:
    Cost <- Cost + 
      Cancerround1*(Survivalyr[i + 1]-Survivalyr[i + 2])*Costdeath/discc^(i + .5) +
      Notdeathfromothercauses5yr*Cancerround2old*(Survivalyr2[i + 1]-Survivalyr2[i + 2])*Costdeath/discc^(t_reset + i +.5) + 
      Notdeathfromothercauses5yr*Cancerround2new*(Survivalyr[i + 1]-Survivalyr[i + 2])*Costdeath/discc^(5 + i + .5)
  }
  
  QALYloss <- NoCIN2round1*0.5*LossCIN01 + # Note: LOSSCIN01 is always set to zero except for one scenario in the sensitivity analysis
    ((CIN2round1 + CIN3round1)*0.5*LossCIN23/disce^0.25) + #  0.25 year discounting (because it lasts 0.5 year).
    ((CIN2round2 + CIN3round2)*0.5*LossCIN23/disce^(5 + 0.25))
  
  # cancer round1 ----
  QALYloss <- QALYloss +
    Probsurvival*Cancerround1*0.5*(Lossdiagnosiscancer) + # cancer survivors
    Probsurvival*Cancerround1*3.5*(Lossmonitor)/disce^2.25 + # cancer survivors
    sum(sapply(1:9, function(i) (Survivalyr[i+1]-Survivalyr[i+2])*Cancerround1*(0.5)*Lossdiagnosiscancer/disce^(0.25))) +    # cancer non-survivors: loss diagnosis, does not apply for 1st year deaths
    sum(sapply(1:9, function(i) (Survivalyr[i+1]-Survivalyr[i+2])*Cancerround1*(i-0.5)*Lossmonitor/disce^((i+0.5)/2))) +     # cancer non-survivors: loss remission, does not apply for 1st year deaths
    sum(sapply(0:9, function(i) (Survivalyr[i+1]-Survivalyr[i+2])*Cancerround1*(0.5)*Losspreterminal/disce^(i+0.25) )) +     # cancer non-survivors: loss (pre)terminal, applies for all years
    sum(sapply(0:9, function(i) (Survivalyr[i+1]-Survivalyr[i+2])*Cancerround1*((1 - (1/disce)^(Maxlifeyearslost-i))/(1-(1/disce)))/disce^(i+0.5))) # cancer non-survivors: life years lost (dead), applies for all years
  
  
  # cancer round2 old ----
  # t_reset = 3 (because some cases are detected by symptoms. Calculation in 2_parameters.R file)
  QALYloss <- QALYloss +
    Notdeathfromothercauses5yr*Probsurvival2*Cancerround2old*0.5*(Lossdiagnosiscancer)/disce^t_reset +
    Notdeathfromothercauses5yr*Probsurvival2*Cancerround2old*3.5*(Lossmonitor)/disce^(t_reset + 2.25) +
    sum(sapply(1:9, function(i) Notdeathfromothercauses5yr*(Survivalyr2[i+1]-Survivalyr2[i+2])*Cancerround2old*(0.5)*Lossdiagnosiscancer/disce^(0.25))) +
    sum(sapply(1:9, function(i) Notdeathfromothercauses5yr*(Survivalyr2[i+1]-Survivalyr2[i+2])*Cancerround2old*(i-0.5)*Lossmonitor/disce^(t_reset + (i+0.5)/2))) +
    sum(sapply(0:9, function(i) Notdeathfromothercauses5yr*(Survivalyr2[i+1]-Survivalyr2[i+2])*Cancerround2old*(0.5)*Losspreterminal/disce^(i+t_reset+0.25))) +
    sum(sapply(0:9, function(i) Notdeathfromothercauses5yr*(Survivalyr2[i+1]-Survivalyr2[i+2])*Cancerround2old*((1 - (1/disce)^(Maxlifeyearslost-t_reset-i))/(1-(1/disce)))/disce^(i+t_reset +.5)))
  
  # cancer round2 new ----
  QALYloss <- QALYloss +
    Notdeathfromothercauses5yr*Probsurvival*Cancerround2new*0.5*(Lossdiagnosiscancer)/disce^5 +
    Notdeathfromothercauses5yr*Probsurvival*Cancerround2new*3.5*(Lossmonitor)/disce^7.25 +
    sum(sapply(1:9, function(i) Notdeathfromothercauses5yr*(Survivalyr[i+1]-Survivalyr[i+2])*Cancerround2new*(0.5)*Lossdiagnosiscancer/disce^(0.25))) +
    sum(sapply(1:9, function(i) Notdeathfromothercauses5yr*(Survivalyr[i+1]-Survivalyr[i+2])*Cancerround2new*(i-0.5)*Lossmonitor/disce^(5 + (i+0.5)/2))) +
    sum(sapply(0:9, function(i) Notdeathfromothercauses5yr*(Survivalyr[i+1]-Survivalyr[i+2])*Cancerround2new*(0.5)*Losspreterminal/disce^(i+5.25))) +
    sum(sapply(0:9, function(i) Notdeathfromothercauses5yr*(Survivalyr[i+1]-Survivalyr[i+2])*Cancerround2new*((1 - (1/disce)^(Maxlifeyearslost-5-i))/(1-(1/disce)))/disce^(i+5.5)))
  
  return(c(cost=Cost, 
           qaly=QALYloss, 
           referred = input_data[8], 
           cin2_2rounds = CIN2round1 + CIN2round2, 
           cin3_2rounds = CIN3round1 + CIN3round2, 
           cxca_2rounds = Cancerround1 + Cancerround2old + Cancerround2new,
           CIN2round1)) 
}

#------------------------------------------
# Function to calculate crude number of referrals for a specific strategy based on cytology & genotype result 
# possible cytology strategies are: referAll, refer1618, refer7, repeat, repeat7, nextround
crudeCases <- function(data, strategy){
  if (strategy == "referAll"){
    baseline.referred <- data$id # everyone is immediately referred
    repeat.tested <- c() # so no one is repeat tested 
    repeat.referred <- c() # no one is referred after a repeat test 
    
  }else if (strategy =='refer1618'){
    baseline.referred <- data$id[which(data$hpv.1618pos)] # 16/18-pos are immediately referred 
    repeat.tested <- data$id[which(!data$hpv.1618pos)] # non-16/18 are repeated 
    repeat.referred <- data$id[which(!data$hpv.1618pos & data$first.repeat.cyt =='>=BMD')] # 16/18-negs >=BMD are referred after repeat test 
    
  }else if (strategy =='repeat'){
    baseline.referred <- c() # no one is immediately referred 
    repeat.tested <- data$id # everyone is repeated
    repeat.referred <- data$id[which(data$first.repeat.cyt =='>=BMD')] # >=BMD referred after repeat test 
    
  }else if (strategy == "repeat7"){
    baseline.referred <- c() # all positive for 7 types are repeated, otherwise go to next round 
    repeat.tested <- data$id[which(data$hpv.genotype.pos)] # only 7 types positive repeated
    repeat.referred <- data$id[which(data$hpv.genotype.pos & data$first.repeat.cyt =='>=BMD')] # 7-types-pos >=BMD are referred after repeat test 
    
  }else if (strategy == 'nextround'){
    baseline.referred <- c() # everyone goes to next round: no-one is immediately referred 
    repeat.referred <- c() # everyone goes to next round: no-one is repeated
    repeat.tested <- c() # everyone goes to next round: no-one is referred after the repeat test (no-one is repeated) 
    
  }else if (strategy == 'refer7'){
    baseline.referred <- data$id[which(data$hpv.genotype.pos)] # immediately refer 7-types positive 
    repeat.tested <- data$id[which(!data$hpv.genotype.pos)] # 7-types negative are repeated
    repeat.referred <- data$id[which(!data$hpv.genotype.pos & data$first.repeat.cyt =='>=BMD')] # 7-types-neg >=BMD are referred after repeat test 
  }
  first.round.referred <- c(baseline.referred, repeat.referred)
  results <- list(baseline.referred, repeat.referred, first.round.referred, repeat.tested)
  names(results) <- c("baseline.referred", "repeat.referred", "first.round.referred", "repeat.tested")
  
  return(results)
}

#------------------------------------------
# Function that combines results (i.e. number referred, number of CIN2/3+ cases) from seperate BMD and normal strategies for ***HPV first strategies***
# possible HPV stratgies are: "partial", "extended"
# possible cytology strategies are: referAll, refer1618, refer7, repeat, repeat7, nextround
inputNumbers_HPVfirst <-  function(data, hpvStrategy, ascusStrategy, normalStrategy){
  Numsubjects <- total <- dim(data)[1]
  
  hpv.column <- case_when(
    hpvStrategy == "partial" ~ "hpv.1618pos",
    hpvStrategy == "extended" ~ "hpv.genotype.pos",
    T ~ "none")
  
  if (hpv.column != 'none'){
    # first we refer based on genotype group: (check if genotype is positive)
    genotype.referred <- data$id[data[[hpv.column]] == 1] # if genotype is positive then these are counted as "genotype-referred"
    not.genotype.referred <- data[data[[hpv.column]] != 1,] # if genotype is negative then these patients are saved in new data called not.genotype.referred 
    
  }else{
    # if hpv.column == 'none' then none are referred based on HPV and all get cytology 
    genotype.referred <- c() # no one is genotype referred
    not.genotype.referred <- data # whole data set is "not-genotype-referred"
  }
  
  # using only "not.genotype.referred" we:
  # calculate numbers for those not referred based on HPV test seperate for each cytology group
  hsilRes <- crudeCases(not.genotype.referred[not.genotype.referred$baseline.cyt==">BMD", ], "referAll") # all remainingHSIL cases are all immediately referred at baseline 
  ascusRes <- crudeCases(not.genotype.referred[not.genotype.referred$baseline.cyt=="BMD", ], ascusStrategy) # count crude referrals for specific BMD strategy
  normalRes <- crudeCases(not.genotype.referred[not.genotype.referred$baseline.cyt=="normal", ], normalStrategy) # count crude referrals for specific normal strategy
  
  # add together all who were referred at baseline 
  baseline.referred <- length(c(genotype.referred, ascusRes$baseline.referred, normalRes$baseline.referred, hsilRes$baseline.referred)) # get total number of baseline referred
  referred <- c(genotype.referred, ascusRes$first.round.referred, normalRes$first.round.referred, hsilRes$first.round.referred) # combine Id's of those who were referred 
  
  CIN3round1_1618 <- sum(data[data$id %in% referred & data$genotype.data =='HPV 16/18', ]$first.round.cin3plus) # calculate CIN3+ 16/18 in first round  
  CIN2round1_1618 <- sum(data[data$id %in% referred & data$genotype.data =='HPV 16/18', ]$first.round.cin2plus) - CIN3round1_1618 # calculate CIN2 16/18 in first round by subtracting CIN3+ 16/18 from CIN2+ 16/18 
  
  CIN3round1_no1618 <- sum(data[data$id %in% referred & !data$genotype.data =='HPV 16/18', ]$first.round.cin3plus) # repeat for non-16/18
  CIN2round1_no1618 <- sum(data[data$id %in% referred & !data$genotype.data =='HPV 16/18', ]$first.round.cin2plus) - CIN3round1_no1618
  
  NoCIN2round1 <- length(referred) - CIN2round1_1618 - CIN3round1_1618 - CIN2round1_no1618 - CIN3round1_no1618 # find the number that were referred that had no CIN2 in round 1
  
  # number who had baseline cytology test --> those who were negative for initial HPV genotyping, or if strategy was referAll BMD & referAll Normal then it is zero
  Firsttest <- case_when(
    ascusStrategy =='referAll' & normalStrategy =='referAll' ~ 0,  # --> 0 if all HPV-positive were referred then no baseline cytology test is done on anyone 
    T ~ dim(not.genotype.referred)[1]  # --> otherwise everyone who wasn't positive for the genotype gets baseline cytology 
  )
  
  # repeat test: count number who had a repeat test (only normal or ascus who were not refered based on HPV genotype)
  Repeattest <- length(c(hsilRes$repeat.tested, ascusRes$repeat.tested, normalRes$repeat.tested))
  numReferred <- length(referred)
  return(c(CIN2round1_1618, CIN2round1_no1618, CIN3round1_1618, CIN3round1_no1618, NoCIN2round1, Firsttest, Repeattest, numReferred))
}


#------------------------------------------
# Function that creates the input data from POBASCAM required in the cost/QALY function, this function will be called in the bootstrap function 
pobInput <- function(dat){ # **needs full data not just intervention group**
  ### DATA POBASCAM BASELINE --> will use non-parametric bootstrap to recalculate this section
  
  ### Regressive CIN2: numbers in women with normal cytology  
  ### (rows: intervention , control // columns: number of subjects, number of cases over two rounds)
  RegressiveCIN2 <- dat[dat$baseline.cyt =='normal',] %>% group_by(grp) %>% summarise(n.sub = n(), n.case = sum(cin2plus.cases) - sum(cin3plus.cases)) %>% select(!"grp")
  
  ProbregressiveCIN2_1618 <- unlist(-1*diff(as.matrix(RegressiveCIN2[,2]/RegressiveCIN2[,1]))/(RegressiveCIN2[1,2]/RegressiveCIN2[1,1]))
  ProbregressiveCIN2_no1618 <- unlist(ProbregressiveCIN2_1618)
  
  return(list(ProbregressiveCIN2_1618 = ProbregressiveCIN2_1618, ProbregressiveCIN2_no1618 = ProbregressiveCIN2_no1618))
}



#------------------------------------------
# function to find the intersection between two lines using linear interpolation
find_intersection <- function(x1, y1, x2, y2) {
  intersections <- data.frame(x = numeric(0), y = numeric(0))
  
  # Create a common x-grid for interpolation
  x_common <- sort(unique(c(x1, x2)))  
  
  # Interpolate y-values on the common x-grid
  y1_interp <- approx(x1, y1, xout = x_common, rule = 2)$y
  y2_interp <- approx(x2, y2, xout = x_common, rule = 2)$y
  
  # Find intersections
  for (i in 1:(length(x_common) - 1)) {
    if ((y1_interp[i] - y2_interp[i]) * (y1_interp[i + 1] - y2_interp[i + 1]) < 0) {
      # Linear interpolation to find exact x where y1 = y2
      x_intersection <- x_common[i] - (y1_interp[i] - y2_interp[i]) * (x_common[i + 1] - x_common[i]) / (y1_interp[i + 1] - y2_interp[i + 1] - (y1_interp[i] - y2_interp[i]))
      
      y_intersection <- approx(x_common, y1_interp, xout = x_intersection)$y  # Interpolated y
      
      intersections <- rbind(intersections, data.frame(x = x_intersection, y = y_intersection))
    }
  }
  
  return(intersections)
}

#------------------------------------------
# Function that plots the NMB for each time point and labels the best strategy at each point in a different colour
plotCEA <- function(results, cuttime=14, label='long', printplot=T, maximumY = 16.5, labelLines = T, sensPlot=F){ # plot: 1100x600 size (or 2200 x 1200) (or 1550 x 850)
  results$time_orig <- results$time
  
  results$hpv <- results$hpv.strategies
  results$strategy2 <- results$strategy #rep(1:length(hpv.strategies), each=length(onsetTimes))
  results$hpv[results$hpv == 'extended'] <- '7-types'
  results$hpv[results$hpv == 'partial'] <- '16/18'
  results$hpv <- factor(results$hpv, levels = c("16/18", "7-types",  "cyt&genotype", "cytOnly", "allHPV+"))
  
  # find which strategy is best at each time point 
  #best.at.time <- results %>% group_by(time) %>% filter(NMB_d == max(NMB_d)) %>% ungroup() %>% arrange(time) %>% pull(strategy2)
  #results$best.at.time <- rep(best.at.time, length(hpv.strategies))
  best.at.time <-  results %>% group_by(time_orig) %>% filter(NMB_d == max(NMB_d)) %>% ungroup() %>% arrange(time) %>% pull(strategy2)
  
  results$strategy_long <-rep(c("**Best strategy: 1**<br>Refer all",
                                "**Best strategy: 2**<br>Refer 7-types-pos or HSIL/ASCUS 7-types-neg<br>Repeat all NILM", 
                                "**Best strategy: 3**<br>Refer 7-types-pos or HSIL/ASCUS 7-types-neg<br>Repeat only NILM 7-types-pos",
                                "**Best strategy: 4**<br>Refer 7-types-pos or HSIL 7-types-neg<br>Repeat ASCUS & NILM",
                                "**Best strategy: 5**<br>Refer 7-types-pos or HSIL 7-types-neg<br>Repeat ASCUS & only NILM 7-types-pos",
                                "**Best strategy: 6**<br>Refer all 16/18-pos or HSIL/ASCUS 16/18-neg<br>Repeat all NILM",
                                "**Best strategy: 7**<br>Refer all 16/18-pos or HSIL/ASCUS 16/18-neg<br>Repeat only NILM 7-types-pos",
                                "**Best strategy: 8**<br>Refer all 16/18-pos or HSIL 16/18-neg or ASCUS 7-types-pos<br>Repeat ASCUS 7-types-neg & NILM 16/18-neg",
                                "**Best strategy: 9**<br>Refer all 16/18-pos or HSIL 16/18-neg or ASCUS 7-types-pos<br>Repeat ASCUS 7-types-neg & NILM 7-types-pos",
                                "**Best strategy: 10**<br>Refer all 16/18-pos or HSIL 16/18-neg<br>Repeat ASCUS 16/18-neg & NILM 16/18-neg",
                                "**Best strategy: 11**<br>Refer all 16/18-pos or HSIL 16/18-neg<br>Repeat ASCUS 16/18-neg & NILM 7-types-pos",
                                "**Best strategy: 12**<br>Refer HSIL/ASCUS<br>Repeat all NILM",
                                "**Best strategy: 13**<br>Refer HSIL/ASCUS<br>Repeat only NILM 7-types-pos", 
                                "**Best strategy: 14**<br>Refer HSIL or ASCUS 7-types-pos<br>Repeat ASCUS 7-types-neg & all NILM",
                                "**Best strategy: 15**<br>Refer HSIL or ASCUS 7-types-pos<br>Repeat ASCUS 7-types-neg & NILM 7-types-pos",
                                "**Best strategy: 16**<br>Refer HSIL or ASCUS 16/18-pos<br>Repeat ASCUS 16/18-neg & all NILM", 
                                "**Best strategy: 17**<br>Refer HSIL or ASCUS 16/18-pos<br>Repeat ASCUS 16/18-neg & NILM 7-types-pos",
                                "**Best strategy: 18**<br>Refer HSIL<br>Repeat all ASCUS & NILM",
                                "**Best strategy: 19**<br>Refer HSIL<br>Repeat ASCUS & NILM 7-types-pos"), each = length(onsetTimes))
  
  if (any(results$NMB_d < 0)){
    # Detect where y crosses below 0 and compute interpolated x for y = 0
    y_fix <- results %>% group_by(strategy2) %>%
      mutate(y = NMB_d, y_next = lead(NMB_d), x=time, x_next = lead(x)) %>%  # Look ahead to next y-value
      dplyr::filter(!is.na(y_next)) %>%
      rowwise() %>%
      mutate(
        # Check if y is positive and y_next is negative -> interpolate x for y = 0
        x_zero = ifelse(y > 0 & y_next < 0, time - (y / (y_next - y)), NA),
        y_zero = ifelse(!is.na(x_zero), 0, NA)
      ) %>% ungroup() %>%
      filter(!is.na(x_zero)) %>% mutate(time = x_zero, NMB_d = y_zero)
    
    results <- bind_rows(results[results$NMB_d >= 0,], y_fix[,1:15]) %>% distinct() %>% arrange(strategy2,time)
    
  }
  
  strategy.counts <- table(best.at.time[1:length(onsetTimes)])
  best.at.time[best.at.time %in% as.numeric(names(strategy.counts)[strategy.counts==1])] <- NA # remove if they were only best for 1 year
  best.strategies <- best.at.time %>% unique() %>% na.omit() %>% as.numeric()
  
  results$col <- 'grey'
  extraColours <- rcartocolor::carto_pal(n=6, "Bold")
  for (i in 1:length(best.strategies)){
    results$col[results$strategy2 == best.strategies[i]] <- extraColours[i]
  }
  
  results$inBest <- ifelse(results$strategy2 %in% best.strategies,T, F)
  results$y <- results$NMB_d
  
  results <- subset(results, strategy2!=1) # remove reference strategy from plot as this is just a line at y=0
  
  # these were customised to previous final plot from last analysis to make nice labels 
  results$strategy3 <- paste0(results$strategy2, ", ", results$strategy2 - 1)
  # results$strategy3[results$strategy2 %in% c(1, seq(2, 20, 2))] <- ""
  # results$strategy3[results$inBest] <- results$strategy2[results$inBest]
  # results$strategy3[results$strategy2==16] <- "16, 13, 12"
  # results$strategy3[results$strategy2==8] <- "8"
  # results$strategy3[results$strategy3 == '13, 12'] <-''
  
  results$strategy3 <- case_when( # labels made specifically for these results
    results$strategy2 == 17 ~ "17",
    results$strategy2 == 15 ~ "11, 9, 15",
    results$strategy2 == 10 ~ "10, 8, 16, 14, 13, 7",
    results$strategy2 == 19 ~ "19",
    #results$strategy2 == 13 ~ "13, 7, 19",
    
    results$strategy2 == 18 ~ "18",
    results$strategy2 == 6 ~ "6, 12",
    
    results$strategy2 == 2 ~ "2",
    results$strategy2 == 3 ~ "3",
    results$strategy2 == 4 ~ "4",
    results$strategy2 == 5 ~ "5",
    T ~ ""
  )
  
  
  cea_plot <- ggplot(results, aes(x=time, y=y, group=strategy2, col=col, label=strategy2, lty=hpv)) + theme_classic() + 
    geom_line(data=results[! results$inBest & results$time < cuttime + 0.05,], lwd=0.4) + 
    geom_line(data=results[results$inBest & results$time < cuttime + 0.05,], lwd=.8) +
    xlab("Time since CIN2/3 onset (years)") + 
    ylab("NMB (\u20AC/woman) (relative to reference strategy: refer all)") + 
    theme(legend.position='inside', legend.position.inside = c(0.30, 0.13), axis.text = element_text(size=15), 
          text = element_text(size=13), plot.margin = unit(c(1,3,1,0.5), "cm")) + 
    scale_x_continuous(expand = c(0,0), limits=c(0, cuttime), breaks=seq(0, cuttime, 1))  + 
    scale_y_continuous(expand = c(0, 0), limits=c(-0.0001, maximumY +5), breaks=seq(0, maximumY, 100)) +
    scale_color_manual(values = setNames(results$col, results$col))   +
    scale_linetype_discrete(name="Triage strategy group", labels = c('Immediate 16/18 referral', "Immediate 7-types referral", "Cytology & genotype triage")) +
    guides(color='none') 
  
  if (labelLines){
    cea_plot <- cea_plot  +
      geom_text_repel(data=results[results$time_orig == cuttime & results$strategy3!="", ], aes(color = col, label = strategy3),
                      direction = "y", xlim = c(cuttime+0.05, Inf), size=3, hjust = 0, segment.size = 0.25, segment.alpha =.25, max.overlaps = 100,
                      box.padding = 0.05, show.legend = F, segment.curvature = 0.1, segment.ncp = 0, segment.angle = 10)
  }
  
  if (sensPlot){
    cea_plot <- cea_plot + 
      theme(legend.position='none', axis.text = element_text(size=12),  plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
            text = element_text(size=10))
  }
  
  # calculate where "best" lines cross 
  num.best <- length(best.strategies)
  last_change <- 0 
  changes <- rep(NA, num.best)
  for (i in 1:(num.best)){
    if (i == num.best){ # if this is the last best one 
      change_index <- cuttime 
      if (cuttime - last_change < 2) change_index <- cuttime + 5
      xline <- F
    } else{
      xline <- T
      change_index <- unlist(find_intersection(results$time[results$strategy2 == best.strategies[i]], #x1
                                               results$y[results$strategy2 == best.strategies[i]], #y1
                                               results$time[results$strategy2 == best.strategies[i+1]], #x2 
                                               results$y[results$strategy2 == best.strategies[i+1]])[1]) #y2
      change_index <- change_index[length(change_index)]
    }
    changes[i] <- change_index
    if (is.na(change_index)) next  
    
    temp_data <- data.frame(x = last_change + ((change_index - last_change - 1) / 2), 
                            y = maximumY-30, 
                            label = results$strategy_long[results$strategy2 == best.strategies[i]][1], # Find the strategy name for the label
                            strategy2 = best.strategies[i], 
                            hpv = results$hpv[results$strategy2 == best.strategies[i]][1]) 
    
    if (xline) cea_plot <- cea_plot + geom_segment(x = change_index, xend = change_index, y = 0, yend = maximumY, linetype = "dotted", linewidth = 0.25, color = 'black')
    
    if (label=='long'){ # label
      #cea_plot <- cea_plot + geom_richtext(data = temp_data, aes(x = x, y = y, label = paste0("**Best: ", strategy2, "**")),color = 'black', fill = 'white', label.color = NA, size = 5)
      
    }else{
      cea_plot <- cea_plot + geom_richtext(data = temp_data, aes(x = x + 0.5, y = y, label = paste0("**Best: ", strategy2, "**")),color = 'black', fill = 'white', label.color = NA, size = 3) +
        ylab("NMB (\u20AC/woman) (relative to Strategy 1)") 
    }
    last_change <- change_index  
  }
  cea_plot <- cea_plot + coord_cartesian(clip = "off") 
  print(changes)
  if(printplot) suppressWarnings(print(cea_plot))
  return(cea_plot)
}



#------------------------------------------
# Function that calculates the marginal 10-year if cancer was detected 5 years later (including symptoms and stage shift)
survival5yearsLater <- function(time1to23 = 8){
  roman_to_numeric <- function(roman) {
    roman_dict <- c("I" = 1, "II" = 2, "III" = 3, "IV" = 4)  # Extend if needed
    return(as.numeric(roman_dict[roman]))
  }
  
  # nki.dat <- read.xlsx(file='~/Desktop/PhD/Projects/Colposcopy CEA/NCR-export-2025_02_25_12_54_07.xlsx', startRow = 10, sheetIndex = 1 ) # Load data (skip first 10 rows because it is not data)
  # colnames(nki.dat) <- c("Years", "Stage", "Surv", "Ncases")
  # nki.dat <- nki.dat[nki.dat$Stage != 'Unknown', 1:4]
  # nki.dat$Years <- as.numeric(nki.dat$Years)
  # save(nki.dat, file='nki.dat.RData')
  load("data/NKI_data.RData")
  
  Fconv_sympt <- function(x, rate1, rate2, sympt1, sympt2){
    return((rate1*rate2 / (rate2 + sympt2 - rate1 - sympt1))*((1 - exp(-x*(rate1 + sympt1)))/(rate1 + sympt1) - (1 - exp(-x*(rate2 + sympt2)))/(rate2 + sympt2)))
  }
  
  sympt.dat <- nki.dat %>% mutate(Stage = roman_to_numeric(Stage)) %>% # convert from roman numerals to numbers, more natural to work with 
    mutate(Group = ifelse(Stage %in% c("2", "3"), "2/3", Stage)) %>% # add a new variable group which recodes 2 or 3 as 2/3 combined 
    group_by(Years, Group) %>%
    summarise(
      Surv = sum(Surv * Ncases) / sum(Ncases), # Weighted survival of stage 2 and 3 when combining them
      Ncases = sum(Ncases), .groups = "drop") %>% # sum number of cases in stage 2 and 3
    ungroup() %>% rename(Stage = Group) %>% as.data.frame() # rename the group variable as Stage and then only keep this variable 
  
  Ncases_5years <- c(33, NA) # empty vector for number of cases, one entry for each stage (1, 2/3, or 4) and each time since diagnosis (0 to 10)
  
  rate1 <- 1/time1to23 # rate1 is 1/8 (using the model from Tiago) (rate progress from stage 1 to 2/3)
  rate2 <- 1/4  # this stays fixed at 1/4 (informed from the data (rate progress from stage 2/3 to 4)
  sympt1 <- -log(1-4*0.15/time1to23) #Y years is (4/Y)x15% per year.
  sympt2 <- -log(1-0.4)
  
  prob_s1to1 <- (1 - (rate1/(rate1 + sympt1))*(1 - exp(-5 * (rate1 + sympt1)))) # 1 - P(progress to stage 2/3 and dont get symptoms)  
  prob_s1to2 <- (rate1/(rate1 + sympt1))*(1 - exp(-5 * (rate1 + sympt1))) # P(progress to stage 2/3 and dont get symptoms) 
  prob_s1to4 <- Fconv_sympt(5, rate1, rate2, sympt1, sympt2) # convolution of 2 exponentials 
  prob_s2to2 <- (1 - (rate2/(rate2 + sympt2))*(1 - exp(-5 * (rate2 + sympt2)))) # 1 - P(progress to stage 4 and dont get symptoms) 
  prob_s2to4 <- (rate2/(rate2 + sympt2))*(1 - exp(-5 * (rate2 + sympt2))) # P(progress to stage 4 and dont get symptoms) 
  prob_s4to4 <- 1
  
  i <- 1
  for (year in seq(0, 10)){ # loop through all 10 years 
    temp <- sympt.dat[sympt.dat$Years == year,] # keep only three values of Ncases for current time since diagnosis 
    
    newS1 <- temp$Ncases[1] * prob_s1to1  
    newS2 <- temp$Ncases[1] * prob_s1to2 + temp$Ncases[2] * prob_s2to2
    newS4 <- temp$Ncases[1] * prob_s1to4 + temp$Ncases[2] * prob_s2to4
    
    Ncases_5years[i:(i+2)] <- round(c(newS1, newS2, newS4)) # add 3 new stage values to the Ncases 5 years later variable 
    i <- i+3 # move index of Ncases_5years forward 
  }
  
  sympt.dat[["Ncases_5years"]] <- Ncases_5years
  
  Survivalyr <- sympt.dat %>% group_by(Years) %>% summarise(Survivalyr = weighted.mean(Surv, Ncases)/100) %>% select(Survivalyr) %>% unlist() %>% as.numeric()
  Survivalyr2 <- sympt.dat %>% group_by(Years) %>% summarise(Survivalyr2 = weighted.mean(Surv, Ncases_5years)/100) %>% select(Survivalyr2) %>% unlist() %>% as.numeric()
  
  return(data.frame(Survivalyr, Survivalyr2))
}


