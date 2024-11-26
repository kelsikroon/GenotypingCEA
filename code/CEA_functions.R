options(dplyr.summarise.inform = FALSE) # keep dplyr messages quiet 

#------------------------------------------
# Function that performs non-parametric bootstrap of POBASCAM data and then calculates
# the net health benefit (NHB) for each resample so that we can get confidence intervals
#------------------------------------------
npBS_NHB <- function(dat, bmd.strategy, normal.strategy, n.reps, timesinceonset, trichot, sample.vink=F, bootstrap=T){
  if (!bootstrap) n.reps <- 1
  nhbs <- rep(NA, n.reps) # make empty vector to store bootstrap estimates 
  
  for (i in 1:n.reps){
    cumulativeprogcin3for1618_lifetime <- 0.5
    
    ### cumulative progprob cin3 for HPV16/18 lifetime translates into cumulative progprob cin23 lifetime by multiplying by (1-0.2) = 0.8
    ## cumulative progprob cin3 lifetime of 0.5 corresponds to 0.5*0.8 = 0.4 which corresponds with article by Vink et al.
    ## However, this ignores the relation between progprob and timesinceonsetcin23. It is fair to say that cumulativeprogcin3for1618_lifetime lies between 0.5 and 1.
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
    
    # resample from POBASCAM with replacement 
    if (bootstrap){
      resample_dat <- dat[sample(1:dim(dat)[1], replace = T),] 
    }else{
      resample_dat <- dat
    }
    
    pob_params <- pobInput(resample_dat) # calculate POBASCAM parameters from the resampled data set (proportion of regressive CIN)
    
    # calculate the costs and QALY's for the resampled data using only the intervention group 
    temp_res <- computeCostsandQALYloss(resample_dat[resample_dat$grp =='i',], bmd.strategy, normal.strategy, pob_params, probcancerParams, trichot) 
    
    nhbs[i] <- -temp_res[2] - temp_res[1]/ ICERthreshold # compute the NHB
  }
  print(temp_res)
  nhb <- mean(nhbs) # take mean of all the NHB's 
  nhb.lower <- quantile(nhbs, 0.025) # take 2.5th percentile as lower bound for confidence interval
  nhb.upper <- quantile(nhbs, 0.975) # take 97.5th percentile as uppper bound for confidence interval
  return(c(bmd.strategy, normal.strategy, timesinceonset, nhb, nhb.lower, nhb.upper))
}


#------------------------------------------
# Function that computes costs and QALYS for a specific strategy (edited from Hans' code)
#------------------------------------------
# --> Input parameters from "CEA_parameters.R" file 
# Costrepeattest,CostnoCIN2,CostCIN2,CostCIN3,CostCancer,Costdeath,disce,discc,
# LossCIN23,Lossdiagnosiscancer,Lossmonitor,Losspreterminal,Lossterminal,Numsubjects,
# Durationcancer,Maxlifeyearslost,Notdeathfromothercauses5yr

# --> Input parameters that are obtained from POBASCAM data (recalculated for each resampled dataset)
# Probcancersdetected_1618,Probcancersdetected_no1618,ProbregressiveCIN2_1618,ProbregressiveCIN2_no1618,

# --> Input parameters that are obtained by running the most aggressive strategy (refer everyone)
# MaxnumberCIN2round1_1618,MaxnumberCIN2round1_no1618,MaxnumberCIN3round1_1618,MaxnumberCIN3round1_no1618,

# --> Input parameters that are obtained for the specific strategy under consideration by running the inputNumbers() function on the resampled data
# CIN2round1_1618,CIN2round1_no1618,CIN3round1_1618,CIN3round1_no1618,NoCIN2round1, 
# Firsttest,Repeattest,Probcancer5years,Survivalyr,Survivalyr2,

computeCostsandQALYloss <- function(dat, bmd.strategy, normal.strategy, pobParams, probCancerParams, trichot=F)   {
  ProbregressiveCIN2_no1618 <- pobParams$ProbregressiveCIN2_no1618
  ProbregressiveCIN2_1618 <- pobParams$ProbregressiveCIN2_1618
  
  Probcancer5years_1618 <- probCancerParams[1]
  Probcancer5years_no1618 <- probCancerParams[2]
  
  Numsubjects <- dim(dat)[1]
  referAll_res <- unname(inputNumbers(dat, "referAll", "referAll", trichot=F)) # get input numbers for most aggressive strategy (used for maximum number of cases)
  
  MaxnumberCIN2round1_1618 <- referAll_res[1] 
  MaxnumberCIN2round1_no1618 <- referAll_res[2] 
  MaxnumberCIN3round1_1618 <- referAll_res[3] 
  MaxnumberCIN3round1_no1618 <- referAll_res[4] 
  
  input_data <- unname(inputNumbers(dat, bmd.strategy, normal.strategy, trichot)) # get input data for the specified strategies at the start of the function 
  CIN2round1_1618 <- input_data[1] # assign values 
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
  CIN3round1 <- (1 - Probcancersdetected_1618)*CIN3round1_1618 + (1 - Probcancersdetected_1618)*CIN3round1_no1618
  Cancerround1 <- Probcancersdetected_1618*CIN3round1_1618 + Probcancersdetected_no1618*CIN3round1_no1618
  
  Survivalyr <- c(1, Survivalyr)
  Survivalyr2 <- c(1, Survivalyr2)


  Cost <- NoCIN2round1*CostnoCIN2 + CIN2round1*CostCIN2 + CIN3round1*CostCIN3 + Firsttest*Costfirsttest + Repeattest*Costrepeattest + Cancerround1*Costcancer +
    (CIN2round2*CostCIN2 + CIN3round2*CostCIN3)/discc^5 +
    Notdeathfromothercauses5yr*Costcancer*(Cancerround2old)/discc^5 +
    Notdeathfromothercauses5yr*Costcancer*(Cancerround2new)/discc^5
  for (i in 0:9){
    Cost <- Cost + Cancerround1*(Survivalyr[i + 1]-Survivalyr[i + 2])*Costdeath/discc^(i + .5) +
      Notdeathfromothercauses5yr*Cancerround2old*(Survivalyr2[i + 1]-Survivalyr2[i + 2])*Costdeath/discc^(5 + i +.5) +
      Notdeathfromothercauses5yr*Cancerround2new*(Survivalyr[i + 1]-Survivalyr[i + 2])*Costdeath/discc^(5 + i + .5)
  }

  QALYloss <- (CIN2round1+CIN3round1)*0.5*LossCIN23
  QALYloss <- QALYloss + (CIN2round2+CIN3round2)*0.5*LossCIN23/disce^5
  # cancer round1 ----
  QALYloss <- QALYloss +
    Probsurvival*Cancerround1*0.5*(Lossdiagnosiscancer) +
    Probsurvival*Cancerround1*3.5*(Lossmonitor)/disce^2.25 +
    sum(sapply(0:9, function(i) (Survivalyr[i+1]-Survivalyr[i+2])*Cancerround1*0.5*Lossterminal/disce^(i+0.5))) +
    sum(sapply(0:9, function(i) (Survivalyr[i+1]-Survivalyr[i+2])*Cancerround1*(i+0.125)*Losspreterminal/disce^(i+0.125) )) +
    sum(sapply(0:9, function(i) (Survivalyr[i+1]-Survivalyr[i+2])*Cancerround1*((1 - (1/disce)^(Maxlifeyearslost-i))/(1-(1/disce)))/disce^(i+0.5)))

  # cancer round2 old ----
  QALYloss <- QALYloss +
    Notdeathfromothercauses5yr*Probsurvival2*Cancerround2old*0.5*(Lossdiagnosiscancer)/disce^5 +
    Notdeathfromothercauses5yr*Probsurvival2*Cancerround2old*3.5*(Lossmonitor)/disce^7.25 +
    sum(sapply(0:9, function(i) Notdeathfromothercauses5yr*(Survivalyr2[i+1]-Survivalyr2[i+2])*Cancerround2old*0.5*Lossterminal/disce^(i+5.5))) +
    sum(sapply(0:9, function(i) Notdeathfromothercauses5yr*(Survivalyr2[i+1]-Survivalyr2[i+2])*Cancerround2old*(i+0.125)*Losspreterminal/disce^(i+5.125))) +
    sum(sapply(0:9, function(i) Notdeathfromothercauses5yr*(Survivalyr2[i+1]-Survivalyr2[i+2])*Cancerround2old*((1 - (1/disce)^(Maxlifeyearslost-5-i))/(1-(1/disce)))/disce^(i+5.5)))

  # cancer round2 new ----
  QALYloss <- QALYloss +
    Notdeathfromothercauses5yr*Probsurvival*Cancerround2new*0.5*(Lossdiagnosiscancer)/disce^5 +
    Notdeathfromothercauses5yr*Probsurvival*Cancerround2new*3.5*(Lossmonitor)/disce^7.25 +
    sum(sapply(0:9, function(i) Notdeathfromothercauses5yr*(Survivalyr[i+1]-Survivalyr[i+2])*Cancerround2new*0.5*Lossterminal/disce^(i+5.5))) +
    sum(sapply(0:9, function(i) Notdeathfromothercauses5yr*(Survivalyr[i+1]-Survivalyr[i+2])*Cancerround2new*(i+0.125)*Losspreterminal/disce^(i+5.125))) +
    sum(sapply(0:9, function(i) Notdeathfromothercauses5yr*(Survivalyr[i+1]-Survivalyr[i+2])*Cancerround2new*((1 - (1/disce)^(Maxlifeyearslost-5-i))/(1-(1/disce)))/disce^(i+5.5)))

  return(c(cost=Cost,qaly=QALYloss)) 
}

#------------------------------------------
# Function to calculate crude number of referrals for a specific strategy based on cytology & genotype result 
#------------------------------------------
# possible cytology strategies are: referAll, refer1618, refer7, repeat, repeat7, nextround
crudeCases <- function(data, strategy){
  if (strategy == "referAll"){
    baseline.referred <- data$id
    repeat.referred <- c()
    repeat.tested <- c()
    
  }else if (strategy =='refer1618'){
    baseline.referred <- data$id[which(data$hpv.1618pos)]
    repeat.tested <- data$id[which(!data$hpv.1618pos)]
    repeat.referred <- data$id[which(!data$hpv.1618pos & data$first.repeat.cyt =='>=BMD')]
    
  }else if (strategy =='repeat'){
    baseline.referred <- c()
    repeat.tested <- data$id
    repeat.referred <- data$id[which(data$first.repeat.cyt =='>=BMD')]
    
  }else if (strategy == "repeat7"){
    baseline.referred <- c()
    repeat.tested <- data$id[which(data$hpv.genotype.pos)]
    repeat.referred <- data$id[which(data$hpv.genotype.pos & data$first.repeat.cyt =='>=BMD')]
    
  }else if (strategy == 'nextround'){
    baseline.referred <- c()
    repeat.referred <- c()
    repeat.tested <- c()
    
  }else if (strategy == 'refer7'){
    baseline.referred <- data$id[which(data$hpv.genotype.pos)]
    repeat.referred <- data$id[which(!data$hpv.genotype.pos & data$first.repeat.cyt =='>=BMD')]
    repeat.tested <- data$id[which(!data$hpv.genotype.pos)]
  }
  first.round.referred <- c(baseline.referred, repeat.referred)
  results <- list(baseline.referred, repeat.referred, first.round.referred, repeat.tested)
  names(results) <- c("baseline.referred", "repeat.referred", "first.round.referred", "repeat.tested")
  
  return(results)
}

#------------------------------------------
# Function that combines results (i.e. number referred, number of CIN2/3+ cases) from seperate BMD and normal strategies 
#------------------------------------------
# possible cytology strategies are: referAll, refer1618, refer7, repeat, repeat7, nextround
inputNumbers <-  function(data, bmd.strategy, normal.strategy, trichot=F){
  Numsubjects <- total <- dim(data)[1]

  # calculate numbers for those not referred based on HPV test seperate for each cytology group
  hsil.results <- crudeCases(data[data$baseline.cyt==">BMD", ], "referAll") # all HSIL cases are all immediately referred at baseline 
  
  bmd.results <- crudeCases(data[data$baseline.cyt=="BMD", ], bmd.strategy) # count crude referrals for specific BMD strategy
  normal.results <- crudeCases(data[data$baseline.cyt=="normal", ], normal.strategy) # count crude referrals for specific normal strategy
  
  # add together all who were referred at baseline 
  baseline.referred <- length(c(bmd.results$baseline.referred, normal.results$baseline.referred, hsil.results$baseline.referred)) # get total number of baseline referred
  referred <- c(bmd.results$first.round.referred, normal.results$first.round.referred, hsil.results$first.round.referred) # combine Id's of those who were referred 
  
  CIN3round1_1618 <- sum(data[data$id %in% referred & data$genotype.data =='HPV 16/18', ]$first.round.cin3plus)
  CIN2round1_1618 <- sum(data[data$id %in% referred & data$genotype.data =='HPV 16/18', ]$first.round.cin2plus) - CIN3round1_1618
  
  CIN3round1_no1618 <- sum(data[data$id %in% referred & ! data$genotype.data =='HPV 16/18', ]$first.round.cin3plus)
  CIN2round1_no1618 <- sum(data[data$id %in% referred & ! data$genotype.data =='HPV 16/18', ]$first.round.cin2plus) - CIN3round1_no1618
  
  NoCIN2round1 <- length(referred) - CIN2round1_1618 - CIN3round1_1618 - CIN2round1_no1618 - CIN3round1_no1618
  
  # number who had baseline cytology test:
  Firsttest <- case_when(
    bmd.strategy =='referAll' & normal.strategy =='referAll' ~ 0,  # --> 0 if all were referred because then no test is needed
    bmd.strategy =='refer1618' & normal.strategy=='refer1618' ~ dim(data[data$genotype.data!= 'HPV 16/18',])[1],
    bmd.strategy =='refer7' & normal.strategy=='refer7' ~ dim(data[!data$hpv.genotype.pos,])[1],
    T ~ dim(data)[1]  # --> otherwise everyone who wasn't positive for the genotype gets baseline cytology 
  )
  # count number who had a repeat test (only normal or bmd who were not refered based on HPV genotype)
  Repeattest <- length(c(bmd.results$repeat.tested, normal.results$repeat.tested, hsil.results$repeat.tested))
  
  res <- c(CIN2round1_1618, CIN2round1_no1618, CIN3round1_1618, CIN3round1_no1618, NoCIN2round1, Firsttest, Repeattest)
  names(res) <- c("CIN2round1_1618", "CIN2round1_no1618", "CIN3round1_1618", "CIN3round1_no1618", "NoCIN2round1", "Firsttest", "Repeattest")
  return(res)
}



#------------------------------------------
# Function that creates the input data from POBASCAM required in the cost/QALY function, this function will be called in the bootstrap function 
#------------------------------------------
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
# Function that organise the output from the list of each strategy and time point into a nice data frame and take difference from most aggressive strategy (refer all)
#------------------------------------------
organiseAndSave <- function(data, name){
  res <- lapply(1:length(bmd.strategies), function(x) {matrix(unlist(data[[x]]), nrow=length(onsetTimes), byrow=T) %>% 
      as.data.frame() %>%
      `colnames<-`(c("BMDStrategy", "NormalStrategy", "time", "NHB", "NHBlower", "NHBupper")) %>%
      mutate_at(c("NHB", "NHBlower", "NHBupper"), as.numeric) }) %>% 
    bind_rows() %>% 
    group_by(time) %>% 
    mutate_if(is.numeric, list(d = ~. - first(.))) %>% 
    ungroup() %>% 
    mutate(time = as.numeric(time), strategy = paste0(BMDStrategy, "-", NormalStrategy)) 
  
  save(res, file=paste0("resultsCEA_", name, "_", file.end, ".RData"))
  
  # organise for excel file for manuscript
  if (name != 'onerun'){
    res  %>% as.data.frame() %>% mutate_if(is.numeric, round, 2) %>%
      mutate(time = as.factor(time), strategy = as.factor(strategy), combined = paste0(NHB_d, " (", NHBlower_d, ", ", NHBupper_d, ")")) %>%
      select(c(combined, time, strategy)) %>%
      reshape(., idvar = "strategy", timevar = "time", direction = "wide") %>% write.xlsx(paste0("resultsCEA_", name, "_", file.end, ".xlsx"))
  }else{
    res  %>% as.data.frame() %>% mutate_if(is.numeric, round, 2) %>%
      mutate(time = as.factor(time), strategy = as.factor(strategy)) %>%
      select(c(NHB_d, time, strategy)) %>%
      reshape(., idvar = "strategy", timevar = "time", direction = "wide")  %>% write.xlsx(paste0("resultsCEA_", name, "_", file.end, ".xlsx"))
  }
  return(res)
}



#------------------------------------------
# Function that creates the plot of the NHB for each time point for each strategy 
#------------------------------------------
plotCEA <- function(resultsCEA, bootstrap=T, sample.vink=F){
  max.y <- max(resultsCEA$NHB_d)
  title.temp <- ifelse(bootstrap, 
                       'Net health benefit (NHB) of genotyping strategies for triage of HPV-positive women (using bootstrapping for NHB)', 
                       'Net health benefit (NHB) of genotyping strategies for triage of HPV-positive women')
  subtitle.temp <- ifelse(sample.vink, 
                          "Using sampled estimates from Vink et al. (2013). Strategies are coded 'BMD-Normal' and all >BMD are immediately referred",
                          "Using fixed estimates from Vink et al. (2013). Strategies are coded 'BMD-Normal' and all >BMD are immediately referred")
  plot <- ggplot(resultsCEA, aes(x=time, y=NHB_d, group=strategy, col=strategy)) + theme_classic() + 
    labs(title=title.temp,subtitle=subtitle.temp) +
    geom_line() + xlab("Time since CIN2/3 onset (years)") + ylab("NHB (difference to reference strategy: refer all)") + 
    coord_cartesian(clip = "off") + theme(legend.position = 'none') +  theme(plot.margin = unit(c(0.1, 4, 0.1, 0.1), "cm")) + 
    geom_text_repel(data = . %>% filter(time == 10),aes(label = strategy, col=strategy), 
                    direction = "y", segment.linetype = "dotted", segment.color='grey',
                    hjust = 0, segment.size = 0.5, xlim = c(10.05, Inf)) + 
    scale_x_continuous(expand = c(0, 0), limits=c(0, 10), breaks=seq(0, 10, 1))  + 
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, round(max(resultsCEA_onerun$NHB_d)/2.5)* 2.5, by = 2.5), 
                       limits=c(-0.05, ceiling(max(resultsCEA_onerun$NHB_d)) + 0.5)) 
  return(plot)
}