library(readxl)

# Paramters

# Discount rates the Netherlands ------------------------------------------
disce <- 1.015
discc <- 1.03
ICERthreshold <- 20000

# Survival data -----------------------------------------------------------
load("~/Projects/Colposcopy CEA/Analysis/data/relative_survivals.RData")
Survivalyr2 <- Survivalyr2[-1] # relative survival 5 years later for estimating stage shift effect
Survivalyr <- Survivalyr[-1]

Notdeathfromothercauses5yr <- 1 # Death from other causes set at 0, value is 0.004
Maxlifeyearslost <- 43.35 # Based on conditional life expectancy of women aged 40

Probcancersdetected_1618 <-  53 / 1038 # Probability of cancer detected, stratified by 16/18 or other 
Probcancersdetected_no1618 <-  19 / 638  
cumulativeprogcin3for1618_lifetime <- 0.5 # cumulative probability of CIN3+ progressing

# Utilities ---------------------------------------------------------------
LossCIN01 <- 0
LossCIN23 <- 0.07
Lossdiagnosiscancer <- 0.43
Lossmonitor <- 0.2
Losspreterminal <- 0.8

# Cost data ---------------------------------------------------------------
# (https://opendata.cbs.nl/#/CBS/nl/dataset/83131NED/table indexed to 2024)
Costfirsttest <- 26 / 106.16 * 130.31  #Jansen et al. BJOG 2020 --> increases when self-sampling 
Costrepeattest <- 53 / 106.16 * 130.31 # Jansen et al. BJOG 2020
CostnoCIN2 <- 599 / 126.09 * 130.31 
CostCIN2 <- 1820 / 126.09 * 130.31
CostCIN3 <- 2183 / 126.09 * 130.31

# cancer/age distribution in NCR data
ageDist_cancer <- read_excel("~/Projects/Colposcopy CEA/Analysis/data/NCR-export-29_01_2025_11_06_33.xlsx", skip = 9)[,c(1:5)]
colnames(ageDist_cancer) <- c("Year", "ageGroup", "Stage", "Dist", "N")
NLperc <- ageDist_cancer %>% filter(ageGroup %in% sapply(4:13, function(i) paste0((i-1)*5, "-", i*5-1))) %>% 
  group_by(ageGroup) %>% summarise(tot = sum(N)) %>% ungroup() %>% mutate(perc = tot/sum(tot)) %>% select(perc)

costsGER <- c(939.67, 2332.42, 2695.65, 2729.11, 2820.06, 2968.87, 2962.30, 2821.22, 2384.68, 1147.39) # cost data from Germany 
Costcancer_indirect <- (sum(costsGER * NLperc) * 0.744/0.711) / 91.59 * 130.31 # indexed from 2010 to 2024  (PPP from https://databank.worldbank.org/source/icp-2021/)
Costcancer_direct <- 8000 / 93.73 * 130.31 # de Kok EJC 2011
Costcancer <- Costcancer_direct + Costcancer_indirect #  15094.61
 
ageDist_death <- read_excel("~/Projects/Colposcopy CEA/Analysis/data/NCR-export-2025_04_07_01_43_21.xlsx", skip = 9)[,c(1:3)] #https://nkr-cijfers.iknl.nl/viewer/sterfte-per-jaar?language=en_GB&viewerId=7dbce7a9-666c-411b-8c24-96c9291c078c
colnames(ageDist_death) <- c("Year", "ageGroup", "N")
NLperc_death <- ageDist_death %>% group_by(ageGroup) %>% summarise(tot = sum(N)) %>% ungroup() %>% mutate(perc = tot/sum(tot)) %>% select(perc)

costsGER_death <- c(1333.31, 3309.52, 3824.91, 3872.38, 4001.44, 4212.58, 4203.26, 4003.08, 3383.67, 1628.06) # cost data from Germany 
Costdeath_indirect <- (sum(costsGER_death * NLperc_death) * 0.744/0.711) / 91.59 * 130.31 # indexed from 2010 to 2024  (PPP from https://databank.worldbank.org/source/icp-2021/)
Costdeath_direct <- 19600 / 93.73 * 130.31 #27249.29
Costdeath <- Costdeath_direct + Costdeath_indirect #  557789.5

# Simulation study for average time detection -----------------------------
# set.seed(1234)
# # stage 1 to stage 2/3, and then stage2/3 to stage 4 if time < 5 years
# move_time1 <- rexp(1000, rate=1/8) # time to move to stage 2/3 from stage 1
# sympt_time1 <- rexp(1000, rate=-log(1-0.075)) # time to get detected by symptoms in stage 1
# time1 <- pmin(move_time1, sympt_time1) # time to first event
# time1[time1 < 5] <- time1[time1 < 5] + pmin(rexp(sum(time1<5), rate=1/4), rexp(sum(time1<5), rate=-log(1-0.4))) # stage2/3 to stage 4 if time < 5 years
# time1[time1>5] <- 5 # all detected at 5 year screen anyway
# 
# # stage 2/3 to 4
# move_time2 <- rexp(1000, rate=1/4) # time to move to stage 4 from stage 2/3
# sympt_time2 <- rexp(1000, rate=-log(1-0.4)) # time to get detected by symptoms in stage 2/3
# time2 <- pmin(move_time2, sympt_time2) # time to first event
# time2[time2>5] <- 5 # all detected at 5 year screen anyway
# 
# load("NCR_stage_props.RData")
#t_reset <- stage_props[1]/(stage_props[1] + stage_props[2])*mean(time1) + stage_props[2]/(stage_props[1] + stage_props[2])*mean(time2) # overall mean time = weighted average = 2.946024
t_reset <- 3 # round the overall mean time to 3