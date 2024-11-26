# CEA_triage_Feb2024


########### GENERAL PARAMETERS ####################
#### Cost data from Sollie et al. MedArxiv and Jansen et al. BJOG 2020 ###
Costfirsttest <- 26 / 106.16 * 121.42  #https://opendata.cbs.nl/#/CBS/nl/dataset/83131NED/table indexed to 2022
Costrepeattest <- 53/ 106.16 * 121.42 #https://opendata.cbs.nl/#/CBS/nl/dataset/83131NED/table indexed to 2022
CostnoCIN2 <- 609
CostCIN2 <- 1820
CostCIN3 <- 2183
Costcancer <- 10364
Costdeath <- 25392
###### Utilities from McDonald et al. Cancer Causes Control ###
LossCIN23 <- 0.07
Lossdiagnosiscancer <- 0.43
Lossmonitor <- 0.2
Losspreterminal <- 0.75
Lossterminal <- 0.93
### Discount rates the Netherlands ####
disce <- 1.015
discc <- 1.03
###### Survival rate age-adjusted according to CIN2/3 age distribution ####
Survivalyr <- rep(0,10)
Survivalyr[1] <- 0.951
Survivalyr[2] <- 0.906
Survivalyr[3] <- 0.871
Survivalyr[4] <- 0.857
Survivalyr[5] <- 0.843
Survivalyr[6] <- 0.833
Survivalyr[7] <- 0.823
Survivalyr[8] <- 0.820
Survivalyr[9] <- 0.810
Survivalyr[10] <- 0.810
#### Define relative survival 5 years later for estimating stage shift effect
relativesurvivalafter5years <- .9
#### A value of .9 corresponds with survival per age group, 5 year older, relative survival 0.9
### But effect may be stronger. However, unlikely to be stronger than 0.7 because that corresponds with 
### relative survival when shifting from stage 1 to 2. I suggest putting relative survival between 0.7 and 1 with .9 as base-case.

Survivalyr2 <- Survivalyr*relativesurvivalafter5years
#### Duration cancer until utility is equal to maximum value is set to 4 years, see McDonald et al. ####
Durationcancer <- 4
##### Death from other causes set at 0, value is 0.004 ###
Notdeathfromothercauses5yr <- 1
### Based on conditional life expectancy of women aged 40 ####
Maxlifeyearslost <- 43.35   


ICERthreshold <- 20000

################### DATA progression / missed cancer NATIONAL registry ########## 
### Cancersdetected : number of cancers versus number of CIN2/3 in national program, year 2017 (Inturrisi et al. Lancet Regional Health 2021)
Probcancersdetected_1618 <- 53 / 1038
Probcancersdetected_no1618 <- 32 / 847

cumulativeprogcin3for1618_lifetime <- 0.5