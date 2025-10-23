# required libraries 
library(purrr)
library(stringr)
library(dplyr)
library(foreign)
library(haven)
library(lubridate)
library(Rlab)
#-----

# Create a data frame which contains all the required information: baseline cytology, HPV genotypes, 
# repeat cytology and second round HPV test result (easier because POBASCAM has so many columns)
create.referral.data<- function(){
  pob.full <- read.spss("/Users/kelsi/Desktop/PhD/Projects/Colposcopy Referral Rates/Analysis/Data/44102_10yrs_fup with age groups.sav", to.data.frame=T) # load the data 
  pob.full <- pob.full[pob.full$lftd < 59 & trimws(pob.full$cyt1_dia) != "Pap 0", ] 
  pob.full <- pob.full[pob.full$hr %in% c("+", "-"), ] # remove "bad" HPV test results 
  pob.full$cyt1_hpvtype <- trimws(pob.full$cyt1_hpvtype) # remove white space from HPV type columns
  pob.pos <- pob.full[pob.full$cyt1_hpv %in% c("*h+", "hr+"),] # filter to include only HPV positives 
  
  pob.pos$cyt5j_dia <- trimws(pob.pos$cyt5j_dia)
  id <- pob.pos$bonr # store the id numbers 
  
  grp <- pob.pos$gr # POBASCAM intervention or control group 
  
  first.screen <- pob.pos$previoussmear
  # baseline cytology result 
  baseline.cyt <- case_when(
    grepl(paste(c('Pap 1'), collapse="|"), pob.pos$cyt1_dia) ~"normal",
    grepl(paste(c('Pap 3a2', 'Pap 3b', '4', '5'), collapse="|"), pob.pos$cyt1_dia) ~ ">BMD",
    grepl(paste(c('Pap 2', 'Pap 3a', 'Pap 3a1'), collapse="|"), pob.pos$cyt1_dia) ~"BMD")
  
  # date of baseline cytology (entry date)
  baseline.date <- as.Date(as.POSIXct(pob.pos$cyt1_dat, origin='1582-10-14', tz='GMT'))
  
  baseline.age <- pob.pos$lftd  # baseline age 
  
  hpv.genotype <- pob.pos$cyt1_hpvtype # HPV genotype 
  hpv.1618pos <- str_detect(hpv.genotype, "16|18") # indicator for HPV16/18+
  hpv.genotype.pos <- str_detect(hpv.genotype, "16|18|31|33|45|52|58") # indicator for HPV16/18/31/33/45/52/58+
  
  first.repeat.cyt <- pob.pos[, c("cyt06_dia", "cyt18_dia", "cyt3j_dia")]
  
  # trim white space in all the rows
  first.repeat.cyt <-  as.data.frame(apply(first.repeat.cyt, 2, function(x) trimws(x)))
  
  # fill Pap 0 with blank space because Pap 0 result got repeat smear due to poor quality smear  
  first.repeat.cyt[first.repeat.cyt =='Pap 0' | is.na(first.repeat.cyt)] <- ""
  # find the index of the first repeat cytology by finding the first non-empty entry in each row 
  first.cyt.indexes <- apply(first.repeat.cyt, 1, function(x) detect_index(x, function(i) trimws(i)!=""))
  
  # create empty vector which will store the values of the first repeat cytology 
  first.repeat.cyt.result <- rep(NA, dim(pob.pos)[1])
  
  for (i in 1:dim(pob.pos)[1]){
    if (first.cyt.indexes[i] ==0){
      # if the index of first cytology = 0 then that means there were no cytology results 
      first.repeat.cyt.result[i] <- 'Missing'
    }else if (first.cyt.indexes[i] >0){
      # store the cytology result using the row number and index of the repeat smear 
      first.repeat.cyt.result[i] <- first.repeat.cyt[i,first.cyt.indexes[i]] #levels(first.repeat.cyt[,1])[]
    }
  }
  
  # store the 'category' of the first repeat cytology result (NILM/BMD/HSIL/Missing)
  # first.repeat.cyt.cat <- case_when(
  #   trimws(first.repeat.cyt.result) %in% c("Pap 0", "Pap 1") ~ "NILM",
  #   trimws(first.repeat.cyt.result) > "Pap 1" ~ "ASC-US/LSIL+",
  #   TRUE ~'Missing')
  
  first.repeat.cyt.cat <- case_when(
    grepl(paste(c('Pap 1', 'Pap 0'), collapse="|"), first.repeat.cyt.result) ~"normal",
    grepl(paste(c('Pap 2', 'Pap 3a', 'Pap 3a1', 'Pap 3a2', 'Pap 3b', 'Pap 4', 'Pap 5'), collapse="|"), first.repeat.cyt.result) ~">=BMD",
    TRUE ~ "Missing")
  
  # keep the dates of the repeat cytology information
  first.cyt.dat <- pob.pos[, c("cyt06_dat", "cyt18_dat", "cyt3j_dat")] 
  # convert the dates from SAS format to dd-MM-YYYY format (SAS stores dates with the origin 1582-10-14)
  first.cyt.dat$cyt06_dat <-  as.Date(as.POSIXct(first.cyt.dat$cyt06_dat, origin='1582-10-14', tz='GMT'))
  first.cyt.dat$cyt18_dat <-  as.Date(as.POSIXct(first.cyt.dat$cyt18_dat, origin='1582-10-14', tz='GMT'))
  first.cyt.dat$cyt3j_dat <-  as.Date(as.POSIXct(first.cyt.dat$cyt3j_dat, origin='1582-10-14', tz='GMT'))
  
  # create empty vector which will store the dates of the first repeat cytology 
  first.repeat.cyt.date <- rep(NA, dim(pob.pos)[1])
  for (i in 1:dim(pob.pos)[1]){
    if (first.cyt.indexes[i] ==0){
      # if the index =0 then they didnt have repeat cytology so date has to be NA
      first.repeat.cyt.date[i] <- NA
    }else if (first.cyt.indexes[i] >0){
      # store the date of first repeat cytology using the row number and index of the repeat smear found above
      first.repeat.cyt.date[i] <- first.cyt.dat[i,first.cyt.indexes[i]]
    }
  }
  # ensure that the dates are in the right format after adding some NA values 
  first.repeat.cyt.date <- as.Date(first.repeat.cyt.date, origin='1970-01-01')
  
  #--- first repeat HPV
  first.repeat.hpv <- pob.pos[, c("cyt06_hr", "cyt18_hr", "cyt3j_hr")]
  
  # trim white space in all the rows
  first.repeat.hpv <-  as.data.frame(apply(first.repeat.hpv, 2, function(x) trimws(x)))
  first.repeat.hpv[is.na(first.repeat.hpv)] <- ""
  # find the index of the first repeat cytology by finding the first non-empty entry in each row 
  first.hpv.indexes <- apply(first.repeat.hpv, 1, function(x) detect_index(x, function(i) trimws(i)!=""))
  
  # create empty vector which will store the values of the first repeat cytology 
  first.repeat.hpv.result <- rep(NA, dim(pob.pos)[1])
  
  for (i in 1:dim(pob.pos)[1]){
    if (first.hpv.indexes[i] ==0){
      # if the index of first cytology = 0 then that means there were no cytology results 
      first.repeat.hpv.result[i] <- 'Missing'
    }else if (first.hpv.indexes[i] >0){
      # store the cytology result using the row number and index of the repeat smear 
      first.repeat.hpv.result[i] <- first.repeat.hpv[i,first.hpv.indexes[i]]
    }
  }
  
  # store the 'category' of the first repeat cytology result (NILM/BMD/HSIL/Missing)
  first.repeat.hpv.cat <- case_when(
    first.repeat.hpv.result =="+" ~ 1,
    first.repeat.hpv.result > "-" ~ 0,
    TRUE ~0)
  
  # keep the dates of the repeat cytology information
  first.hpv.dat <- pob.pos[, c("cyt06_dat", "cyt18_dat", "cyt3j_dat")] 
  # convert the dates from SAS format to dd-MM-YYYY format (SAS stores dates with the origin 1582-10-14)
  first.hpv.dat$cyt06_dat <-  as.Date(as.POSIXct(first.hpv.dat$cyt06_dat, origin='1582-10-14', tz='GMT'))
  first.hpv.dat$cyt18_dat <-  as.Date(as.POSIXct(first.hpv.dat$cyt18_dat, origin='1582-10-14', tz='GMT'))
  first.hpv.dat$cyt3j_dat <-  as.Date(as.POSIXct(first.hpv.dat$cyt3j_dat, origin='1582-10-14', tz='GMT'))
  
  # create empty vector which will store the dates of the first repeat cytology 
  first.repeat.hpv.date <- rep(NA, dim(pob.pos)[1])
  for (i in 1:dim(pob.pos)[1]){
    if (first.hpv.indexes[i] ==0){
      # if the index =0 then they didnt have repeat cytology so date has to be NA
      first.repeat.hpv.date[i] <- NA
    }else if (first.hpv.indexes[i] >0){
      # store the date of first repeat cytology using the row number and index of the repeat smear found above
      first.repeat.hpv.date[i] <- first.hpv.dat[i,first.hpv.indexes[i]]
    }
  }
  # ensure that the dates are in the right format after adding some NA values 
  first.repeat.hpv.date <- as.Date(first.repeat.hpv.date, origin='1970-01-01')
  
  # columns of repeat cytology results 
  second.repeat.cyt <- pob.pos[, c("c5t06_dia", "c5t18_dia", "c5t3j_dia")]
  
  # trim white space in all the rows
  second.repeat.cyt <-  as.data.frame(apply(second.repeat.cyt, 2, function(x) trimws(x)))
  
  # fill Pap 0 with blank space because Pap 0 result got repeat smear due to poor quality smear  
  second.repeat.cyt[second.repeat.cyt =='Pap 0' | is.na(second.repeat.cyt)] <- ""
  
  # find the index of the first repeat cytology in the SECOND round by finding the first non-empty entry in each row 
  second.cyt.indexes <- apply(second.repeat.cyt, 1, function(x) detect_index(x, function(i) trimws(i)!=""))
  
  # create empty vector which will store the values of the first repeat cytology 
  second.repeat.cyt.result <- rep(NA, dim(pob.pos)[1])
  
  for (i in 1:dim(pob.pos)[1]){
    if (second.cyt.indexes[i] ==0){
      # if the index of first cytology = 0 then that means there were no cytology results 
      second.repeat.cyt.result[i] <- 'Missing'
    }else if (second.cyt.indexes[i] >0){
      # store the cytology result using the row number and index of the repeat smear 
      second.repeat.cyt.result[i] <- second.repeat.cyt[i,second.cyt.indexes[i]] #levels(second.repeat.cyt[,1])[second.repeat.cyt[i,second.cyt.indexes[i]]]
    }
  }
  
  # store the 'category' of the first repeat cytology result (NILM/BMD/HSIL/Missing)
  second.repeat.cyt.cat <- case_when(
    second.repeat.cyt.result =="Pap 1" ~ "normal",
    second.repeat.cyt.result > "Pap 1" ~ ">=BMD",
    TRUE ~'Missing')
  second.repeat.cyt <- ifelse(second.repeat.cyt.cat ==">=BMD", TRUE, FALSE)
  
  
  # keep the dates of the repeat cytology information
  second.cyt.dat <- pob.pos[, c("c5t06_dat", "c5t18_dat", "c5t3j_dat")] 
  # convert the dates from SAS format to dd-MM-YYYY format (SAS stores dates with the origin 1582-10-14)
  second.cyt.dat$c5t06_dat <-  as.Date(as.POSIXct(second.cyt.dat$c5t06_dat, origin='1582-10-14', tz='GMT'))
  second.cyt.dat$c5t18_dat <-  as.Date(as.POSIXct(second.cyt.dat$c5t18_dat, origin='1582-10-14', tz='GMT'))
  second.cyt.dat$c5t3j_dat <-  as.Date(as.POSIXct(second.cyt.dat$c5t3j_dat, origin='1582-10-14', tz='GMT'))
  
  # create empty vector which will store the dates of the first repeat cytology 
  second.repeat.cyt.date <- rep(NA, dim(pob.pos)[1])
  for (i in 1:dim(pob.pos)[1]){
    if (second.cyt.indexes[i] ==0){
      # if the index =0 then they didnt have repeat cytology so date has to be NA
      second.repeat.cyt.date[i] <- NA
    }else if (second.cyt.indexes[i] >0){
      # store the date of first repeat cytology using the row number and index of the repeat smear found above
      second.repeat.cyt.date[i] <- second.cyt.dat[i,second.cyt.indexes[i]]
    }
  }
  # ensure that the dates are in the right format after adding some NA values 
  second.repeat.cyt.date <- as.Date(second.repeat.cyt.date, origin='1970-01-01')
  
  
  second.round.hpv <- trimws(pob.pos$cyt5j_hr) # store first HPV test of the second round 
  
  ue.list <- c("UE", "UE No CIN", "UE CIN I", "UE Met.Sqm")
  ue.detected <- ifelse(trimws(pob.pos$hist_dia) %in% ue.list, 1, 0)
  
  # store histology results that are classified as CIN3+
  cin3plus <- c("Aden.Ca", "Aden.Cis", "CIN III", "Plav.Ca", "UE CIN III", "ACIS", "ACIS/CIN-3",
                "PlavCa", "UE AdenCa", "UE Carcino", "UE PlavCa", "AdenCa", "UE ACIS")
  
  cin2plus <- c("Aden.Ca", "Aden.Cis", "CIN III", "Plav.Ca", "UE CIN III", "ACIS", "ACIS/CIN-3",
                "PlavCa", "UE AdenCa", "UE Carcino", "UE PlavCa", "AdenCa", "UE ACIS",
                "CIN II", "UE CIN II")
  # check whether any of the histology results reported are CIN3+
  
  # cin3+ cases & dates 
  cin3plus.cases.table <- ifelse(trimws(pob.pos$hist_dia) %in% cin3plus & pob.pos$hist_tim/365 < 9, 1, 
                                 ifelse(trimws(pob.pos$hist_dia) =='', "missing", 0))
  
  cin3plus.cases <- ifelse(trimws(pob.pos$hist_dia) %in% cin3plus & pob.pos$hist_tim/365 < 9, 1, 0)
  # store the time (in years) of the HPV test at the start of the 2nd round so we can check when the CIN3+ developed
  # --> if there was no test at the start of the second round (NA value) then we set the time to 5 years 
  #second.round.times <- ifelse(is.na(pob.pos$cyt5j_tim/365.25), 5, pob.pos$cyt5j_tim/365.25)
  # keep only CIN3+ cases that were detected in the first round (0-4 years)
  first.round.cin3plus <- ifelse(pob.pos$hist_tim/365  < 4,  cin3plus.cases, 0)
  first.round.cin3plus[is.na(first.round.cin3plus)] <- 0 # no histology result means that test was not done because dr saw that it was not needed
  # keep only CIN3+ cases that were detected in the second round (4-9 years)
  second.round.cin3plus <- ifelse(pob.pos$hist_tim/365 >= 4 & pob.pos$hist_tim/365 < 9, cin3plus.cases, 0)
  second.round.cin3plus[is.na(second.round.cin3plus)] <- 0 # same reason as above about NA histology results 
  hist.time <- ifelse(cin3plus.cases.table=='missing', NA, pob.pos$hist_tim)
  
  # cin2+ cases & dates 
  cin2plus.cases <- ifelse(trimws(pob.pos$hist_dia) %in% cin2plus & pob.pos$hist_tim/365 < 9, 1, 0)
  # keep only CIN3+ cases that were detected in the first round (0-4 years)
  first.round.cin2plus <- ifelse(pob.pos$hist_tim/365 < 4,  cin2plus.cases, 0)
  first.round.cin2plus[is.na(first.round.cin2plus)] <- 0 # no histology result means that test was not done because dr saw that it was not needed
  # keep only CIN3+ cases that were detected in the second round (4-9 years)
  second.round.cin2plus <- ifelse(pob.pos$hist_tim/365 >= 4 & pob.pos$hist_tim/365 < 9, cin2plus.cases, 0)
  second.round.cin2plus[is.na(second.round.cin2plus)] <- 0 # same reason as above about NA histology results 
  
  
  second.round.cyt <- pob.pos$cyt5j_dia
  
  second.round.cyt.cat <- case_when(
    grepl(paste(c('Pap 0','Pap 1'), collapse="|"), second.round.cyt) ~"normal",
    grepl(paste(c('Pap 3a2', 'Pap 3b', '4', '5'), collapse="|"), second.round.cyt) ~ ">BMD",
    grepl(paste(c('Pap 2', 'Pap 3a', 'Pap 3a1'), collapse="|"),second.round.cyt) ~"BMD", 
    TRUE~"Missing")
  
  # hierarchical genotype grouping 
  genotype.data <- case_when(
    hpv.1618pos ~ "HPV 16/18",
    hpv.genotype.pos ~ "HPV 31/33/45/52/58", # if they weren't 16/18 then they would be caught in this group
    TRUE ~ "Other"
  )
  age.group <- pob.pos$lftdgr
  # finally combine everything in a dataframe called 'referral.data' to be used to check screening programs 
  referral.data <- data.frame(id, grp, first.screen, baseline.age, baseline.date, baseline.cyt, hpv.1618pos,
                              hpv.genotype.pos, first.repeat.cyt.date, first.repeat.cyt = first.repeat.cyt.cat, 
                              first.repeat.hpv.date, first.repeat.hpv = first.repeat.hpv.cat, 
                              second.round.hpv, cin3plus.cases, first.round.cin3plus, second.round.cin3plus, ue.detected,
                              hist.time,
                              second.round.cyt = second.round.cyt, second.round.cyt.cat = second.round.cyt.cat,
                              second.round.hpv1618 = str_detect(pob.pos$cyt5j_hpvtype, "16|18"),
                              second.round.hpv.genotype = str_detect(pob.pos$cyt5j_hpvtype, "16|18|31|33|45|52|58"),
                              second.repeat.cyt = second.repeat.cyt,
                              second.repeat.cyt.cat = second.repeat.cyt.cat,
                              second.repeat.cyt.date, 
                              cin2plus.cases, first.round.cin2plus, second.round.cin2plus, genotype.data, age.group, cin3plus.cases.table,
                              hist.dia = pob.pos$hist_dia)
  return(referral.data)
}

referral.data.raw <- create.referral.data()


#-----
# flowchart calculations... 
hist.groups <- function(x){
  return(case_when(trimws(x) %in% cin3plus.codes ~ "CIN3+",
                   trimws(x) %in% c("CIN II", "UE CIN II") ~ "CIN2",
                   trimws(x) =='' ~ "Missing",
                   T ~ "<=CIN1"))
}

cin3plus.codes <- c("Aden.Ca", "Aden.Cis", "CIN III", "Plav.Ca", "UE CIN III", "ACIS", "ACIS/CIN-3",
                    "PlavCa", "UE AdenCa", "UE Carcino", "UE PlavCa", "AdenCa", "UE ACIS")


referral.data.raw$baseline.group <- case_when(
  referral.data.raw$baseline.cyt =='BMD' ~'BMD',
  (referral.data.raw$baseline.cyt=='normal' & referral.data.raw$grp=='i') ~ "Normal (i)",
  (referral.data.raw$baseline.cyt=='normal' & referral.data.raw$grp=='c' & !referral.data.raw$hpv.genotype.pos) ~ "Normal (c)",
  T ~ NA
)

referral.data.raw$hist.groups <- hist.groups(referral.data.raw$hist.dia)
referral.data.raw$first.round.hist <- ifelse(referral.data.raw$hist.time/365.25 <= 4, "first", "second")

# referral.data.raw[! is.na(referral.data.raw$baseline.group),] %>% 
#   group_by(baseline.cyt, baseline.group, first.repeat.cyt) %>% summarise(n=n()) %>% print(n=51)
# 
# 
# referral.data.raw[! is.na(referral.data.raw$baseline.group),] %>% 
#   group_by(baseline.cyt, baseline.group, first.repeat.cyt, first.round.hist, hist.groups) %>% summarise(n=n()) %>% print(n=51)
# 
# referral.data.raw[! is.na(referral.data.raw$baseline.group),] %>% 
#   group_by(baseline.cyt, first.round.hist, hist.groups) %>% summarise(n=n()) %>% print(n=51)

#-----

# impute repeat cytology test for CIN3+ cases in women with BMD cytology that violated protocol (either CIN3+ detected before first repeat test or no repeat test available)
illegal.bmd.cases <- referral.data.raw[referral.data.raw$first.repeat.cyt=='Missing' & referral.data.raw$first.round.cin3plus==1 & referral.data.raw$baseline.cyt =='BMD',]$id
legal.bmd.cases <- referral.data.raw[referral.data.raw$first.repeat.cyt!='Missing' & referral.data.raw$first.round.cin3plus==1 & referral.data.raw$baseline.cyt =='BMD',] # adhere protocol
referral.data.raw[referral.data.raw$id %in% illegal.bmd.cases,]$first.repeat.cyt <- ifelse(rbern(length(illegal.bmd.cases), prop.table(table(legal.bmd.cases$first.repeat.cyt))[1]), ">=BMD", "normal")

illegal.nilm.cases <-  referral.data.raw[referral.data.raw$first.repeat.cyt=='Missing' & referral.data.raw$first.round.cin3plus==1 & referral.data.raw$baseline.cyt =='normal',]$id

# impute second round HPV test for those first round HPV cases that will be missed so we will count them as a referral and detection in the 2nd round
referral.data <- referral.data.raw
referral.data$imputed.2nd.hpv <- rep( '', dim(referral.data)[1])
referral.data$imputed.2nd.hpv[referral.data$first.round.cin3plus ==1] <- '+'

# censor after UE or CIN2:
# if you had CIN2 or UE reported in the 1st round then censor at the time of histology test (make all further tests NA)
referral.data$time2repeat <- difftime(referral.data$first.repeat.cyt.date, referral.data$baseline.date, units = 'days')

referral.data$ue.cin2plus.first.round <- ifelse((referral.data$ue.detected==1 & referral.data$hist.time/365.25 <= 4) | (referral.data$first.round.cin2plus), T, F)

ue.cin2.first.round <- which((referral.data$ue.detected==1 & referral.data$hist.time/365.25 <= 4) | (referral.data$first.round.cin2plus & !referral.data$first.round.cin3plus))

# cens at repeat test (remove repeat and second round tests) if UE or CIN2 was reported before the first repeat test
cens.at.repeat <- which(referral.data[ue.cin2.first.round,]$hist.tim  < referral.data[ue.cin2.first.round,]$time2repeat)
# cens after repeat test (remove second round test) if UE or CIN2 was reported after the first repeat test
cens.after.repeat <- which(referral.data[ue.cin2.first.round,]$hist.tim  >= referral.data[ue.cin2.first.round,]$time2repeat)

ref.data.cens <- referral.data
ref.data.cens[cens.at.repeat,]$first.repeat.cyt <- 'Missing' # remove future test results
ref.data.cens[cens.at.repeat,]$second.round.hpv <- '' # remove second round test results

ref.data.cens[cens.after.repeat,]$second.round.hpv <- '' # remove second round test results
ref.data.cens$grp <- factor(ref.data.cens$grp, levels=c("i", "c"))
ref.data.cens$baseline.cyt <- factor(ref.data.cens$baseline.cyt, levels=c(">BMD", "BMD", "normal"))

#ref.data.intv <- ref.data.cens[ref.data.cens$grp =='i',]

rm(list=setdiff(ls(), "ref.data.cens"))
save(ref.data.cens, file='ref.data.cens.RData')
