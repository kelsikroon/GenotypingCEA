
load('data/ref.data.cens.RData') # data already organised/prepared for analysis

# 1. TIME BETWEEN FIGO STAGE 1 TO 2/3 CHANGES FROM 8 YEARS TO 4 YEARS
source("code/2_parameters.R") # reset parameters 
Survivalyr2 <- survival5yearsLater(time1to23 = 4)[,2][-1] # update survival 5 years later
res_SA1 <- do.call(rbind, lapply(onsetTimes, function(x) one_rep(ref.data.cens, x)))

# 2. PROGRESSION PROBABILITY CHANGES FROM 0.5 TO 0.7
source("code/2_parameters.R") # reset parameters 
cumulativeprogcin3for1618_lifetime <- 0.7 # update progression parameter
res_SA2 <- do.call(rbind, lapply(onsetTimes, function(x) one_rep(ref.data.cens, x)))

# 3. ADD DISUTILITY FOR CIN0/1 DIAGNOSIS
source("code/2_parameters.R") # reset parameters 
LossCIN01 <- 0.01 # add disutility for CIN0/1
res_SA3 <- do.call(rbind, lapply(onsetTimes, function(x) one_rep(ref.data.cens, x)))

# 4. ICER THRESHOLD CHANGES TO 50,000EUR
source("code/2_parameters.R")
ICERthreshold <- 50000 # update ICER threshold
res_SA4 <- do.call(rbind, lapply(onsetTimes, function(x) one_rep(ref.data.cens, x)))

# 5. ONLY HEALTHCARE COSTS (DIRECT COSTS
source("code/2_parameters.R")
CostnoCIN2 <- 589 / 126.09 * 130.31 # update costs
CostCIN2 <- 1795 / 126.09 * 130.31
CostCIN3 <- 2154 / 126.09 * 130.31
Costcancer <- 8000 / 93.73 * 130.31 # de Kok EJC 2011
Costdeath <- 19600 / 93.73 * 130.31 #27249.29
res_SA5 <- do.call(rbind, lapply(onsetTimes, function(x) one_rep(ref.data.cens, x)))

# 6. 100% PRIMARY SELF-SAMPLING
source("code/2_parameters.R")
Costfirsttest <- 53 / 106.16 * 130.31 # update cost of first test --> increased to same as repeat test when considering self-sampling 
res_SA6 <- do.call(rbind, lapply(onsetTimes, function(x) one_rep(ref.data.cens, x)))

#  Create plots for manuscript
SA1_plot <- plotCEA(res_SA1, cuttime=10, label='short', labelLines = F, sensPlot = T, maximumY =400)
SA2_plot <- plotCEA(res_SA2, cuttime=10, label='short', labelLines = F, sensPlot = T, maximumY =400)
SA3_plot <- plotCEA(res_SA3, cuttime=10, label='short', labelLines = F, sensPlot = T, maximumY =400)
SA4_plot <- plotCEA(res_SA4, cuttime=10, label='short', labelLines = F, sensPlot = T, maximumY =400)
SA5_plot <- plotCEA(res_SA5, cuttime=10, label='short', labelLines = F, sensPlot = T, maximumY =400)
SA6_plot <- plotCEA(res_SA6, cuttime=10, label='short', labelLines = F, sensPlot = T, maximumY =400)

#------------- Figure 2: sensitivity analysis plots -------------# 
plot_grid(SA1_plot + ggtitle("Time between FIGO 1 and FIGO 2/3 decreased to 4 years"), #A
          SA2_plot + ggtitle("CIN3 progression probability increased to 0.7"), #B
          SA3_plot + ggtitle("Included disutility for CIN0/1 diagnosis"), #C
          SA4_plot + ggtitle("ICER threshold increased to 50,000EUR"), #D
          SA5_plot + ggtitle("Only healthcare costs included"), #E
          SA6_plot + ggtitle("100% primary self-sampling assumed"),#F
          nrow=2, ncol=3, 
          labels= paste0("(", LETTERS[1:6], ")"))

ggsave("~/Projects/Colposcopy CEA/IJC Submission/Figure3.pdf", dpi = 700, width=20, height=11) # (2200x1200)


