# CEA strategies
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
library(latex2exp)
library(cowplot)
library(gridExtra)

# Analysis preparation ----------------------------------------------------

load("data/ref.data.cens.RData")
source("1_functions.R")
source("2_parameters.R") # ensure parameters are set

set.seed(130601)
file.end <- format(Sys.Date(), "%Y%m%d")


# Create strategy dataframe -----------------------------------------------
# (dataframe where each row is one strategy under consideration by making combinations of HPV, ascus and normal triage options)
hpv.strategies <-  rep(c("allHPV+", "extended", "partial", "cyt&genotype"), times=c(1, 4, 6, 8))
ascus.strategies <- c("referAll", "referAll", "referAll", "repeat", "repeat", "referAll", "referAll", "refer7", "refer7", "repeat", 
                      "repeat", "referAll", "referAll", "refer7", "refer7", "refer1618", "refer1618", "repeat", "repeat")
normal.strategies <- c("referAll", rep(c("repeat", "repeat7"), 9))
strategies <- data.frame(hpv.strategies, ascus.strategies, normal.strategies)
onsetTimes <- seq(0, 10, 1)

# Basecase analysis (no bootstrap) ----------------------------------------
basecase_res <- do.call(rbind, lapply(onsetTimes, function(x) one_rep(ref.data.cens, x)))

# Figure 1a: NMB of each strategy at each time point compared to the reference
nmb_plot <- plotCEA(basecase_res, cuttime=10, maximumY = 400) 

# Table S2: number of referrals, cancers at 5 and 10 years per 100,000 woman in each strategy 
basecase_res %>% select(strategy, time, referred, NMB_d, cin23_2rounds, cxca_2rounds) %>% subset(time %in% c(5, 10)) %>% 
  pivot_wider(names_from = time, values_from = c(NMB_d, cin23_2rounds, cxca_2rounds)) %>% 
  mutate(across(starts_with("c"), ~ round(., 0))) %>%
  mutate(across(starts_with("N"), ~ round(., 1))) %>%  
  mutate(across(starts_with("r"), ~ round(., 0))) %>%  
  mutate(across(everything(), as.numeric)) %>% 
  select(strategy, referred, NMB_d_5, cin23_2rounds_5, cxca_2rounds_5, NMB_d_10, cin23_2rounds_10, cxca_2rounds_10) %>% 
  write.csv("results/Table3.csv", row.names = F)


# Bootstrap analysis (for 95% CI's & for probabilistic analysis) ----------
myCluster <- makeCluster(4) # type of cluster for parallel programming to speed up computation 
registerDoParallel(myCluster)
n.reps <- 500 
res_temp1 <- foreach (i = 1:n.reps) %dopar% {
  library(dplyr)
  load('data/ref.data.cens.RData') # data already organised/prepared for analysis
  source("code/1_functions.R") # functions for CEA and to get number of referrals
  source("code/2_parameters.R")# reset parameters
  
  dat <- ref.data.cens[sample(1:dim(ref.data.cens)[1], replace = T),] # for each rep resample POBASCAM data
  temp <- lapply(onsetTimes, function(x) one_rep(dat, x)) # loop through the time points and store results for each strategy at that time since onset
  temp
}
stopCluster(myCluster)  ; save(res_temp1, file='results/bootstrap_results_20251023.Rdata')

res_temp1 <- do.call(rbind, unlist(res_temp1, recursive=F))
res_temp1$rep <- rep(1:n.reps, each=length(onsetTimes)*length(ascus.strategies))

res_temp1 <- subset(res_temp1, time<=10)

# Values for results section 
cbind(res_temp1, NMB_bc = rep(basecase_res$NMB_d, n.reps)) %>% group_by(strategy, time) %>%
  summarise(nmb = paste0(spec_dec(NMB_bc[1], 1), " (", spec_dec(quantile(NMB_d, 0.025), 1), ", ", spec_dec(quantile(NMB_d, 0.975),1), ")")) %>%
  pivot_wider(names_from = time, values_from = nmb) %>% select(`0`, `10`)

# Probabilistic analysis --------------------------------------------------
# for each time since onset, find the probability a strategy comes out best in the bootstrap reps 
prob.df <- matrix(NA, nrow=19, ncol=length(onsetTimes))
for (t in onsetTimes){
  prob.best <- rep(0, 19)
  bests <- sapply(1:n.reps, function(x) which.max(res_temp1$NMB_d[res_temp1$time == t & res_temp1$rep == x]))
  prob.best <- tabulate(bests, nbins = 19)
  prob.df[, t + 1] <- (prob.best / n.reps) * 100  # Convert to percentage
}
prob.df <- cbind(1:19, prob.df) ; colnames(prob.df) <- c("Strategy",  0:10)

# create data for the plot
long <- prob.df %>% as.data.frame() %>% pivot_longer(cols = `0`:`10`, names_to = "time", values_to = "prob_best")
long$prob_best <- as.numeric(long$prob_best)
long$time <- as.numeric(long$time)
long$Strategy <- droplevels(factor(long$Strategy, levels=19:1))

long <- long %>% group_by(Strategy) %>% mutate(sum.strat =sum(prob_best)) %>% subset(sum.strat !=0) %>% as.data.frame() # remove strategies that are always 0%
long$best <- rep("grey", nrow(long))
extraColours <- rcartocolor::carto_pal(n=6, "Bold")
long$best[long$Strategy == 17] <- extraColours[1]
long$best[long$Strategy == 9] <- extraColours[2]
long$best[long$Strategy == 11] <- extraColours[3]

long$rightLabel <- case_when( # labels made specifically for these results
  long$Strategy == 17 ~ "17",
  long$Strategy == 11 ~ "11",
  long$Strategy == 9 ~ "9",
  long$Strategy == 19 ~ "19",
  long$Strategy == 15 ~ "15",
  long$Strategy == 16 ~ "16",
  long$Strategy == 10 ~ "10, 14, 8",

  T ~ ""
)

# Figure 1b: probabilistic analysis plot (using bootstrap samples)
prob_plot <- ggplot(long[long$time!=15,], aes(x=time, y=prob_best, color=best, group=Strategy)) + 
  theme_classic() + 
  geom_line(data=long[long$best =='grey' & long$time!=15,], lwd=0.5) +
  geom_line(data=long[long$best !='grey' & long$time!=15,], lwd=0.8) +
  ylab("Probability of a strategy having the highest NMB (%)") + 
  xlab("Time since CIN2/3 onset (years)") + 
  scale_y_continuous(expand=c(0,0), limits = c(0, 100), breaks=seq(0, 100, 20)) + 
  scale_x_continuous(expand=c(0,0), breaks=seq(0, 10, 1), limit=c(0, 10)) +
  theme(legend.position = 'none', plot.margin = unit(c(1,2,1,0.5), "cm"), axis.text = element_text(size=15), text = element_text(size=13)) + 
  coord_cartesian(clip = "off") +
  scale_color_manual(values = setNames(long$best, long$best)) +
  geom_text_repel(data=long[long$time == 10 & long$rightLabel!='', ], aes(color = best, label = rightLabel),
                  direction = "y", xlim = c(10.05, Inf), size=3, hjust = 0, segment.size = 0.1,
                  segment.alpha =0.1, max.overlaps = 100, box.padding = 0.05, show.legend = F,  segment.angle = 9)




# Combine plots for manuscript --------------------------------------------
nmb_plot <- plotCEA(basecase_res, cuttime=10, maximumY = 350) 

plot_grid(nmb_plot, prob_plot, ncol=2, labels=c("(A)", "(B)"))
ggsave("~/Projects/Colposcopy CEA/IJC Submission/Figure2.pdf", dpi=600, width=18, height=10) # save plot 


