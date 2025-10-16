# GenotypingCEA
Code for analysis and results in paper "Cost-effectiveness analysis of human papillomavirus (HPV) genotyping strategies for management of HPV-positive women in cervical cancer screening" (under review)

- `codes/`: codes used in this analysis
- `plots/`: plots presented in the manuscript
- `results/`: output csv files

## codes
- `0_data.R`: data preparation file
- `1_functions.R`: file with functions for main analysis
- `2_parameters.R`: file with parameter inputs for main analysis
- `3_analysis.R`: main analysis file (base case analysis) and results (Figure 1, Table 3)
- `4_sensitivity_analysis.R`: sensitivity analysis and results (Figure 2)

## plots
- `Figure1_NMB.png`: Results under the base-case assumptions. For times since CIN2/3 from 0 to 10 years, panel (A) shows the net monetary benefit of nineteen HPV genotyping strategies and (B) shows the probabilistic analysis results (note: strategies not labelled had 0% probability of having the highest NMB across all time points). 
- `Figure2_NMB.png`: Results of one-way sensitivity analyses for the net monetary benefit (NMB) of nineteen strategies for times since onset CIN2/3 from 0 to 10 years.

## results
- `Table_3.csv`: Total number of referrals per 100,000 HPV-positive women in the baseline round, net monetary benefit, number of CIN2/3 and cancers detected per 100,000 HPV-positive women over two rounds of screening if the time since CIN2/3 onset per strategy at baseline is set at 5 and 10 years

