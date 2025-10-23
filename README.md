# GenotypingCEA
Code for analysis and results in paper "Cost-effectiveness analysis of human papillomavirus (HPV) genotyping strategies for management of HPV-positive women in cervical cancer screening" (under review)

- `codes/`: codes used in this analysis
- `data/`: data used in this analysis
  
## codes
- `1_functions.R`: file with functions for main analysis
- `2_parameters.R`: file with parameter inputs for main analysis
- `3_analysis.R`: main analysis file (base case analysis) and results (Figure 1, Table 3)
- `4_sensitivity_analysis.R`: sensitivity analysis and results (Figure 2)

## data
- `NCR-export-2025_04_07_01_43_21.xlsx`: data on mortality age distribution for cervical cancer from the NCR database 
- `NCR-export-29_01_2025_11_06_33.xlsx`: data on the age/cancer stages distribution from the NCR database
- `NKI_data.RData`: 10-year survival data on cervical cancer by cancer stage from the NCR database
- `relative_survivals.RData`: 10-year survival data, column 1 is from the NCR database and column 2 accounts for stage shift if the cancer was detected 5 years later 
