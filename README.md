#######################
REPLICATION MATERIALS
#######################

AUTHOR: Aleksei Opacic

PAPER: "Monotonic Path-Specific Effects: Application to Estimating Educational Returns" (https://alekseiopacic.github.io/papers/JASA_ed_decomp_sm_merged.pdf)

UPDATED: 6 August 2024

CONTACT: aopacic@g.harvard.edu

INSTRUCTIONS:

Please unzip all files and put them in the same folder. 

To reproduce table 1 and figure 1, open ed_decomp.Rproj and run each file in turn: 01_dml.R, 02_rwr.R, and 03_out.R.

Data files provided are included in the Samples folder. Intermediate output created by 01_dml.R and 02_rwr.R will be saved in the main folder.

Figures and Tables created by R scripts are stored in the Output folder. 

SOFTWARE:

To run the R files you will need a current version of R. I used version 2023.06.0+421.

R PACKAGES REQUIRED:

- tidyverse
- survey
- rlang
- caret
- glmnet
- ranger
- SuperLearner
- parallel
- latex2exp
- scales
- xtable

