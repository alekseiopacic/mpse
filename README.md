# Replication materials

This repository contains code to replicate "Monotonic Path-Specific Effects: Application to Estimating Educational Returns" (https://alekseiopacic.github.io/papers/JASA_ed_decomp_sm_merged.pdf). Download this repository by clicking on the green "Clone or download" button above. This code replicates the main figures and tables in the manuscript. Please contact aopacic@g.harvard.edu if you have any questions.

## Instructions:
To reproduce table 1 and figure 1, run each of the following files in turn: 01_dml.R, 02_rwr.R, and 03_out.R.

Data files provided are included in the Samples folder. Intermediate output created by 01_dml.R and 02_rwr.R will be saved in the main folder.

Figures and Tables created by R scripts are stored in the Output folder. 

## Software:
To run the R files you will need a current version of R. I used version 2023.06.0+421.

## R packages required:

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

