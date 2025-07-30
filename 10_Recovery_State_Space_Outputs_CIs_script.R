### Identifying significant variables with confidence interval overlap
# Script for working with slope (beta) confidence intervals

## Load libraries
library(dplyr)
library(tidyverse)
library(rjags)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

## Load model output files
model_params <- read.csv("2025_07_29_all_base_uni_recov_models_param_means.csv")
load("2025_07_29_recov_models_outputs_list.RData")  # object is called model_outputs

## Calculating CIs for all models
model_CIs <- matrix(NA, nrow = length(model_outputs), ncol = 3)
colnames(model_CIs) <- c("model", "")
for (i in 1:length(model_outputs)){
  
}
