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
# getting means, low CI (0.05 quantile, 5%), high CI (0.95 quantile, 95%)
model_means <- matrix(NA, nrow = length(model_outputs), ncol = 6)
model_lows <- matrix(NA, nrow = length(model_outputs), ncol = 6)
model_highs <- matrix(NA, nrow = length(model_outputs), ncol = 6)
for (i in 1:length(model_outputs)){
  # add model name:
  model_CIs[i,1] <- names(model_outputs)[i]
  # get results from list:
  result <- model_outputs[[i]]
  if (TRUE %in% grepl("beta", colnames(result))){
    # separate beta cols
    beta_means <- apply(result[grep("beta", colnames(result))], 2, mean)
    beta_lows <- apply(result[grep("beta", colnames(result))], 2, quantile, c(0.05))
    beta_highs <- apply(result[grep("beta", colnames(result))], 2, quantile, c(0.95))
    for (j in 1:length(beta_means)){
      model_means[i,j+2] <- round(beta_means[j], 3)
      model_lows[i,j+2] <- round(beta_lows[j], 3)
      model_highs[i,j+2] <- round(beta_highs[j], 3)
    } 
  } else {
    model_means[i] <- NA
    model_lows[i] <- NA
    model_highs[i] <- NA
  }
}


