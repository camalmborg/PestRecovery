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
# getting means, low CI (0.025 quantile), high CI (0.975 quantile):
model_means <- matrix(NA, nrow = length(model_outputs), ncol = 6)
model_lows <- matrix(NA, nrow = length(model_outputs), ncol = 6)
model_highs <- matrix(NA, nrow = length(model_outputs), ncol = 6)
for (i in 1:length(model_outputs)){
  # add model name:
  model_means[i,1] <- names(model_outputs)[i]
  model_lows[i,1] <- names(model_outputs)[i]
  model_highs[i,1] <- names(model_outputs)[i]
  # get results from list:
  result <- model_outputs[[i]]
  if (TRUE %in% grepl("beta", colnames(result))){
    # separate beta cols
    beta_means <- apply(result[grep("beta", colnames(result))], 2, mean)
    beta_lows <- apply(result[grep("beta", colnames(result))], 2, quantile, c(0.025))
    beta_highs <- apply(result[grep("beta", colnames(result))], 2, quantile, c(0.975))
    for (j in 1:length(beta_means)){
      model_means[i,j+1] <- round(beta_means[j], 3)
      model_lows[i,j+1] <- round(beta_lows[j], 3)
      model_highs[i,j+1] <- round(beta_highs[j], 3)
    }
  }
}
# remove unnecessary things:
rm(result)

# making one data frame with all the results compiled:
model_results_CIs <- as.data.frame(model_means) %>%
  rename(model = 1) %>%
  rename_with(~ paste0("mean_", seq_along(.)), .cols = -c(model)) %>%
  # add the lower CI
  bind_cols(model_lows[,2:ncol(model_lows)]) %>%
  rename_with(~ paste0("low_", seq_along(.)), .cols = -c(1:6)) %>%
  # add the higher CI
  bind_cols(model_highs[,2:ncol(model_highs)]) %>%
  rename_with(~ paste0("high_", seq_along(.)), .cols = -c(1:11)) #%>%
  # # sort long
  # pivot_longer(
  #   cols = starts_with(c("mean", "high", "low")),
  #   names_to = c("CI", "number"),
  #   names_sep = c("_"),
  #   values_to = "value")

