### Script for testing spatial and temporal autocorrelation in best model residuals ###

## Load libraries
library(dplyr)
library(tidyverse)
library(rjags)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

## Load Best Model
# all model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]
# best model:
best <- 38  # also important: 42, 29, 7
best_model <- models[best] 
# load model:
load(paste0(dir, "model_runs/", best_model))

## Get Model Information and Data
# get parameter values:
best_model_params <- model_params[best,]
# model name:
model_name <- best_model
if (grepl("cov", name) == TRUE){
  model_name <- str_match(name, "cov_(.*?)_data")[,2]
} else {
  model_name <- str_match(name, "_(.*?)_model")[,2]
}
# load model inputs:
model_inputs <- model_info$metadata$model_data

## Extract Model Parameters:
# model output:
jags_out <- model_info$jags_out
# model parameter names:
jags_vars <- varnames(jags_out)
params <- jags_out[,grep("r|R|^tau|^b|^at", jags_vars)]
# remove burn in:
burn_in = 25000
params_burn <- window(params, start = burn_in)
# convert output to dataframe:
params_out <- as.data.frame(as.matrix(params_burn))
# separate out specific params:
x_params <- grep("^x", colnames(out))
beta_params <- grep("^b",colnames(out))
taus <- grep("tau", colnames(out))
r_int <- grep("r0", colnames(out))
r_rate <- grep("R", colnames(out))
# extract:
model_x <- out[,x_params]
x_ci <- apply(model_x, 2, quantile, c(0.05, 0.5, 0.95))
x_cols <- colnames(x_ci)
x_means <- matrix(NA, nrow = 5000, ncol = 7)
for (i in 1:nrow(x_means)){
  cols <- which(!is.na(str_match(x_cols, paste0(paste0("x\\[", i, ",")))))
  col_values <- x_ci[,cols]
}
## 