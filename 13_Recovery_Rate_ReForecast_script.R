### Re-Forecast script for testing recovery rate state space models ###

## Load libraries
library(dplyr)
library(rjags)
library(coda)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

# model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]
# get model:
top_model <- dic_sort[dic_sort$perform == 1,]
top_model <- top_model$model_number
# load the best model:
load(paste0(dir, "model_runs/", models[as.numeric(top_model)]))

## Sample parameters from the model posterior
# get parameters from model output:
posterior <- as.matrix(model_info$jags_out)
# random sample of an entire row:
sample <- sample(1:nrow(posterior), 1)
# make them the right objects to go into the forecast:
# separate out specific params:
beta_params <- grep("^b",colnames(posterior))
taus <- grep("tau", colnames(posterior))
r <- grep("r0", colnames(posterior))
# group them:
params <- cbind(posterior[,beta_params], posterior[,taus], posterior[,r])

## Prepare the forecast
# time steps:
start = 2017
end = 2023
nt = 1:(end-start)

