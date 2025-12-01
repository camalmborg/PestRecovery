### Re-Forecast params prep script for testing recovery rate state space models ###

## Load libraries
library(dplyr)
library(rjags)
library(coda)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)
# load model performance information:
dic_sort <- read.csv("2025_11_30_all_recov_models_dics.csv", row.names = 1)

## Get model
# model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]


## Get model parameters for running forecasts - used to save and keep all sample numbers the same across time horizons
#'@param model_num = numeric, model to upload based on dic performance
#'@param n_ens = number of ensemble members for running re-forecast
get_ens_params <- function(model_num, n_ens){
  # get model:
  top_model <- dic_sort[dic_sort$perform == model_num,]
  top_model <- top_model$model_number
  # load the best model:
  load(paste0(dir, "model_runs/", models[as.numeric(top_model)]))
  ## Sample parameters from the model posterior
  # remove burn in:
  jags_out <- model_info$jags_out
  burn_in = 50000
  jags_burn <- window(jags_out, start = burn_in) 
  # get parameters from model output:
  posterior <- as.matrix(jags_burn)
  # random sample of an entire row:
  sample <- sample(1:nrow(posterior), n_ens, replace = FALSE)
  # make them the right objects to go into the forecast:
  # separate out specific params:
  beta_params <- grep("^b",colnames(posterior))
  taus <- grep("tau", colnames(posterior))
  r <- grep("r0", colnames(posterior))
  tau_time <- grep("at", colnames(posterior))
  tau_site <- grep("as", colnames(posterior))
  x_ic <- grep("^x", colnames(posterior))
  # group them:
  params <- cbind(posterior[,beta_params], posterior[,taus], posterior[,tau_time], posterior[,tau_site], r = posterior[,r], posterior[,x_ic])
  # sample:
  params <- cbind(col_num = sample, params[sample,])
  return(params)
}


model_1_params <- get_ens_params(1, 1500)
#write.csv(model_1_params, file = "Recovery_Forecasts/model_1_params.csv", row.names = FALSE)

model_2_params <- get_ens_params(2, 1500)
#write.csv(model_2_params, file = "Recovery_Forecasts/model_2_params.csv", row.names = FALSE)

model_3_params <- get_ens_params(3, 1500)
#write.csv(model_3_params, file = "Recovery_Forecasts/model_3_params.csv", row.names = FALSE)

model_params <- list()
model_params[[1]] <- model_1_params
model_params[[2]] <- model_2_params
model_params[[3]] <- model_3_params
save(model_params, file = "Recovery_Forecasts/model_params_list.RData")

