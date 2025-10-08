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
model_num = 1

## Get model parameters for running forecast
get_params <- function(model_num){
  # get model:
  top_model <- dic_sort[dic_sort$perform == model_num,]
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
  tau_time <- grep("at", colnames(posterior))
  # group them:
  params <- cbind(posterior[,beta_params], posterior[,taus], posterior[,tau_time], r = posterior[,r])
  # sample:
  params <- params[sample,]
  return(params)
}
# model params:
params <- get_params(model_num)

## Prepare the forecast
# time steps:
start = 2017
end = 2023
nt = 1:(end-start)
# model covariates inputs:
covs <- grep("cov", names(model_info$metadata$model_data))
covariates <- do.call(cbind, model_info$metadata$model_data[covs])
ncovs = length(covs)
# initial conditions - 2017 tcg:
ic <- model_info$metadata$model_data$x_ic


## Prep params and inputs
# covariates:
cov_groups <- unique(sub("\\..*", "", colnames(covariates)))
cov_list <- list()
for (i in cov_groups){
  cov_list[[i]] <- covariates[,grep(i, colnames(covariates))]
}
# assign:
for (name in names(cov_list)) {
  assign(name, cov_list[[name]])
}
# model parameters:
# betas:
betas <- params[grep("beta", names(params))]
# taus - convert to SD from 1/precision:
tau_obs <- 1/params[grep("obs", names(params))]
tau_add <- 1/params[grep("add", names(params))]
a_time <- c(0, 1/params[grep("atime", names(params))])
a_time[is.infinite(a_time)] <- 0  # fix Inf value from conversion to SD from precisions
# r0:
r <- params[grep("r", names(params))]

## Run the forecast
# for single site
s = 1
ic = ic[s]
cov_one = cov_one[s,]
cov_two = cov_two[s,]
cov_three = cov_three[s]

# matrix to hold results:
N <- matrix(NA, nrow = 1, ncol = length(nt)+1)
N[,1] <- ic
# loop over time:
for (t in 2:length(N)){
  R <- r + a_time[t-1] + (betas[1]*cov_one[,t-1]) + (betas[2]*cov_two[,t-1]) + betas[3]*cov_three
  N[,t] <- rnorm(1, mean = N[,t-1]*R, sd = tau_add)
}



## Forecast function
#'@param nt = vector - timesteps
#'@param n_ens = numeric - number of ensemble members
#'@param ic = vector - initial condition, in this case 2017 TCG
#'@param params = vector - parameters from get_params function
#'@param covariates = covariates for model
# rss_forecast <- function(nt, n_ens, ic, params, covariates){
#   ## Prep params and inputs
#   # covariates:
#   cov_groups <- unique(sub("\\..*", "", colnames(covariates)))
#   cov_list <- list()
#   for (i in cov_groups){
#     cov_list[[i]] <- covariates[,grep(i, colnames(covariates))]
#   }
#   # assign:
#   for (name in names(cov_list)) {
#     assign(name, cov_list[[name]])
#   }
#   # model parameters:
#   # betas:
#   betas <- params[grep("beta", names(params))]
#   # taus - convert to SD from 1/precision:
#   tau_obs <- 1/params[grep("obs", names(params))]
#   tau_add <- 1/params[grep("add", names(params))]
#   tau_time <- 1/params[grep("atime", names(params))]
#   tau_time[is.infinite(tau_time)] <- 0  # fix Inf value from conversion to SD from precisions
#   # r0:
#   r <- params[grep("r", names(params))]
#   
#   ## Run the forecast
#   # matrix to hold results:
#   N <- matrix(NA, nrow = 1, ncol = last(nt))
#   # loop over time:
#   for (t in nt){
#     R[,t] <- r + tau_time[t] + (betas[1]*cov_one[,t]) + (betas[2]*cov_two[,t]) + betas[3]*cov_three[s]
#     N[,t] <- rnorm(N[,t-1] * R[,t], tau_add)
#     # x[s,t] <- dnorm(mu[s,t], tau_add)
#     # x[s,1] <- ic[s]
#     }
# }
# 
# 
