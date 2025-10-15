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
# get model:
top_model <- dic_sort[dic_sort$perform == model_num,]
top_model <- top_model$model_number
# load the best model:
load(paste0(dir, "model_runs/", models[as.numeric(top_model)]))

## Get model parameters for running forecast
get_params <- function(model_info, yr){
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
  x_ic <- grep(paste0("^x\\[[0-9+,", as.character(yr),"\\]"), colnames(posterior))
  # group them:
  params <- cbind(posterior[,beta_params], posterior[,taus], posterior[,tau_time], r = posterior[,r], posterior[,x_ic])
  # sample:
  params <- params[sample,]
  return(params)
}


## Prepare the forecast
# model covariates inputs:
covs <- grep("cov", names(model_info$metadata$model_data))
covariates <- do.call(cbind, model_info$metadata$model_data[covs])
ncovs = length(covs)
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

# time steps:
yr = 2
start = 2018
end = 2023
nt = 1:(end-start)
# number of ensemble members:
n_ens = 2
# number of sites:
ns = 2
# # matrix to hold results:
# N <- matrix(NA, nrow = n_ens, ncol = length(nt)+1)


## Run the forecast
# function to run across all sites and ensemble members:
run_forecast <- function(nt, ns, n_ens){
  # empty list to hold result:
  forecast_result <- list()
  # loop over sites:
  for (s in 1:ns){
    # matrix to hold results:
    N <- matrix(NA, nrow = n_ens, ncol = length(nt)+1)
    # loop over ensemble members:
    for (i in 1:n_ens){
      # model params:
      params <- get_params(model_info, yr)
      
      ## Prep params and inputs
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
      c_one = cov_one[s,]
      c_two = cov_two[s,]
      c_three = cov_three[s]
      # x_ic:
      N[,1] <- params[grep("^x", names(params))][s]
      
      # loop over time:
      for (t in 2:ncol(N)){
        R <- r + a_time[t-1] + (betas[1]*c_one[,t-1]) + (betas[2]*c_two[,t-1]) + betas[3]*c_three
        N[i,t] <- rnorm(1, mean = N[i,t-1]*R, sd = tau_add)
      }
    }
    forecast_result[[s]] <- as.data.frame(N)
    print(s)
  }
  # bind together:
  forecast_result <- bind_rows(forecast_result, .id = "site")
  return(forecast_result)
}

test <- run_forecast(nt = nt, ns = ns, n_ens = n_ens)
