### Re-Forecast script for testing recovery rate state space models ###

## Load libraries
library(dplyr)
library(rjags)
library(coda)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)
# load model performance information:
dic_sort <- read.csv("2025_10_06_all_recov_models_dics.csv", row.names = 1)

## Get model
# model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]
# setting task id for cluster runs:
model_num <- as.numeric(Sys.getenv("SGE_TASK_ID"))
#model_num = 2
# get model:
top_model <- dic_sort[dic_sort$perform == model_num,]
top_model <- top_model$model_number
# load the best model:
load(paste0(dir, "model_runs/", models[as.numeric(top_model)]))


## Prepare the forecast model covariates
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

## Get model parameters for running forecast
#'@param model_info = model output that loads from folders including model posterior and input data saved in model runs
#'@param n_ens = number of ensemble members for running re-forecast
get_params <- function(model_info, n_ens){
  ## Sample parameters from the model posterior
  # get parameters from model output:
  posterior <- as.matrix(model_info$jags_out)
  # random sample of an entire row:
  sample <- sample(1:nrow(posterior), n_ens, replace = FALSE)
  # make them the right objects to go into the forecast:
  # separate out specific params:
  beta_params <- grep("^b",colnames(posterior))
  taus <- grep("tau", colnames(posterior))
  r <- grep("r0", colnames(posterior))
  tau_time <- grep("at", colnames(posterior))
  x_ic <- grep("^x\\[[0-9+,1\\]", colnames(posterior))
  # group them:
  params <- cbind(posterior[,beta_params], posterior[,taus], posterior[,tau_time], r = posterior[,r], posterior[,x_ic])
  # sample:
  params <- cbind(col_num = sample, params[sample,])
  return(params)
}

## Forecast functions to run across all sites and ensemble members:
#'@param start = numeric, year prior to recovery initiation: 2017
#'@param end = numeric, year end of forecast: 2023
#'@param ns = numeric, number of sites: 5000
#'@param n_ens = numeric, number of ensemble members
#'@param params = matrix/df, sampled parameters from get_params function
#'@param yr = numeric, year starting place of forecast: e.g. 2017 = 1, 2018 = 2, etc

# for models with 3 variables (2 time-varying, 1 static) - e.g. models 1, 3 (also 5)
run_forecast_3_var <- function(start, end, ns, n_ens, params, yr){
  # time steps:
  nt = 1:(end-start)
  # empty list to hold result:
  forecast_result <- list()
  # loop over sites:
  for (s in 1:ns){
    # matrix to hold results:
    N <- matrix(NA, nrow = n_ens, ncol = length(nt)+1)
    # loop over ensemble members:
    for (i in 1:n_ens){
      ## Prep params and inputs
      # model parameters:
      # betas:
      betas <- params[n_ens, c(grep("beta", colnames(params)))]
      # taus - convert to SD from 1/precision:
      tau_obs <- 1/params[n_ens, c(grep("obs", colnames(params)))]
      tau_add <- 1/params[n_ens, c(grep("add", colnames(params)))]
      a_time <- c(0, 1/params[n_ens, c(grep("atime", colnames(params)))])
      a_time[is.infinite(a_time)] <- 0  # fix Inf value from conversion to SD from precisions
      # r0:
      r <- params[n_ens, c(grep("r", colnames(params)))]
      
      ## Run the forecast
      c_one = cov_one[s,]
      c_two = cov_two[s,]
      c_three = cov_three[s]
      # x_ic:
      x_ic <- grep(paste0("^x\\[[0-9+,", as.character(yr), "\\]"), colnames(params))
      N[,1] <- params[n_ens, x_ic][s]
      
      # loop over time:
      for (t in 2:ncol(N)){
        big_R <- r + a_time[t-1] + (betas[1]*c_one[,t-1]) + (betas[2]*c_two[,t-1]) + betas[3]*c_three
        N[i,t] <- rnorm(1, mean = N[i,t-1]*big_R, sd = tau_add)
      }
    }
    forecast_result[[s]] <- as.data.frame(N)
    print(s)
  }
  # bind together:
  forecast_result <- bind_rows(forecast_result, .id = "site")
  return(forecast_result)
}

## For models with 4 variables (3 time-varying, one static) - e.g. model 2
run_forecast_4_var <- function(start, end, ns, n_ens, params, yr){
  # time steps:
  nt = 1:(end-start)
  # empty list to hold result:
  forecast_result <- list()
  # loop over sites:
  for (s in 1:ns){
    # matrix to hold results:
    N <- matrix(NA, nrow = n_ens, ncol = length(nt)+1)
    # loop over ensemble members:
    for (i in 1:n_ens){
      ## Prep params and inputs
      # model parameters:
      # betas:
      betas <- params[n_ens, c(grep("beta", colnames(params)))]
      # taus - convert to SD from 1/precision:
      tau_obs <- 1/params[n_ens, c(grep("obs", colnames(params)))]
      tau_add <- 1/params[n_ens, c(grep("add", colnames(params)))]
      a_time <- c(0, 1/params[n_ens, c(grep("atime", colnames(params)))])
      a_time[is.infinite(a_time)] <- 0  # fix Inf value from conversion to SD from precisions
      # r0:
      r <- params[n_ens, c(grep("r", colnames(params)))]
      
      ## Run the forecast
      c_one = cov_one[s,]
      c_two = cov_two[s,]
      c_three = cov_three[s,]
      c_four = cov_four[s]
      # x_ic:
      x_ic <- grep(paste0("^x\\[[0-9+,", as.character(yr), "\\]"), colnames(params))
      N[,1] <- params[n_ens, x_ic][s]
      
      # loop over time:
      for (t in 2:ncol(N)){
        big_R <- r + a_time[t-1] + (betas[1]*c_one[,t-1]) + (betas[2]*c_two[,t-1]) + (betas[3]*c_three[,t-1]) + betas[4]*c_four
        N[i,t] <- rnorm(1, mean = N[i,t-1]*big_R, sd = tau_add)
      }
    }
    forecast_result[[s]] <- as.data.frame(N)
    print(s)
  }
  # bind together:
  forecast_result <- bind_rows(forecast_result, .id = "site")
  return(forecast_result)
}

## Running for different models
# preparing inputs:
# number of sites:
ns = 5000
# number of ensemble members:
n_ens = 2500

# get sampled parameters: nrow = n_ens
ens_params <- get_params(model_info, n_ens)
# years of analysis:
years <- 2017:2023
n_yr <- 1:(length(years)-1)

# making list for results from every year being "re-forecast":
if (model_num != 2){
  reforecast_list <- list()
  for (i in n_yr){
    reforecast_list[[i]] <- run_forecast_3_var(start = years[i], end = last(years), ns = ns, n_ens = n_ens, params = ens_params, yr = i)
  }
  save(reforecast_list, file = paste0("Recovery_Forecasts/", Sys.Date(), "_model_", as.character(model_num), "_reforecast_result.RData"))
} else if (model_num == 2){
  reforecast_list <- list()
  for (i in n_yr){
    reforecast_list[[i]] <- run_forecast_4_var(start = years[i], end = last(years), ns = ns, n_ens = n_ens, params = ens_params, yr = i)
  }
  save(reforecast_list, file = paste0("Recovery_Forecasts/", Sys.Date(), "_model_", as.character(model_num), "_reforecast_result.RData"))
}


#test <- run_forecast_4_var(start = 2017, end = 2023, ns = ns, n_ens = n_ens, params = ens_params, yr = 1)
#write.csv(test, file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/2025_10_16_re_forecast_test.csv")