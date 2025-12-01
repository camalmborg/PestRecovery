### Re-Forecast script for testing recovery rate state space models ###

## Load libraries
library(dplyr)
library(rjags)
library(coda)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

# load environment if needed:
load("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Environments/2025_07_07_environment.RData")

## Input data:
# first x:
x_miss <- mean(tcg[,grep("^2017",names(tcg))], na.rm = T)
x_prec_miss <- sd(tcg[,grep("^2017", names(tcg))], na.rm = T)
x1 <- tcg[,grep("2017", names(tcg))]
x1[which(is.na(x1))] <- x_miss
# x1 precision:
x_prec <- 1/(sd_tcg[,grep("2017", names(sd_tcg))]^2)
x_prec[which(is.na(x_prec))] <- 1/x_prec_miss
x_prec[which(x_prec == Inf)] <- 1/x_prec_miss

# forecast for base (suggests next year is just mean of last year plus some sd):
run_base_forecast <- function(start, end, ns, n_ens, xs, xprec, yr){
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
      # x_ic:
      x_ic <- xs[s]
      N[,1] <- x_ic
      
      # tau_add:
      tau_add <- sqrt(1/xprec[s,])
      
      # loop over time:
      for (t in 2:ncol(N)){
        N[i,t] <- rnorm(1, mean = N[i,t-1], sd = tau_add[s,t])
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
n_ens = 1500

# get sampled parameters: 
ens_params <- model_params[[model_num]]  # for these runs, using pre-sampled to make sure all time horizons have the same parameters for each model
#ens_params <- read.csv(paste0("Recovery_Forecasts/model_", as.character(model_num), "_params.csv"))
#ens_params <- get_params(model_info, n_ens)
# years of analysis:
years <- 2017:2023
n_yr <- 1:(length(years)-1)
# use taskid to select proper years for model run:
y = model_jobs$year_run[task_id]
