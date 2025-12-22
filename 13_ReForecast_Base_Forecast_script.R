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

## Input data
# get mean tcg decay baseline 2017 data:
tcg_base <- tcg_ts %>%
  # pre-disturbance 5 year baseline:
  mutate(baseline = rowMeans(select(., `2010-05-01`:`2015-05-01`), na.rm = TRUE), .before = 1) %>%
  # pre-disturbance 5 year sd:
  mutate(sd = apply(pick(`2010-05-01`:`2015-05-01`), 1, sd, na.rm = TRUE), .before = 2)
# extract them as vectors:
means <- tcg_base$baseline
sds <- tcg_base$sd

# forecast for base (suggests next year is just mean of last year plus some sd as prec):
run_base_forecast <- function(start, end, ns, n_ens, means, sds, yr){
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
      # loop over time:
      for (t in 1:ncol(N)){
        N[i,t] <- N[i,t] <- rnorm(1, mean = means[s], sd = sds[s])
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
# data.frame of model jobs:
model_jobs <- data.frame(task_id = 1:6, 
                         year_run = rep(c(1:6)))

# setting task id for cluster runs:
task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# preparing inputs:
# number of sites:
ns = 5000
# number of ensemble members:
n_ens = 1500
# years of analysis:
years <- 2017:2023
n_yr <- 1:(length(years)-1)
# use taskid to select proper years for model run:
y = model_jobs$year_run[task_id]

## Run base forecast
base_forecast <- run_base_forecast(start = years[y], 
                                   end = last(years), 
                                   ns = ns, 
                                   n_ens = n_ens, 
                                   means = means,
                                   sds = sds,
                                   yr = y)
write.csv(base_forecast, file = paste0("Recovery_Forecasts/", Sys.Date(),
                                    "_ens_", as.character(n_ens),
                                    "_base_model_", as.character(task_id),
                                    "_start_year_", as.character(years[y]),
                                    "_reforecast_result.csv"),
          row.names = FALSE)

