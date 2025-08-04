### Sample Time Series plots ###

## Load libraries and necessary environments
librarian::shelf(dplyr, tidyverse, rjags, coda)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

## pull in model output files if not in the environment already
model_params <- read.csv("2025_07_31_all_base_uni_recov_models_param_means.csv")
load("2025_07_31_recov_models_outputs_list.RData")  # object is called model_outputs

## Make sample time series for best performing models:
