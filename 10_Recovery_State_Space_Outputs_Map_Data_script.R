### Recovery Rate Across Sites script

## Load libraries and necessary environments
librarian::shelf(dplyr, tidyverse, rjags, coda)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

## pull in model output files if not in the environment already
model_params <- read.csv("2025_07_31_all_base_uni_recov_models_param_means.csv")
load("2025_07_31_recov_models_outputs_list.RData")  # object is called model_outputs
# model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]

# choose model:
m_num <- 1  # change this when changing models
model_pick <- models[m_num]
# load model_info object
load(paste0(dir, "model_runs/", model_pick))
model_inputs <- model_info$metadata$model_data
# model parameters:
out <- as.matrix(model_info$jags_out)
# separate out specific params:
x_params <- grep("^x", colnames(out))
beta_params <- grep("^b",colnames(out))
taus <- grep("tau", colnames(out))
r <- grep("r0", colnames(out))