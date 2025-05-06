# This script is for checking out model results prior to making presentation
# and publication figures

# load libraries
librarian::shelf(tidyverse, dplyr, rjags, coda)

# navigate to folder:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Mortality_Model_Runs/"
setwd(dir)

## model_list:
models <- list.files(paste0(dir, "model_runs"))
uni_model_list <- models[grep("log", models)]
multi_model_list <- models[grep("covs", models)]

## loading results:
load("2025_05_06_univar_model_results_list.RData")    # univariate results
load("2025_05_01_multivar_model_results_list.RData")  # multivariate results

## loading dic results:
all_dics <- read.csv("2024_04_29_all_mort_model_dics.csv")[-1]

### Visualizations and parsing outputs:
# selecting model:
dir = dir
num = 1
model <- all_dics$model[all_dics$perform == num]
if(all_dics[num,]$type == "multi"){
  model_name <- multi_model_list[grep(paste0("modelrun_", model), multi_model_list)]
} else if (all_dics[num,]$type == "uni" & all_dics[num,]$mod_info == "log") {
  model_name <- uni_model_list[grep(paste0("modelrun_", model, "_log"), uni_model_list)]
} else if (all_dics[num,]$type == "uni" & all_dics[num,]$mod_info == "nolog") {
  model_name <- uni_model_list[grep(paste0("modelrun_", model, "_nolog"), uni_model_list)]
}

load(paste0(dir, "model_runs/", model_name))  # will load in environment as "model_info"

# visual inspections:
jags_out <- model_info$jags_out
vars <- varnames(jags_out)
params <- jags_out[,grep("^alpha|b|q|tau", vars)]

plot(params)
gelman.diag(params)
gelman.plot(params)
effectiveSize(params)
