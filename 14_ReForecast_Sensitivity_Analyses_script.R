### Script for sensitivities to drivers analyses ###

## Load libraries
library(dplyr)
library(tidyverse)
library(rjags)
library(coda)
library(ggplot2)
library(tigris)

## Loading TCG and TCG baselines:
# get baselines:
tcg_base <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  # rename:
  rename_with(~ str_replace_all(.x, c("^\\s*X" = "", "\\." = "-"))) %>%
  # get baseline for anomolies:
  mutate(baseline = rowMeans(select(., `2010-05-01`:`2015-05-01`), na.rm = TRUE), .before = 1) %>%
  # create anomalies from baseline:
  mutate(across(!baseline, ~ baseline - .x))

tcg <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  # rename:
  rename_with(~ str_replace_all(.x, c("^\\s*X" = "", "\\." = "-"))) #%>%


## Loading parameter values from best models:
# load model performance information:
dic_sort <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/2025_11_30_all_recov_models_dics.csv", row.names = 1)

## Load Best Model
# all model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]
# best model:
best <- dic_sort$model_number[1]
best_model <- models[best] 
# load model:
load(paste0(dir, "model_runs/", best_model))
# load model_params:
model_params <- read.csv(file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/2025_11_30_all_recov_models_param_means.csv")

# collect site random effect model params:
jags_out <- model_info$jags_out
vars <- varnames(jags_out)
site_params <- jags_out[,grep("^as", vars)]
# remove burn in:
burn_in = 25000
site_params_burn <- window(site_params, start = burn_in)
out_site <- as.matrix(site_params_burn)
space_re <- apply(out_site, 2, mean, na.rm = T)
# quickly attach to lat/lon for making a map later:


## Get Model Information and Data:
# get parameter values:
best_params <- model_params[best,]
# separate parameters for calling:
time_re <- best_params[grep("atime", names(best_params))]

