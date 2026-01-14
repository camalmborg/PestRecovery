### Script for looking at the outputs of the forecasts ###

## Load libraries
library(dplyr)
library(readr)
library(stringr)
library(scoringRules)

## Load forecasts
# set working directory:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/"
setwd(dir)

## Load model forecast files:
files <- list.files(pattern = "result.csv$")
model_num = as.numeric(Sys.getenv("SGE_TASK_ID"))
# get years of analysis:
start_year <- as.numeric(stringr::str_extract(files[model_num], "(?<=start_year_)\\d+"))
years <- start_year:2023
# loading predicted forecast values file:
model_out <- read.csv(files[model_num])

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
  rename_with(~ str_replace_all(.x, c("^\\s*X" = "", "\\." = "-"))) %>%
  # select forecast years:
  select(matches(as.character(years))) %>%
  # rename columns for years:
  rename_with(~ as.character(years)[seq_along(.)])

# prepare for residual calculation:
pred <- model_out %>%
  # add baseline to predictions:
  mutate(across(!site, ~ tcg_base$baseline - .x)) %>%
  # ensemble means:
  mutate(site = as.factor(site)) %>%
  group_by(site) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  select(-c(site)) %>%
  # rename columns for years:
  rename_with(~ as.character(years)[seq_along(.)])
  

## Get RMSE, Mean Absolute Error (MAE) and bias for each model
# calculate model residuals:
resids <- matrix(nrow = nrow(tcg), ncol = ncol(tcg))
colnames(resids) <- as.character(years)
site = 1:nrow(tcg)
# run scores for every site/year:
for (s in site){
  for (i in years){
    # get observation:
    obs <- tcg[s, as.character(i)]
    # get ensemble for each site:
    ens <- as.numeric(pred[s, as.character(i)])
    # get the difference between model and observation:
    if (is.na(obs)|any(is.na(ens))){
      diff <- NA
    } else {
      diff <- ens - obs
    }
    # residuals:
    resids[s, as.character(i)] <- diff
  }
}

# make results table:
resid_calls <- c("RMSE", "MAE", "Bias")
resid_result <- matrix(nrow = length(resid_calls), ncol = ncol(tcg))
rownames(resid_result) <- resid_calls
colnames(resid_result) <- as.character(years)
# loop over years:
for (i in 1:ncol(resids)){
  resid_result[1, i] <- sqrt(mean(resids[,i]^2, na.rm = TRUE))
  resid_result[2, i] <- mean(abs(resids[,i]), na.rm = TRUE)
  resid_result[3, i] <- mean(resids[,i], na.rm = TRUE)
}

## Saving results
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/RMSE_Bias/"
# select information from model name to match file names from forecasts:
model <- stringr::str_extract(files[model_num], "(?<=model_)\\d+")
# is base?
if (grepl("base", files[model_num]) == TRUE){
  # add base to model name:
  model <- paste0("base_", model)
}
# file name and save:
write.csv(resids, paste0(save_dir, Sys.Date(), "_model_", model, "_start_year_", as.character(start_year), "_resids_raw_all_sites.csv"))
write.csv(resid_result, paste0(save_dir, Sys.Date(), "_model_", model, "_start_year_", as.character(start_year), "_rmse_mae_bias.csv"))


### Archive ###
# test <- sqrt(mean(resids[1,]^2))
# test <- mean(abs(resids[1,]))
# test <- mean(resids[1,])
