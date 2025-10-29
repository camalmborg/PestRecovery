### CRPS and Bias Results Script ###

## Load libraries
library(dplyr)
library(stringr)

## Loading CRPS and Bias run results
# working directories:
crps_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/CRPS/"
rmse_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/RMSE_Bias/"
# files:
crps_files <- list.files(crps_dir)[grep("across_site", list.files(crps_dir))]
rmse_files <- list.files(rmse_dir)[grep("rmse_mae_bias", list.files(rmse_dir))]

## CRPS
# loading all as list:
crps_load <- lapply(paste0(crps_dir, crps_files), read.csv)
# identify the maximum col num:
max_col <- max(sapply(crps_load, ncol))
# colnames: 
columns <- c("Model", as.character(2017:2023))

# how to make it into one dataframe:
crps_all <- matrix(data = NA, nrow = length(crps_load), ncol = max_col)
model_names <- c()
colnames(crps_all) <- columns
# loop to fill in matrix:
for (i in 1:length(crps_load)){
  # select name:
  name <- crps_files[i]
  name_extract <- str_extract(name, "(?<=_model_).*?(?=_crps)")
  # add name to column
  model_names[i] <- name_extract
  # load member:
  crps <- crps_load[[i]] %>%
    rename(model = 1) %>%
    select(-model) %>%
    rename_with(~ str_replace_all(., c("X" = "")))
  # put it in the matrix:
  for (j in colnames(crps)){
    crps_all[i, j] <- as.numeric(crps[,j])
  }
}
crps_all <- as.data.frame(crps_all)
crps_all$Model <- model_names

## Plotting:


