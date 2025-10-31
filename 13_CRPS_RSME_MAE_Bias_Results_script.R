### CRPS and Bias Results Script ###

## Load libraries
library(dplyr)
library(tidyverse)
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

## Plotting CRPS
# x <- colnames(crps_all)[-1]
# y <- log10(crps_all[1, 2:ncol(crps_all)])
# plot(x, y, type = "l")
# for (i in 2:nrow(crps_all)){
#   lines(x, log10(crps_all[i, 2:ncol(crps_all)]))
# }

## RMSE, MAE, and Bias
# load files:
bias_load <- lapply(paste0(rmse_dir, rmse_files), read.csv)
# identify the maximum col num:
max_col <- max(sapply(bias_load, ncol))
# colnames: 
columns <- c("Model", "Metric", as.character(2017:2023))
# metric:
metric <- c("RMSE", "MAE", "Bias")

# how to make it into one dataframe:
rmse_all <- matrix(data = NA, nrow = length(bias_load), ncol = max_col + 1)
mae_all <- matrix(data = NA, nrow = length(bias_load), ncol = max_col + 1)
bias_all <- matrix(data = NA, nrow = length(bias_load), ncol = max_col + 1)
colnames(rmse_all) <- columns
colnames(mae_all) <- columns
colnames(bias_all) <- columns
model_names <- c()

# loading one:
#i = 1
# loop to fill in matrices:
for (i in 1:length(bias_load)){
  # select name:
  name <- rmse_files[i]
  name_extract <- str_extract(name, "(?<=_model_).*?(?=_rmse)")
  # add name to column
  model_names[i] <- name_extract
  # separate all bias data:
  all_bias <- bias_load[[i]] %>%
    rename(Metric = 1) %>%
    rename_with(~ str_replace_all(., c("X" = ""))) %>%
    pivot_longer(
      cols = starts_with("2"),
      names_to = "Year",
      values_to = "Score")
  # select each group:
  rmse <- all_bias %>%
    filter(Metric == "RMSE") %>%
    pivot_wider(names_from = Year,
                values_from = Score) %>%
    select(-Metric)
  mae <- all_bias %>%
    filter(Metric == "MAE") %>%
    pivot_wider(names_from = Year,
                values_from = Score) %>%
    select(-Metric)
  bias <- all_bias %>%
    filter(Metric == "Bias") %>%
    pivot_wider(names_from = Year,
                values_from = Score) %>%
    select(-Metric)
  # put it in the matrix:
  for (j in colnames(rmse)){
    rmse_all[i, j] <- suppressWarnings(as.numeric(rmse[,j]))
    mae_all[i, j] <- suppressWarnings(as.numeric(mae[,j]))
    bias_all[i, j] <- suppressWarnings(as.numeric(bias[,j]))
  }
}
rmse_all <- as.data.frame(rmse_all)
mae_all <- as.data.frame(mae_all)
bias_all <- as.data.frame(bias_all)
rmse_all$Model <- model_names
rmse_all$Metric <- "RMSE"
mae_all$Model <- model_names
mae_all$Metric <- "MAE"
bias_all$Model <- model_names
bias_all$Metric <- "Bias"

## Plotting RSME
# x <- colnames(rmse_all)[-c(1,2)]
# y <- log10(rmse_all[1, 3:ncol(rmse_all)])
# plot(x, y, type = "l")
# for (i in 2:nrow(rmse_all)){
#   lines(x, log10(rmse_all[i, 3:ncol(rmse_all)]))
# }

## Preparing data for score vs lead time and score vs 1-year lag plots

# function for getting results:
#'@param test = crps/rmse/mae/bias data
#'@param starts = character object with years being forecast 
get_plot_data <- function(test, starts){
  result <- test %>%
    # add column for model number for indexing:
    mutate(model_num = as.numeric(substr(result$Model, start = 1, stop = 1)), .after = 1)
  # make empty matrix to fill with diagonals:
  plots <- matrix(NA, nrow = length(starts)*3, ncol = length(starts))
  for (i in 1:nrow(plots)){
    model_num <- result$model_num[i]
    get_years <- result[i, grep("2", colnames(result))]
    nas <- which(is.na(get_years)) 
    if (length(nas) > 0){
      casts <- result[result$model_num == model_num, starts[-nas]]
    } else {
      casts <- result[result$model_num == model_num, starts]
    }
    diag <- diag(as.matrix(casts))
    plots[i,1:length(diag)] <- diag
  }
  colnames(plots) <- starts
  # get diagonal means and one-year lags:
  plot_data <- as.data.frame(plots) %>%
    mutate(model_num = result$model_num, .before = 1) %>%
    mutate(diag_mean = rowMeans(select(., -1), na.rm = TRUE))
  nums <- which(complete.cases(plot_data))
  plot_data$yr_one_lag <- c(plots[nums[1],], plots[nums[2],], plots[nums[3],])
  plot_data <- plot_data[,-grep("2", colnames(plot_data))]
}

# Run a test version:
# years of forecasts:
starts <- as.character(2018:2023)
test <- crps_all
crps_plot <- get_plot_data(test, starts)
