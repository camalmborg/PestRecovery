### Script for looking at the outputs of the forecasts ###

## Load libraries
library(dplyr)
library(readr)
library(scoringRules)

## Load forecasts
# set working directory:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/"
setwd(dir)
# load observation data:
years <- 2017:2023
tcg <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  # select columns with observations for 2017-2023:
  select(matches(as.character(years))) 
  
# load model forecast files:
files <- list.files(pattern = "result.csv$")
#forecast_list <- lapply(files, read.delim)

# loading test file:
test <- read.csv(files[1])

testing <- test %>%
  # rename columns with years:
  rename_with(~ as.character(years)[seq_along(.)], .cols = -1) %>%
  # remove 2017 to get only the predicted values:
  select(-"2017")

  
  
# testing <- test %>%
#   group_by(site) %>% 
#  # summarise_all(., mean, na.rm = TRUE) %>%
#   summarise_all(list(min, max, mean, na.rm = TRUE))

# means <- apply(test, 2, mean, na.rm = TRUE)
# lower <- apply(test, 2, quantile, probs = c(0.05), na.rm = TRUE)
# upper <- apply(test, 2, quantile, probs = c(0.95), na.rm = TRUE)
