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
  select(matches(as.character(years))) %>%
  # rename columns for years:
  rename_with(~ as.character(years)[seq_along(.)])
  
# load model forecast files:
files <- list.files(pattern = "result.csv$")
model_num = as.numeric(Sys.getenv("SGE_TASK_ID"))
# loading test file:
model_out <- read.csv(files[model_num])
pred <- model_out %>%
  # rename columns with years:
  rename_with(~ as.character(years)[seq_along(.)], .cols = -1)

## Get CRPS scores for ensemble model
# make matrix to hold scores:
crps_scores <- matrix(nrow = nrow(tcg), ncol = ncol(tcg))
colnames(crps_scores) <- as.character(years)
site = 1:nrow(tcg)
# run scores for every site/year:
for (s in site){
  for (i in years){
    # get observation:
    obs <- tcg[s, as.character(i)]
    # get ensemble for each site:
    ens <- pred[pred$site == s, as.character(i)]
    # if NA:0
    if (is.na(obs)|any(is.na(ens))){
      crps_scores[s, as.character(i)] <- NA
    } else {
    # calculate crps:
    crps <- scoringRules::crps_sample(y = obs, dat = ens)
    # fill in table:
    crps_scores[s, as.character(i)] <- crps
    }
  }
}


# # get observed value for each site:
# year = years[1]
# site = 1
# obs <- tcg[site, as.character(year)]
# ens <- testing[testing$site == 1, as.character(year)]
# crps <- scoringRules::crps_sample(y = obs, dat = ens)
  
# testing <- test %>%
#   group_by(site) %>% 
#  # summarise_all(., mean, na.rm = TRUE) %>%
#   summarise_all(list(min, max, mean, na.rm = TRUE))

# means <- apply(test, 2, mean, na.rm = TRUE)
# lower <- apply(test, 2, quantile, probs = c(0.05), na.rm = TRUE)
# upper <- apply(test, 2, quantile, probs = c(0.95), na.rm = TRUE)
