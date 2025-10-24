### Script for looking at the outputs of the forecasts ###

## Load libraries
library(dplyr)
library(readr)
library(scoringRules)

## Load forecasts
# set working directory:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/"
setwd(dir)
  
# load model forecast files:
files <- list.files(pattern = "result.csv$")
model_num = as.numeric(Sys.getenv("SGE_TASK_ID"))
# get years of analysis:
start_year <- as.numeric(stringr::str_extract(files[model_num], "(?<=start_year_)\\d+"))
years <- start_year:2023
# loading predicted forecast values file:
model_out <- read.csv(files[model_num])
pred <- model_out %>%
  # rename columns with years:
  rename_with(~ as.character(years)[seq_along(.)], .cols = -1)

# load observation data:
tcg <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  # select columns with observations for 2017-2023:
  select(matches(as.character(years))) %>%
  # rename columns for years:
  rename_with(~ as.character(years)[seq_along(.)])

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

crps_means <- as.data.frame(crps_scores) %>%
  # take mean across sites:
  summarise_all(., mean, na.rm = TRUE)

## Saving results
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/CRPS/"
# select information from model name to match file names from forecasts:
model <- stringr::str_extract(files[model_num], "(?<=model_)\\d+")
# file name:
write.csv(crps_scores, paste0(save_dir, Sys.Date(), "_model_", model, "_start_year_", as.character(start_year), "_crps_result_all_sites.csv"))
write.csv(crps_means, paste0(save_dir, Sys.Date(), "_model_", model, "_start_year_", as.character(start_year), "_crps_result_across_site_means.csv"))



### ARCHIVE ###
# plot(1:ncol(crps_scores), log10(crps_scores[1,]), type = "l")
# for (i in 2:nrow(crps_scores)){
#   lines(1:ncol(crps_scores), log10(crps_scores[i,]), type = "l")
# }
  

