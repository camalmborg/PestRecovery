### Script for looking at the outputs of the forecasts ###

## Load libraries
library(dplyr)
library(readr)

## Load forecasts
# set working directory:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/"
setwd(dir)
# load files:
