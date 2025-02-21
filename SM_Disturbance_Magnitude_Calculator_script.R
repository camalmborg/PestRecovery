### This is the script for Harvard Forest 2022 field site disturbance magnitude and disturbance occurrance
### The script calculates a disturbance magnitude for the years of a disturbance, in this case 2016-2017 (into 2018)

#### ----- Libraries ----- ####
#install.packages("librarian")
librarian::shelf(tidyverse, dplyr, googledrive)

#### ----- Load and set up data ----- ####
# if not already in environment:
hf_scores <- read.csv("Data/hf_plot_scores_clean.csv")
hf_tcg <- read.csv("Data/hf_plot_tcg_clean.csv")

# make the calculator function:
#'@param ts = time series data for detecting disturbance in canopy greenness
#'@param distyr = onset year of disturbance
dist_mag_calc <- function(ts, distyr){
  # separate columns with just the canopy observation data:
  series <- ts[,grep("^2", names(ts))]
  # get disturbance onset column:
  dist <- grep(as.character(distyr), colnames(series))
  # calculate a pre-disturbance steady state to compare disturbance condition:
  steady <- apply(series[,(dist-6):(dist-1)], 1, mean, na.rm = T)
  # calculate disturbance magnitude:
  
}

# Update timeline
# 2025-02-20 created script