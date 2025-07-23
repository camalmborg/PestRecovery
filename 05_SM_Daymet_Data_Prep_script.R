# Daymet data preparing for recovery state space models

### libraries
librarian::shelf(tidyverse, dplyr, stringr)

### load data
data_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/daymet"
setwd(data_dir)
hf_daym <- read.csv("hf_daymet.csv")
samp_daym <- read.csv("samp_daymet.csv")
# coordinates:
# coords object loaded from 03_SM_Sample_Pixels_GEE_Data_cleaning

### Data for "static" covariates - covs not changing in time
# set the growing season months (may-sep):
gs <- c(5:9)

# get data:
seasonal_daym <- samp_daym %>%
  # rename columns:
  rename_with(~ str_extract(., "^[^\\.]+")) %>%
  # filter growing season months
  filter(month %in% gs) %>%
  # take mean of each variable during growing season:
  group_by(site, year) %>% summarise(across(c(prcp, tmax, tmin, vp), mean, na.rm = TRUE))

# separate into static and changing through time groups:
static_daym <- seasonal_daym %>% filter(year == c(2014, 2015))
time_daym <- seasonal_daym %>% filter(year > 2017)

# pivot time_daym wider to get time series:
daym_time_series <- time_daym %>%
  pivot_wider(names_from = c(year),
              values_from = c(prcp, tmax, tmin, vp))
