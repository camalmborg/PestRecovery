# Daymet data preparing for recovery state space models

### libraries
librarian::shelf(dplyr)

### load data
data_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/daymet"
setwd(data_dir)
hf_daym <- read.csv("hf_daymet.csv")
samp_daym <- read.csv("samp_daymet.csv")

### Data for "static" covariates - covs not changing in time
# set the growing season months:

static_daym <- samp_daym %>%
  # get 