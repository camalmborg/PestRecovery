### Multivariate Recovery State Space JAGS model script ###

dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)

## Load libraries
librarian::shelf(rjags, coda, dplyr, tidyverse)

## Load data
# load environment if needed:
load("Environments/2025_07_07_environment.RData")

## Load models

## Prepare data for models
# time series data:
tcg <- tcg_ts 

# make the disturbance magnitude version:
start <- grep("^2013", names(tcg))
end <- grep("^2015", names(tcg))
steady_state <- apply(tcg[,start:end], 1, mean)
dist <- grep("^2017", names(tcg))
post_dist <- tcg[,(dist + 1):ncol(tcg)]
dm_post_dist <- steady_state - post_dist

# data for model:
recov_data <- as.matrix(dm_post_dist)
# time series length:
time = 1:ncol(recov_data)
sites = 1:nrow(recov_data)
# first x:
x1 <- mean(tcg[,grep("^2017",names(tcg))], na.rm = T)

## Prepare covariate data:
pre_dist_covs <- data.frame(#lat = coords$lat, lon = coords$lon,
  # adding disturbance history:
  dmags_tcg_y2 = dist_hist_tcg$`2017-05-01`,
  dmag_tcg_sum = dist_hist_tcg$dist_2_yrs,
  # adding environmental variables:
  precip_2015 = static_daym$prcp[which(static_daym$year == 2015)]) %>%
  # z-score normalizing (value-mean/sd): 
  mutate(across(-(grep("^cat_", colnames(covs))), ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)))
  
post_dist_covs <- time_daym %>%
  select(-tmin) %>%
  # z-score normalizing (value-mean/sd): 
  mutate(across(-any_of(c("site", "year")), ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE))) %>%
  # pivot variables to make a variable column:
  pivot_longer(cols = c(prcp, tmax, vp),
               names_to = "variable",
               values_to = "value") %>%
  # pivot wide to make time series for each row
  pivot_wider(names_from = year,
              values_from = value) %>%
  # arrange by variable:
  arrange(variable, site)

## Model runs
# setting task id for cluster runs:
task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
