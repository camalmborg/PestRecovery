# recovery rate state space model script for non-static (time-varying) recovery models
# this script contains the code for running a jags model for estimating recovery rates

# there will be 4 models, one for time-varying variable for univariate runs - growing season
# means of each environmental variable: precip, max temp, min temp, vpd
# task_ids will be 1-4

dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)

# setting task id for running array job on SCC:
task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# load environment if needed:
load("Environments/2025_07_07_environment.RData")

# libraries:
librarian::shelf(rjags, coda, dplyr)

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

# covariates:
covs <- time_daym %>%
  # pivot variables to make a variable column:
  pivot_longer(cols = c(prcp, tmax, tmin, vp), 
               names_to = "variable", 
               values_to = "value") %>%
  # pivot wide to make time series for each row
  pivot_wider(names_from = year, 
              values_from = value) %>%
  # arrange by variable:
  arrange(variable, site)

# choosing variable time series function:
#'@param cov_df = dataframe object of covariate time series sorted by variable and site
#'@param var = variable you would like to use for analysis
choose_covs <- function(cov_df, var){
  cov_df <- cov_df %>%
    filter(variable == var) %>%
    ungroup() %>%
    select(starts_with("2"))
}


