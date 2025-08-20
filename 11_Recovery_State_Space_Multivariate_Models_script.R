### Multivariate Recovery State Space JAGS model script ###

dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)

## Load libraries
librarian::shelf(rjags, coda, dplyr, tidyverse)

## Load data
# load environment if needed:
load("Environments/2025_07_07_environment.RData")

## Load models
tv_model 

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

# choosing time-varying variable time series function:
#'@param cov_df = dataframe object of covariate time series sorted by variable and site
#'@param var = variable you would like to use for analysis
choose_covs <- function(cov_df, var){
  cov_df <- cov_df %>%
    filter(variable == var) %>%
    ungroup() %>%
    select(starts_with("2"))
}
# time-varying variables:
tv_vars <- unique(covs$variable)

## Make list of model_data objects for entering in models:
# make list object
model_data <- list(y = recov_data,
                   #cov = cov_ts, ### CHANGE THIS
                   nt = length(time),
                   sites = sites, 
                   t_obs = 0.001, a_obs = 0.001,
                   t_add = 0.001, a_add = 0.001,
                   r_ic = 1, r_prec = 0.001,
                   x_ic = x1, t_ic = 0.01,
                   b0 = 0, Vb = 0.001)

## Making function for running models
# function for running model:----
#'@param model_data = list object with data for jags model
#'@param model = character - jags model
#'@param model_name = character - model name selected with task_id from variable list
state_space_model_run <- function(model_data, model, model_name){
  # model run:
  jags_model <- jags.model(file = textConnection(model),
                           data = model_data,
                           n.chains = 3)
  #model test:
  jags_out <- coda.samples(jags_model,
                           variable.names = c("x", "R",
                                              "tau_obs", "tau_add",
                                              "r0", "atime", "tautime", #"asite",
                                              "beta"),
                           n.iter = 150000,
                           adapt = 50000,
                           thin = 50)
  
  # run DIC
  DIC <- dic.samples(jags_model, n.iter = 50000)
  sum <- sum(DIC$deviance, DIC$penalty)
  
  # Make output list
  # track metadata:
  metadata <- tibble::lst(model, model_data)
  # model selection
  dic <- list(DIC, sum)
  # model output
  out <- as.matrix(jags_out)
  # combine output
  model_output <- tibble::lst(metadata, dic, jags_out, out)
  
  # save base model output:
  out_path <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/model_outputs/"
  run_path <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/model_runs/"
  date <- as.character(Sys.Date())
  filename_outputs <- paste0(out_path, date, "_model_cov_multi_", model_name, "_output.csv")
  filename_runs <- paste0(run_path, date, "_model_cov_multi_", model_name, "_data.RData")
  
  # save output
  write.csv(out, file = filename_outputs)
  
  # save model selection and metadata to folder
  model_info <- model_output[c('jags_out', 'dic', 'metadata')]
  save(model_info, file = filename_runs)
}

## Model runs
# setting task id for cluster runs:
task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# set model name:


state_space_model_run(model_data = model_data,
                      #model = CHANGE THIS!
                      model_name = model_name)

