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
librarian::shelf(rjags, coda, dplyr, tidyverse)

# time series data:
tcg <- tcg_ts 

# make the disturbance magnitude version:
start <- grep("^2013", names(tcg))
end <- grep("^2015", names(tcg))
steady_state <- apply(tcg[,start:end], 1, mean)
dist <- grep("^2017", names(tcg))
post_dist <- tcg[,(dist + 1):(ncol(tcg)-1)]  # minus 1 here to account for match with daymet data 2018-2023 (2024 not fully available)
dm_post_dist <- steady_state - post_dist

# data for model:
recov_data <- dm_post_dist
# time series length:
time = 1:ncol(recov_data)
sites = 1:nrow(recov_data)
# first x:
x1 <- mean(tcg[,grep("^2017",names(tcg))], na.rm = T)

# covariates:
covs <- time_daym %>%
  # z-score normalizing (value-mean/sd): 
  mutate(across(-any_of(c("site", "year")), ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE))) %>%
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

# variables for running models in array:
vars <- unique(covs$variable)

# get time series for model:
cov_ts <- choose_covs(covs, vars[task_id])
model_name <- vars[task_id]


# model with time-varying covariate added: 
recov_state_space_uni_tv <- "model {
for (s in sites){

### Data Model:
  for (t in 1:nt){
    y[s,t] ~ dnorm(x[s,t], tau_obs)
  }

### Process Model:
for (t in 2:nt){
    R[s,t] <- r0 + atime[t-1] + beta*cov[s,t] ##+ asite[s] 
    mu[s,t] <- R[s,t] * x[s,t-1]  
    x[s,t] ~ dnorm(mu[s,t], tau_add)
  }
  x[s,1] ~ dnorm(x_ic, t_ic)
}

# ### Random Effects:
# for (s in sites){
#   asite[s] ~ dnorm(0, tausite)
# }


atime[1] = 0                   # option 2: indexing for atime[0]
for (t in 2:(nt-1)){
  atime[t] ~ dnorm(0, tautime)
}

### Priors:
r0 ~ dnorm(r_ic, r_prec)  # initial condition r
beta ~ dnorm(b0, Vb) #initial beta
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
#tausite ~ dgamma(0.001, 0.001)
tautime ~ dgamma(0.001, 0.001)
}
"

# make list object
model_data <- list(y = recov_data,
                   cov = cov_ts,
                   nt = length(time),
                   sites = sites, 
                   t_obs = 0.001, a_obs = 0.001,
                   t_add = 0.001, a_add = 0.001,
                   r_ic = 1, r_prec = 0.001,
                   x_ic = x1, t_ic = 0.01,
                   b0 = 0, Vb = 0.001)

# model run:
#jags_model <- jags.model(file = textConnection(recov_state_space_uni_tv),
#                         data = model_data,
#                         n.chains = 3)
#model test:
#jags_out <- coda.samples(jags_model,
#                         variable.names = c("x", "R",
#                                            "tau_obs", "tau_add",
#                                            "r0", "atime", #"asite",
#                                            "beta"),
#                         n.iter = 150000,
#                         adapt = 50000,
#                         thin = 50)

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
  filename_outputs <- paste0(out_path, date, "_model_tv_cov_uni_", model_name, "_output.csv")
  filename_runs <- paste0(run_path, date, "_model_tv_cov_uni_", model_name, "_data.RData")
  
  # save output
  write.csv(out, file = filename_outputs)
  
  # save model selection and metadata to folder
  model_info <- model_output[c('jags_out', 'dic', 'metadata')]
  save(model_info, file = filename_runs)
}

state_space_model_run(model_data = model_data,
                      model = recov_state_space_uni_tv,
                      model_name = model_name)
