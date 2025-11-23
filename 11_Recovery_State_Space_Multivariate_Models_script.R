### Multivariate Recovery State Space JAGS model script ###

dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)

## Load libraries
librarian::shelf(rjags, coda, dplyr, tidyverse)

## Load data
# load environment if needed:
load("Environments/2025_07_07_environment.RData")

## Load models
# two time-varying covariates:
tv_model <- "model{
for (s in sites){

### Data Model:
  for (t in 1:nt){
    y[s,t] ~ dnorm(x[s,t], tau_obs)
  }

### Process Model:
for (t in 2:nt){
    R[s,t] <- r0 + atime[t-1] + beta[1]*cov_one[s,t-1] + beta[2]*cov_two[s,t-1]
    mu[s,t] <- R[s,t] * x[s,t-1]  
    x[s,t] ~ dnorm(mu[s,t], tau_add)
  }
  x[s,1] ~ dnorm(x_ic[s], t_ic[s])
}

### Random Effects:
atime[1] = 0                   # option 2: indexing for atime[0]
for (t in 2:(nt-1)){
  atime[t] ~ dnorm(0, tautime)
}

### Priors:
r0 ~ dnorm(r_ic, r_prec)  # initial condition r
beta[1] ~ dnorm(b0, Vb) #initial beta
beta[2] ~ dnorm(b00, Vbb)
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
tautime ~ dgamma(0.001, 0.001)
}"

# one time varying covariate and one static covariate with missing data:
tv_stat_miss_model <- "model{
for (s in sites){

### Data Model:
  for (t in 1:nt){
    y[s,t] ~ dnorm(x[s,t], tau_obs)
  }

### Process Model:
for (t in 2:nt){
    R[s,t] <- r0 + atime[t-1] + beta[1]*cov_one[s,t-1] + beta[2]*cov_two[s]
    mu[s,t] <- R[s,t] * x[s,t-1]  
    x[s,t] ~ dnorm(mu[s,t], tau_add)
  }
  x[s,1] ~ dnorm(x_ic[s], t_ic[s])
}

### Random Effects:
atime[1] = 0                   # option 2: indexing for atime[0]
for (t in 2:(nt-1)){
  atime[t] ~ dnorm(0, tautime)
}

### Priors:
r0 ~ dnorm(r_ic, r_prec)  # initial condition r
beta[1] ~ dnorm(b0, Vb) #initial beta
beta[2] ~ dnorm(b00, Vbb) #initial beta
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
tautime ~ dgamma(0.001, 0.001)
# missing data:
for (s in miss){
 cov_two[s] ~ dnorm(mis_s, mis_t)
  }
}"

## Prepare data for models
# time series data:
tcg <- tcg_ts 

# make the disturbance magnitude version:
start <- grep("^2013", names(tcg))
end <- grep("^2015", names(tcg))
steady_state <- apply(tcg[,start:end], 1, mean)
dist <- grep("^2017", names(tcg))

# data for model:
#recov_data <- as.matrix(dm_post_dist)
r_start <- grep("^2017", names(tcg))
r_end <- grep("^2023", names(tcg))
recov_data <- as.matrix(tcg[,r_start:r_end])
# time series length:
time = 1:ncol(recov_data)
sites = 1:nrow(recov_data)
# first x:
x_miss <- mean(tcg[,grep("^2017",names(tcg))], na.rm = T)
x_prec_miss <- sd(tcg[,grep("^2017", names(tcg))], na.rm = T)
x1 <- tcg[,grep("2017", names(tcg))]
x1[which(is.na(x1))] <- x_miss
# x1 precision:
x_prec <- 1/(sd_tcg[,grep("2017", names(sd_tcg))]^2)
x_prec[which(is.na(x_prec))] <- 1/x_prec_miss
x_prec[which(x_prec == Inf)] <- 1/x_prec_miss

## Prepare covariate data:
pre_dist_covs <- data.frame(#lat = coords$lat, lon = coords$lon,
  # adding disturbance history:
  dmags_tcg_y2 = dist_hist_tcg$`2017-05-01`,
  dmag_tcg_sum = dist_hist_tcg$dist_2_yrs) %>%
  # adding environmental variables:
  #precip_2015 = static_daym$prcp[which(static_daym$year == 2015)]) %>%
  # z-score normalizing (value-mean/sd): 
  mutate(across(-(grep("^cat_", colnames(covs))), ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)))
  
post_dist_covs <- time_daym %>%
  select(-c(tmax, tmin)) %>%
  # z-score normalizing (value-mean/sd): 
  mutate(across(-any_of(c("site", "year")), ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE))) %>%
  # pivot variables to make a variable column:
  pivot_longer(cols = c(prcp, vp),
               names_to = "variable",
               values_to = "value") %>%
  # pivot wide to make time series for each row
  pivot_wider(names_from = year,
              values_from = value) %>%
  # arrange by variable:
  arrange(variable, site) %>%
  # select 2017 onwards:
  select(-c("2016", "2017"))

# year lag covariates:
year_lag_covs <- time_daym %>%
  #select(-c(tmax)) %>%
  rename(., vpd_year_lag = vp) %>%
  rename(., tmax_year_lag = tmax) %>%
  rename(., tmin_year_lag = tmin) %>%
  rename(., prcp_year_lag = prcp) %>%
  # z-score normalizing (value-mean/sd): 
  mutate(across(-any_of(c("site", "year")), ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE))) %>%
  # pivot variables to make a variable column:
  pivot_longer(cols = c(vpd_year_lag, tmax_year_lag, tmin_year_lag, prcp_year_lag),
               names_to = "variable",
               values_to = "value") %>%
  # pivot wide to make time series for each row
  pivot_wider(names_from = year,
              values_from = value) %>%
  # arrange by variable:
  arrange(variable, site) %>%
  # select 2017 onwards:
  select(-c("2016", "2023"))


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
tv_vars <- unique(post_dist_covs$variable)  # 1 = prcp, 2 = tmax, 3 = vp
tv_yl_vars <- unique(year_lag_covs$variable)


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
                                              "r0", "atime", "tautime", 
                                              "beta"),
                           n.iter = 200000,
                           adapt = 75000,
                           thin = 100)
  
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

## Model runs:
# covariate lists:
model_covariates <- list(prcp_year_lag = choose_covs(year_lag_covs, tv_yl_vars[1]),
                         tmax_year_lag = choose_covs(year_lag_covs, tv_yl_vars[2]),
                         tmin_year_lag = choose_covs(year_lag_covs, tv_yl_vars[3]),
                         vpd_year_lag = choose_covs(year_lag_covs, tv_yl_vars[4]),
                         prcp = choose_covs(post_dist_covs, tv_vars[1]),
                         vpd = choose_covs(post_dist_covs, tv_vars[2]),
                         dmagy2 = pre_dist_covs$dmags_tcg_y2,
                         dmagsum = pre_dist_covs$dmag_tcg_sum)
# # model names:
# model_name <- c("prcpyrlag_tmaxyrlag", "prcpyrlag_vpd", "prcpyrlag_tminyrlag", "prcpyrlag_prcp", 
#                 "prcpyrlag_dmagy2", "prcpyrlag_dmagsum", "prcpyrlag_vpdyrlag", "tmaxyrlay_prcp", 
#                 "tmaxyrlag_dmagy2", "tmaxyrlag_vpd", "tmaxyrlag_dmagsum", "tmaxyrlag_vpdyrlag",
#                 "tmaxyrlag_tminyrlag")
# 
# # model data lists:
# input_data_list <- list()
# input_data_list[[1]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$tmax_year_lag)
# input_data_list[[2]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$vpd)
# input_data_list[[3]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$tmin_year_lag)
# input_data_list[[4]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$prcp)
# input_data_list[[5]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$dmagy2)
# input_data_list[[6]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$dmagsum)
# input_data_list[[7]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$vpd_year_lag)
# input_data_list[[8]] <- list(cov_one = model_covariates$tmax_year_lag, cov_two = model_covariates$prcp)
# input_data_list[[9]] <- list(cov_one = model_covariates$tmax_year_lag, cov_two = model_covariates$dmagy2)
# input_data_list[[10]] <- list(cov_one = model_covariates$tmax_year_lag, cov_two = model_covariates$vpd)
# input_data_list[[11]] <- list(cov_one = model_covariates$tmax_year_lag, cov_two = model_covariates$dmagsum)
# input_data_list[[12]] <- list(cov_one = model_covariates$tmax_year_lag, cov_two = model_covariates$vpd_year_lag)
# input_data_list[[13]] <- list(cov_one = model_covariates$tmax_year_lag, cov_two = model_covariates$tmin_year_lag)
# 
# # missing data models:
# missing <- c(5, 6, 9, 11)
# 
# # models:
# model_list <- list()
# model_list[[1]] <- tv_model
# model_list[[2]] <- tv_model
# model_list[[3]] <- tv_model
# model_list[[4]] <- tv_model
# model_list[[5]] <- tv_stat_miss_model
# model_list[[6]] <- tv_stat_miss_model
# model_list[[7]] <- tv_model
# model_list[[8]] <- tv_model
# model_list[[9]] <- tv_stat_miss_model
# model_list[[10]] <- tv_model
# model_list[[11]] <- tv_stat_miss_model
# model_list[[12]] <- tv_model
# model_list[[13]] <- tv_model
# 
# # setting task id for cluster runs:
# task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# 
# # model data object:
# if (task_id %in% missing) {
#   # model data object for missing data:
#   model_data <- list(y = recov_data,
#                      cov_one = input_data_list[[task_id]]$cov_one,
#                      cov_two = input_data_list[[task_id]]$cov_two,
#                      nt = length(time),
#                      sites = sites,
#                      t_obs = 0.001, a_obs = 0.001,
#                      t_add = 0.001, a_add = 0.001,
#                      r_ic = 1, r_prec = 0.001,
#                      x_ic = x1, t_ic = x_prec,
#                      b0 = 0, Vb = 0.001,
#                      b00 = 0, Vbb = 0.001)
#   # missing data:
#   model_data$miss <- which(is.na(pre_dist_covs$dmags_tcg_y2))
#   model_data$mis_s = mean(dmag_data$steady, na.rm = T)
#   model_data$mis_t = 0.01
# } else {
#   model_data <- list(y = recov_data,
#                      cov_one = input_data_list[[task_id]]$cov_one,
#                      cov_two = input_data_list[[task_id]]$cov_two,
#                      nt = length(time),
#                      sites = sites,
#                      t_obs = 0.001, a_obs = 0.001,
#                      t_add = 0.001, a_add = 0.001,
#                      r_ic = 1, r_prec = 0.001,
#                      x_ic = x1, t_ic = x_prec,
#                      b0 = 0, Vb = 0.001,
#                      b00 = 0, Vbb = 0.001)
# }

# jags_run <- model_list[[task_id]]
# 
# # run the models according to task_id:
# state_space_model_run(model_data = model_data,
#                       model = jags_run,
#                       model_name = model_name[task_id])


# ### Run a 3-variable multivariate run:
# tv_tv_tv_model <- "model{
# for (s in sites){
# 
# ### Data Model:
#   for (t in 1:nt){
#     y[s,t] ~ dnorm(x[s,t], tau_obs)
#   }
# 
# ### Process Model:
# for (t in 2:nt){
#     R[s,t] <- r0 + atime[t-1] + beta[1]*cov_one[s,t-1] + beta[2]*cov_two[s, t-1] + beta[3]*cov_three[s, t-1]
#     mu[s,t] <- R[s,t] * x[s,t-1]
#     x[s,t] ~ dnorm(mu[s,t], tau_add)
#   }
#   x[s,1] ~ dnorm(x_ic[s], t_ic[s])
# }
# 
# ## Random Effects:
# atime[1] = 0                   # option 2: indexing for atime[0]
# for (t in 2:(nt-1)){
#   atime[t] ~ dnorm(0, tautime)
# }
# 
# ### Priors:
# r0 ~ dnorm(r_ic, r_prec)  # initial condition r
# beta[1] ~ dnorm(b0, Vb) #initial beta
# beta[2] ~ dnorm(b00, Vbb) #initial beta 2
# beta[3] ~ dnorm(b000, Vbbb) #initial beta 3
# tau_obs ~ dgamma(t_obs, a_obs)
# tau_add ~ dgamma(t_add, a_add)
# tautime ~ dgamma(0.001, 0.001)
# }
# "
# # model names:
# model_name <- c("prcpyrlag_prcp_tmaxyrlag", "prcpyrlag_prcp_vpd", "prcpyrlag_prcp_tminyrlag", "prcpyrlag_tmaxyrlag_vpd")
# 
# # model data lists:
# input_data_list <- list()
# input_data_list[[1]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$prcp, cov_three = model_covariates$tmax_year_lag)
# input_data_list[[2]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$prcp, cov_three = model_covariates$vpd)
# input_data_list[[3]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$prcp, cov_three = model_covariates$tmin_year_lag)
# input_data_list[[4]] <- list(cov_one = model_covariates$prcp_year_lag, cov_two = model_covariates$tmax_year_lag, cov_three = model_covariates$vpd)
# 
# # models:
# model_list <- list()
# model_list[[1]] <- tv_tv_tv_model
# model_list[[2]] <- tv_tv_tv_model
# model_list[[3]] <- tv_tv_tv_model
# model_list[[4]] <- tv_tv_tv_model
# 
# # setting task id for cluster runs:
# #task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# 
# # model data object for missing data (all models in this run):
# model_data <- list(y = recov_data,
#                    cov_one = input_data_list[[task_id]]$cov_one,
#                    cov_two = input_data_list[[task_id]]$cov_two,
#                    cov_three = input_data_list[[task_id]]$cov_three,
#                    nt = length(time),
#                    sites = sites,
#                    t_obs = 0.001, a_obs = 0.001,
#                    t_add = 0.001, a_add = 0.001,
#                    r_ic = 1, r_prec = 0.001,
#                    x_ic = x1, t_ic = x_prec,
#                    b0 = 0, Vb = 0.001,
#                    b00 = 0, Vbb = 0.001,
#                    b000 = 0, Vbbb = 0.001)

# # missing data:
# model_data$miss <- which(is.na(pre_dist_covs$dmags_tcg_y2))
# model_data$mis_s = mean(dmag_data$steady, na.rm = T)
# model_data$mis_t = 0.01

# # model:
# jags_run <- model_list[[task_id]]
# 
# # run the models according to task_id:
# state_space_model_run(model_data = model_data,
#                       model = jags_run,
#                       model_name = model_name[task_id])


# ## Run a 4-var model
# tv_tv_tv_tv_model <- "model{
# for (s in sites){
# 
# ### Data Model:
#   for (t in 1:nt){
#     y[s,t] ~ dnorm(x[s,t], tau_obs)
#   }
# 
# ### Process Model:
# for (t in 2:nt){
#     R[s,t] <- r0 + atime[t-1] + beta[1]*cov_one[s,t-1] + beta[2]*cov_two[s, t-1] + beta[3]*cov_three[s, t-1] + beta[4]*cov_four[s, t-1] + beta[5]*cov_five[s, t-1]
#     mu[s,t] <- R[s,t] * x[s,t-1]
#     x[s,t] ~ dnorm(mu[s,t], tau_add)
#   }
#   x[s,1] ~ dnorm(x_ic[s], t_ic[s])
# }
# 
# ### Random Effects:
# atime[1] = 0                   # option 2: indexing for atime[0]
# for (t in 2:(nt-1)){
#   atime[t] ~ dnorm(0, tautime)
# }
# 
# ### Priors:
# r0 ~ dnorm(r_ic, r_prec)  # initial condition r
# beta[1] ~ dnorm(b0, Vb) #initial beta
# beta[2] ~ dnorm(b00, Vbb) #initial beta 2
# beta[3] ~ dnorm(b000, Vbbb) #initial beta 3
# beta[4] ~dnorm(b0000, Vbbbb) #initial beta 4
# beta[5] ~ dnorm(b00000, Vbbbbb) # initial beta 5
# tau_obs ~ dgamma(t_obs, a_obs)
# tau_add ~ dgamma(t_add, a_add)
# tautime ~ dgamma(0.001, 0.001)
# 
# }
# "
# 
# # model names:
# model_name <- c("prcpyrlag_prcp_tmaxyrlag_vpd_tminyrlag")
# 
# # model data lists:
# input_data_list <- list()
# input_data_list[[1]] <- list(cov_one = model_covariates$prcp_year_lag,
#                              cov_two = model_covariates$prcp,
#                              cov_three = model_covariates$tmax_year_lag,
#                              cov_four = model_covariates$vpd,
#                              cov_five = model_covariates$tmin_year_lag)
# 
# # models:
# model_list <- list()
# model_list[[1]] <- tv_tv_tv_tv_model
# 
# # setting task id for cluster runs:
# #task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# task_id = 1
# 
# # model data object for missing data (all models in this run):
# model_data <- list(y = recov_data,
#                    cov_one = input_data_list[[task_id]]$cov_one,
#                    cov_two = input_data_list[[task_id]]$cov_two,
#                    cov_three = input_data_list[[task_id]]$cov_three,
#                    cov_four = input_data_list[[task_id]]$cov_four,
#                    cov_five = input_data_list[[task_id]]$cov_five,
#                    nt = length(time),
#                    sites = sites,
#                    t_obs = 0.001, a_obs = 0.001,
#                    t_add = 0.001, a_add = 0.001,
#                    r_ic = 1, r_prec = 0.001,
#                    x_ic = x1, t_ic = x_prec,
#                    b0 = 0, Vb = 0.001,
#                    b00 = 0, Vbb = 0.001,
#                    b000 = 0, Vbbb = 0.001,
#                    b0000 = 0, Vbbbb = 0.001,
#                    b00000 = 0, Vbbbbb = 0.001)
# # # missing data:
# # model_data$miss <- which(is.na(pre_dist_covs$dmags_tcg_y2))
# # model_data$mis_s = mean(dmag_data$steady, na.rm = T)
# # model_data$mis_t = 0.01
# 
# # model:
# jags_run <- model_list[[task_id]]
# 
# # run the models according to task_id:
# state_space_model_run(model_data = model_data,
#                       model = jags_run,
#                       model_name = model_name[task_id])




### ARCHIVE ###

# # one time-varying covariate, one static covariate:
# tv_stat_model <- "model{
# for (s in sites){
# 
# ### Data Model:
#   for (t in 1:nt){
#     y[s,t] ~ dnorm(x[s,t], tau_obs)
#   }
# 
# ### Process Model:
# for (t in 2:nt){
#     R[s,t] <- r0 + atime[t-1] + beta[1]*cov_one[s,t-1] + beta[2]*cov_two[s]
#     mu[s,t] <- R[s,t] * x[s,t-1]  
#     x[s,t] ~ dnorm(mu[s,t], tau_add)
#   }
#   x[s,1] ~ dnorm(x_ic[s], t_ic[s])
# }
# 
# ### Random Effects:
# atime[1] = 0                   # option 2: indexing for atime[0]
# for (t in 2:(nt-1)){
#   atime[t] ~ dnorm(0, tautime)
# }
# 
# ### Priors:
# r0 ~ dnorm(r_ic, r_prec)  # initial condition r
# beta[1] ~ dnorm(b0, Vb) #initial beta
# beta[2] ~ dnorm(b00, Vbb)
# tau_obs ~ dgamma(t_obs, a_obs)
# tau_add ~ dgamma(t_add, a_add)
# tautime ~ dgamma(0.001, 0.001)
# }"
# 
# # two static covariates:
# stat_model <- "model{
# for (s in sites){
# 
# ### Data Model:
#   for (t in 1:nt){
#     y[s,t] ~ dnorm(x[s,t], tau_obs)
#   }
# 
# ### Process Model:
# for (t in 2:nt){
#     R[s,t] <- r0 + atime[t-1] + beta[1]*cov_one[s] + beta[2]*cov_two[s]
#     mu[s,t] <- R[s,t] * x[s,t-1]  
#     x[s,t] ~ dnorm(mu[s,t], tau_add)
#   }
#   x[s,1] ~ dnorm(x_ic[s], t_ic[s])
# }
# 
# # Random Effects:
# atime[1] = 0                   # option 2: indexing for atime[0]
# for (t in 2:(nt-1)){
#   atime[t] ~ dnorm(0, tautime)
# }
# 
# ### Priors:
# r0 ~ dnorm(r_ic, r_prec)  # initial condition r
# beta[1] ~ dnorm(b0, Vb) #initial beta
# beta[2] ~ dnorm(b00, Vbb) #initial beta
# tau_obs ~ dgamma(t_obs, a_obs)
# tau_add ~ dgamma(t_add, a_add)
# tautime ~ dgamma(0.001, 0.001)
# }"
# 
# # two static covariates with one having missing data:
# stat_miss_model <- "model{
# for (s in sites){
# 
# ### Data Model:
#   for (t in 1:nt){
#     y[s,t] ~ dnorm(x[s,t], tau_obs)
#   }
# 
# ### Process Model:
# for (t in 2:nt){
#     R[s,t] <- r0 + atime[t-1] + beta[1]*cov_one[s] + beta[2]*cov_two[s]
#     mu[s,t] <- R[s,t] * x[s,t-1]  
#     x[s,t] ~ dnorm(mu[s,t], tau_add)
#   }
#   x[s,1] ~ dnorm(x_ic[s], t_ic[s])
# }
# 
# ### Random Effects:
# atime[1] = 0                   # option 2: indexing for atime[0]
# for (t in 2:(nt-1)){
#   atime[t] ~ dnorm(0, tautime)
# }
# 
# ### Priors:
# r0 ~ dnorm(r_ic, r_prec)  # initial condition r
# beta[1] ~ dnorm(b0, Vb) #initial beta
# beta[2] ~ dnorm(b00, Vbb) #initial beta
# tau_obs ~ dgamma(t_obs, a_obs)
# tau_add ~ dgamma(t_add, a_add)
# tautime ~ dgamma(0.001, 0.001)
# # missing data:
# for (s in miss){
#  cov_one[s] ~ dnorm(mis_s, mis_t)
#   }
# }"
# 
# tv_stat_stat_miss_model <- "model{
# for (s in sites){
# 
# ### Data Model:
#   for (t in 1:nt){
#     y[s,t] ~ dnorm(x[s,t], tau_obs)
#   }
# 
# ### Process Model:
# for (t in 2:nt){
#     R[s,t] <- r0 + atime[t-1] + beta[1]*cov_one[s,t-1] + beta[2]*cov_two[s] + beta[3]*cov_three[s]
#     mu[s,t] <- R[s,t] * x[s,t-1]
#     x[s,t] ~ dnorm(mu[s,t], tau_add)
#   }
#   x[s,1] ~ dnorm(x_ic[s], t_ic[s])
# }
# 
# ### Random Effects:
# atime[1] = 0                   # option 2: indexing for atime[0]
# for (t in 2:(nt-1)){
#   atime[t] ~ dnorm(0, tautime)
# }
# 
# ### Priors:
# r0 ~ dnorm(r_ic, r_prec)  # initial condition r
# beta[1] ~ dnorm(b0, Vb) #initial beta
# beta[2] ~ dnorm(b00, Vbb) #initial beta 2
# beta[3] ~ dnorm(b000, Vbbb) #initial beta 3
# tau_obs ~ dgamma(t_obs, a_obs)
# tau_add ~ dgamma(t_add, a_add)
# tautime ~ dgamma(0.001, 0.001)
# # missing data:
# for (s in miss){
#  cov_two[s] ~ dnorm(mis_s, mis_t)
#   }
# }"

# tv_tv_stat_miss_model <- "model{
# for (s in sites){
# 
# ### Data Model:
#   for (t in 1:nt){
#     y[s,t] ~ dnorm(x[s,t], tau_obs)
#   }
# 
# ### Process Model:
# for (t in 2:nt){
#     R[s,t] <- r0 + atime[t-1] + beta[1]*cov_one[s,t-1] + beta[2]*cov_two[s, t-1] + beta[3]*cov_three[s]
#     mu[s,t] <- R[s,t] * x[s,t-1]
#     x[s,t] ~ dnorm(mu[s,t], tau_add)
#   }
#   x[s,1] ~ dnorm(x_ic[s], t_ic[s])
# }
# 
# ## Random Effects:
# atime[1] = 0                   # option 2: indexing for atime[0]
# for (t in 2:(nt-1)){
#   atime[t] ~ dnorm(0, tautime)
# }
# 
# ### Priors:
# r0 ~ dnorm(r_ic, r_prec)  # initial condition r
# beta[1] ~ dnorm(b0, Vb) #initial beta
# beta[2] ~ dnorm(b00, Vbb) #initial beta 2
# beta[3] ~ dnorm(b000, Vbbb) #initial beta 3
# tau_obs ~ dgamma(t_obs, a_obs)
# tau_add ~ dgamma(t_add, a_add)
# tautime ~ dgamma(0.001, 0.001)
# # missing data:
# for (s in miss){
#  cov_three[s] ~ dnorm(mis_s, mis_t)
#   }
# }
# "