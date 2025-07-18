# recovery rate state space model script for base and random effect models
# this script contains the code for running a jags model for estimating recovery rates
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)

# load environment if needed:
#load("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Environments/2025_06_12_environment.RData")
#load("Environments/2025_06_29_environment.RData")
#load("Environments/2025_07_03_environment.RData")
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

# covariate data:
covs <- data.frame(#lat = coords$lat, lon = coords$lon,
  # adding disturbance history:
  dmags_tcg_y1 = dist_hist_tcg$`2016-05-01`,
  dmags_tcg_y2 = dist_hist_tcg$`2017-05-01`,
  dmag_tcg_sum = dist_hist_tcg$dist_2_yrs,
  # adding environmental variables:
  precip_2014 = static_daym$prcp[which(static_daym$year == 2014)],
  precip_2015 = static_daym$prcp[which(static_daym$year == 2015)],
  temp_mx_2014 = static_daym$tmax[which(static_daym$year == 2014)],
  temp_mx_2015 = static_daym$tmax[which(static_daym$year == 2015)],
  temp_min_2014 = static_daym$tmin[which(static_daym$year == 2014)],
  temp_min_2015 = static_daym$tmin[which(static_daym$year == 2015)],
  vpd_2014 = static_daym$vp[which(static_daym$year == 2014)],
  vpd_2015 = static_daym$vp[which(static_daym$year == 2015)],
  # adding NLCD data:
  pct_tree_cov = nlcd$percent_tree_cover,
  nlcd_cat[,grep("^cat_", colnames(nlcd_cat))]) %>%
  # adding means of multiple pre-disturbance years for environmental variables:
  mutate(precip_mean = rowMeans(select(., precip_2014, precip_2015)), .after = precip_2015) %>%
  mutate(temp_mx_mean = rowMeans(select(., temp_mx_2014, temp_mx_2015)), .after = temp_mx_2015) %>%
  mutate(temp_min_mean = rowMeans(select(., temp_min_2014, temp_min_2015)), .after = temp_min_2015) %>%
  mutate(vpd_mean = rowMeans(select(., vpd_2014, vpd_2015)), .after = vpd_2015) %>%
  # z-score normalizing (value-mean/sd): 
  mutate(across(-c(grep("^cat_", colnames(covs))), ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)))

# continuous variables:
stat_covs <- covs[-c(grep("^cat_", colnames(covs)))]

# separate categorical variables:
cat_covs <- covs[,grep("^cat_", colnames(covs))]

# model with covariate added: 
recov_state_space_uni_static <- "model {
for (s in sites){

### Data Model:
  for (t in 1:nt){
    y[s,t] ~ dnorm(x[s,t], tau_obs)
  }

### Process Model:
for (t in 2:nt){
    R[s,t] <- r0 + asite[s] + atime[t-1] + beta*cov[s]
    mu[s,t] <- R[s,t] * x[s,t-1]  
    x[s,t] ~ dnorm(mu[s,t], tau_add)
  }
  x[s,1] ~ dnorm(x_ic, t_ic)
}

### Random Effects:
for (s in sites){
  asite[s] ~ dnorm(0, tausite)
}


atime[1] = 0                   # option 2: indexing for atime[0]
for (t in 2:(nt-1)){
  atime[t] ~ dnorm(0, tautime)
}

### Priors:
r0 ~ dnorm(r_ic, r_prec)  # initial condition r
beta ~ dnorm(b0, Vb) #initial beta
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
tausite ~ dgamma(0.001, 0.001)
tautime ~ dgamma(0.001, 0.001)
}
"

# model with missing data model: #####
recov_state_space_uni_static_miss <- "model {
for (s in sites){

### Data Model:
  for (t in 1:nt){
    y[s,t] ~ dnorm(x[s,t], tau_obs)
  }

### Process Model:
for (t in 2:nt){
    R[s,t] <- r0 + asite[s] + atime[t-1] + beta*cov[s]
    mu[s,t] <- R[s,t] * x[s,t-1]  
    x[s,t] ~ dnorm(mu[s,t], tau_add)
  }
  x[s,1] ~ dnorm(x_ic, t_ic)
}

### Random Effects:
for (s in sites){
  asite[s] ~ dnorm(0, tausite)
}


atime[1] = 0                   # option 2: indexing for atime[0]
for (t in 2:(nt-1)){
  atime[t] ~ dnorm(0, tautime)
}

### Priors:
r0 ~ dnorm(r_ic, r_prec)  # initial condition r
beta ~ dnorm(b0, Vb) #initial beta
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
tausite ~ dgamma(0.001, 0.001)
tautime ~ dgamma(0.001, 0.001)
# missing data:
for (s in miss){
 cov[s] ~ dnorm(mis_s, mis_t)
  }
}
"

# model for categorical variable:----
recov_state_space_uni_static_cat <- "model {
for (s in sites){

### Data Model:
  for (t in 1:nt){
    y[s,t] ~ dnorm(x[s,t], tau_obs)
  }

### Process Model:
for (t in 2:nt){
    R[s,t] <- r0 + asite[s] + atime[t-1] + inprod(beta[], cov[s,])
    mu[s,t] <- R[s,t] * x[s,t-1]  
    x[s,t] ~ dnorm(mu[s,t], tau_add)
  }
  x[s,1] ~ dnorm(x_ic, t_ic)
}

### Random Effects:
for (s in sites){
  asite[s] ~ dnorm(0, tausite)
}


atime[1] = 0                   # option 2: indexing for atime[0]
for (t in 2:(nt-1)){
  atime[t] ~ dnorm(0, tautime)
}

### Priors:
r0 ~ dnorm(r_ic, r_prec)  # initial condition r
beta ~ dmnorm(b0, Vb) #initial beta
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
tausite ~ dgamma(0.001, 0.001)
tautime ~ dgamma(0.001, 0.001)
}
"

# function for running model:----
#'@param cov_df = data frame object of covariates for particular run
#'@param model_num = numeric - which column of data frame to run: 1-16 for stat_covs, 1-5 for cat_covs
#'@param model = character - version of jags model to run
state_space_model_run <- function(cov_df, model_num){
  # get covs:
  if (ncol(cov_df) > 5){
    model_name = names(cov_df[model_num])
    covs = cov_df[,model_num]
  } else {
    model_name <- "nlcd_cat"
    covs = cov_df
  }
  
  # model choice:
  if(ncol(cov_df) <= 5){
    model <- recov_state_space_uni_static_cat
  } else if (length(which(is.na(covs))) > 0){
    model <- recov_state_space_uni_static_miss
  } else {
    model <- recov_state_space_uni_static
  }
  
  # make list object
  model_data <- list(y = recov_data,
                     cov = covs,
                     nt = length(time),
                     sites = sites, 
                     t_obs = 0.001, a_obs = 0.001,
                     t_add = 0.001, a_add = 0.001,
                     r_ic = 1, r_prec = 0.001,
                     x_ic = x1, t_ic = 0.01,
                     b0 = 0, Vb = 0.001)
  
  # model run:
  jags_model <- jags.model(file = textConnection(model),
                           data = model_data,
                           n.chains = 3)
  
  #model test:
  jags_out <- coda.samples(jags_model,
                           variable.names = c("x", "R",
                                              "tau_obs", "tau_add",
                                              "r0", "atime", "asite",
                                              "beta"),
                           n.iter = 150000,
                           adapt = 50000,
                           thin = 15)
  
  # run DIC
  DIC <- dic.samples(jags_model, n.iter = 50000)
  sum <- sum(DIC$deviance, DIC$penalty)
  
  # Make output list
  # track metadata:
  metadata <- tibble::lst(recov_state_space_uni_static, model_data)
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
  model_name <- names(covs)
  filename_outputs <- paste0(out_path, date, "_model_static_cov_uni_", model_name, "_output.csv")
  filename_runs <- paste0(run_path, date, "_model_static_cov_uni_", model_name, "_data.RData")
  
  # save output
  write.csv(out, file = filename_outputs)
  
  # save model selection and metadata to folder
  model_info <- model_output[c('jags_out', 'dic', 'metadata')]
  save(model_info, file = filename_runs)
}

# Model Run:


##### ARCHIVE #####

# #for running model with a sample of data
# # samp <- sort(sample(1:5000, 10))
# test_samp <- recov_data[samp,]
# cov_samp <- covs[samp,]
# 
# 
# # time series length:
# time = 1:ncol(test_samp)
# sites = 1:nrow(test_samp)

# for (i in 1:ncol(covs)){
#   nas <- which(is.na(covs[,i]))
#   print(nas)
# }
# 
# # run loop for all the models that aren't categorical:
# for (i in 9:ncol(covs[-c(grep("^cat_", colnames(covs)))])){
#   # make list object
#   model_data <- list(y = recov_data,
#                      cov = covs[,i],
#                      nt = length(time),
#                      sites = sites, 
#                      t_obs = 0.001, a_obs = 0.001,
#                      t_add = 0.001, a_add = 0.001,
#                      r_ic = 1, r_prec = 0.001,
#                      x_ic = x1, t_ic = 0.01,
#                      b0 = 0, Vb = 0.001)
#   
#   # model run:
#   jags_model <- jags.model(file = textConnection(recov_state_space_uni_static),
#                            data = model_data,
#                            n.chains = 3)
#   
#   #model test:
#   jags_out <- coda.samples(jags_model,
#                            variable.names = c("x", "R",
#                                               "tau_obs", "tau_add",
#                                               "r0", "atime", "asite",
#                                               "beta"),
#                            n.iter = 150000,
#                            adapt = 50000,
#                            thin = 15)
#   
#   # run DIC
#   DIC <- dic.samples(jags_model, n.iter = 50000)
#   sum <- sum(DIC$deviance, DIC$penalty)
#   
#   # Make output list
#   # track metadata:
#   metadata <- tibble::lst(recov_state_space_uni_static, model_data)
#   # model selection
#   dic <- list(DIC, sum)
#   # model output
#   out <- as.matrix(jags_out)
#   # combine output
#   model_output <- tibble::lst(metadata, dic, jags_out, out)
#   
#   # save base model output:
#   out_path <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/model_outputs/"
#   run_path <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/model_runs/"
#   date <- as.character(Sys.Date())
#   model_name <- names(covs)[i]
#   filename_outputs <- paste0(out_path, date, "_model_static_cov_uni_", model_name, "_output.csv")
#   filename_runs <- paste0(run_path, date, "_model_static_cov_uni_", model_name, "_data.RData")
#   
#   # save output
#   write.csv(out, file = filename_outputs)
#   
#   # save model selection and metadata to folder
#   model_info <- model_output[c('jags_out', 'dic', 'metadata')]
#   save(model_info, file = filename_runs)
# }
# 
# 
# ## Run the categorial:
# # model with covariate added:
# recov_state_space_uni_static_cat <- "model {
# for (s in sites){
# 
# ### Data Model:
#   for (t in 1:nt){
#     y[s,t] ~ dnorm(x[s,t], tau_obs)
#   }
# 
# ### Process Model:
# for (t in 2:nt){
#     R[s,t] <- r0 + asite[s] + atime[t-1] + inprod(beta[], cov[s,])
#     mu[s,t] <- R[s,t] * x[s,t-1]  
#     x[s,t] ~ dnorm(mu[s,t], tau_add)
#   }
#   x[s,1] ~ dnorm(x_ic, t_ic)
# }
# 
# ### Random Effects:
# for (s in sites){
#   asite[s] ~ dnorm(0, tausite)
# }
# 
# 
# atime[1] = 0                   # option 2: indexing for atime[0]
# for (t in 2:(nt-1)){
#   atime[t] ~ dnorm(0, tautime)
# }
# 
# ### Priors:
# r0 ~ dnorm(r_ic, r_prec)  # initial condition r
# beta ~ dmnorm(b0, Vb) #initial beta
# tau_obs ~ dgamma(t_obs, a_obs)
# tau_add ~ dgamma(t_add, a_add)
# tausite ~ dgamma(0.001, 0.001)
# tautime ~ dgamma(0.001, 0.001)
# }
# "
# 
# # separate categorical variables
# cat_covs <- covs[,grep("^cat_", colnames(covs))]
# 
# # make list object
# model_data <- list(y = recov_data,
#                    cov = cat_covs,
#                    nt = length(time),
#                    sites = sites, 
#                    t_obs = 0.001, a_obs = 0.001,
#                    t_add = 0.001, a_add = 0.001,
#                    r_ic = 1, r_prec = 0.001,
#                    x_ic = x1, t_ic = 0.01)
# # beta priors:
# model_data$b0 <- as.vector(c(rep(0, ncol(cat_covs))))      # beta means
# model_data$Vb <- solve(diag(1000, ncol(cat_covs)))   # beta precisions
# 
# # model run:
# jags_model <- jags.model(file = textConnection(recov_state_space_uni_static_cat),
#                          data = model_data,
#                          n.chains = 3)
# # model run:
# jags_out <- coda.samples(jags_model,
#                          variable.names = c("x", "R",
#                                             "tau_obs", "tau_add",
#                                             "r0", "atime", "asite",
#                                             "beta"),
#                          n.iter = 150000,
#                          adapt = 50000,
#                          thin = 15)
# 
# # run DIC
# DIC <- dic.samples(jags_model, n.iter = 50000)
# sum <- sum(DIC$deviance, DIC$penalty)
# 
# # Make output list
# # track metadata:
# metadata <- tibble::lst(recov_state_space_uni_static_cat, model_data)
# # model selection
# dic <- list(DIC, sum)
# # model output
# out <- as.matrix(jags_out)
# # combine output
# model_output <- tibble::lst(metadata, dic, jags_out, out)
# 
# # save base model output:
# out_path <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/model_outputs/"
# run_path <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/model_runs/"
# date <- as.character(Sys.Date())
# model_name <- "forest_type"
# filename_outputs <- paste0(out_path, date, "_model_static_cov_uni_", model_name, "_output.csv")
# filename_runs <- paste0(run_path, date, "_model_static_cov_uni_", model_name, "_data.RData")
# 
# # save output
# write.csv(out, file = filename_outputs)
# 
# # save model selection and metadata to folder
# model_info <- model_output[c('jags_out', 'dic', 'metadata')]
# save(model_info, file = filename_runs)
