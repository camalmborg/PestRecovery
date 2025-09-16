# recovery rate state space model script for base and random effect models
# this script contains the code for running a jags model for estimating recovery rates
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)

# load environment if needed:
#load("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Environments/2025_06_12_environment.RData")
load("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Environments/2025_06_29_environment.RData")

# libraries:
librarian::shelf(rjags, coda)

# time series data:
tcg <- tcg_ts 

# make the disturbance magnitude version:
start <- grep("^2013", names(tcg))
end <- grep("^2015", names(tcg))
steady_state <- apply(tcg[,start:end], 1, mean)
dist <- grep("^2017", names(tcg))

# data for model:
r_start <- grep("^2018", names(tcg))
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
x_prec <- 1/(sd_tcg[,grep("2017", names(sd_tcg))])
x_prec[which(is.na(x_prec))] <- 1/x_prec_miss
x_prec[which(x_prec == Inf)] <- 1/x_prec_miss

# make list object
model_data <- list(y = recov_data,
                   nt = length(time),
                   sites = sites, 
                   t_obs = 0.001, a_obs = 0.001,
                   t_add = 0.001, a_add = 0.001,
                   r_ic = 1, r_prec = 0.001,
                   x_ic = x1, t_ic = x_prec)


### "Base Model" - no covariates or random effects:
recov_state_space <- "model {
for (s in sites){

### Data Model:
  for (t in 1:nt){
    y[s,t] ~ dnorm(x[s,t], tau_obs)
  }

### Process Model:
for (t in 2:nt){
    R[s,t] <- r0
    mu[s,t] <- R[s,t] * x[s,t-1]
    x[s,t] ~ dnorm(mu[s,t], tau_add)
  }
  x[s,1] ~ dnorm(x_ic[s], t_ic[s])
}

### Priors:
r0 ~ dnorm(r_ic, r_prec)  # initial condition r
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
}
"

# base model run:
jags_model <- jags.model(file = textConnection(recov_state_space),
                         data = model_data,
                         n.chains = 3)
jags_out <- coda.samples(jags_model,
                          variable.names = c("x", "R",
                                             "tau_obs", "tau_add",
                                             "r0"),
                          n.iter = 150000,
                          adapt = 50000,
                          thin = 50)

# run DIC
DIC <- dic.samples(jags_model, n.iter = 50000)
sum <- sum(DIC$deviance, DIC$penalty)

# Make output list
# track metadata:
metadata <- tibble::lst(recov_state_space, model_data)
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
filename_outputs <- paste0(out_path, date, "_base_model","_output.csv")
filename_runs <- paste0(run_path, date, "_base_model", "_data.RData")

# save output
write.csv(out, file = filename_outputs)

# save model selection and metadata to folder
model_info <- model_output[c('jags_out', 'dic', 'metadata')]
save(model_info, file = filename_runs)



### "Random Effects Model: option 2" - no covariates but includes random effects for site (pixel) and time (year)
recov_state_space_re_opt2 <- "model {
for (s in sites){

### Data Model:
  for (t in 1:nt){
    y[s,t] ~ dnorm(x[s,t], tau_obs)
  }

### Process Model:
for (t in 2:nt){
    R[s,t] <- r0 + asite[s] + atime[t-1]
    mu[s,t] <- R[s,t] * x[s,t-1]
    x[s,t] ~ dnorm(mu[s,t], tau_add)
  }
  x[s,1] ~ dnorm(x_ic[s], t_ic[s])
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
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
tausite ~ dgamma(0.001, 0.001)
tautime ~ dgamma(0.001, 0.001)
}
"

# base RE model run:
jags_model <- jags.model(file = textConnection(recov_state_space_re_opt2),
                         data = model_data,
                         n.chains = 3)
jags_out <- coda.samples(jags_model,
                         variable.names = c("atime", "asite", "x", "R",
                                            "tau_obs", "tau_add",
                                            "tausite", "tautime",
                                            "r0"),
                         n.iter = 150000,
                         n.adapt = 50000,
                         thin = 50)

# run DIC
DIC <- dic.samples(jags_model, n.iter = 50000)
sum <- sum(DIC$deviance, DIC$penalty)

# Make output list
# track metadata:
metadata <- tibble::lst(recov_state_space_re_opt2, model_data)
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
filename_outputs <- paste0(out_path, date, "_base_RE_model_opt2", "_output.csv")
filename_runs <- paste0(run_path, date, "_base_RE_model_opt2", "_data.RData")

# save output
write.csv(out, file = filename_outputs)

# save model selection and metadata to folder
model_info <- model_output[c('jags_out', 'dic', 'metadata')]
save(model_info, file = filename_runs)



### "Random Effects Model: JUST TIME" - no covariates but includes random effects for time (year)
recov_state_space_re_time <- "model {
for (s in sites){

### Data Model:
  for (t in 1:nt){
    y[s,t] ~ dnorm(x[s,t], tau_obs)
  }

### Process Model:
for (t in 2:nt){
    R[s,t] <- r0 + atime[t-1]  
    mu[s,t] <- R[s,t] * x[s,t-1]  
    x[s,t] ~ dnorm(mu[s,t], tau_add)
  }
  x[s,1] ~ dnorm(x_ic[s], t_ic[s])
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
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
#tausite ~ dgamma(0.001, 0.001)
tautime ~ dgamma(0.001, 0.001)
}
"

# base RE model run:
jags_model <- jags.model(file = textConnection(recov_state_space_re_time),
                         data = model_data, 
                         n.chains = 3)
jags_out <- coda.samples(jags_model,
                         variable.names = c("atime", "asite", "x", "R",
                                            "tau_obs", "tau_add",
                                            "tautime",
                                            "r0"),
                         n.iter = 150000,
                         n.adapt = 50000,
                         thin = 50)

# run DIC
DIC <- dic.samples(jags_model, n.iter = 50000)
sum <- sum(DIC$deviance, DIC$penalty)

# Make output list
# track metadata:
metadata <- tibble::lst(recov_state_space_re_time, model_data)
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
filename_outputs <- paste0(out_path, date, "_base_RE_time_model", "_output.csv")
filename_runs <- paste0(run_path, date, "_base_RE_time_model", "_data.RData")

# save output
write.csv(out, file = filename_outputs)

# save model selection and metadata to folder
model_info <- model_output[c('jags_out', 'dic', 'metadata')]
save(model_info, file = filename_runs)


### "Random Effects Model: SPACE" - no covariates but includes random effects for site (pixel)
recov_state_space_re_space <- "model {
for (s in sites){

### Data Model:
  for (t in 1:nt){
    y[s,t] ~ dnorm(x[s,t], tau_obs)
  }

### Process Model:
for (t in 2:nt){
    R[s,t] <- r0 + asite[s]  
    mu[s,t] <- R[s,t] * x[s,t-1]  
    x[s,t] ~ dnorm(mu[s,t], tau_add)
  }
  x[s,1] ~ dnorm(x_ic[s], t_ic[s])
}

### Random Effects: 
for (s in sites){
  asite[s] ~ dnorm(0, tausite)
}


# atime[1] = 0                   # option 2: indexing for atime[0]
# for (t in 2:(nt-1)){
#   atime[t] ~ dnorm(0, tautime)
# }

### Priors:
r0 ~ dnorm(r_ic, r_prec)  # initial condition r
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
tausite ~ dgamma(0.001, 0.001)
#tautime ~ dgamma(0.001, 0.001)
}
"

# base RE model run:
jags_model <- jags.model(file = textConnection(recov_state_space_re_space),
                         data = model_data, 
                         n.chains = 3)
jags_out <- coda.samples(jags_model,
                         variable.names = c("atime", "asite", "x", "R",
                                            "tau_obs", "tau_add",
                                            "tausite", 
                                            "r0"),
                         n.iter = 150000,
                         n.adapt = 50000,
                         thin = 50)

# run DIC
DIC <- dic.samples(jags_model, n.iter = 50000)
sum <- sum(DIC$deviance, DIC$penalty)

# Make output list
# track metadata:
metadata <- tibble::lst(recov_state_space_re_space, model_data)
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
filename_outputs <- paste0(out_path, date, "_base_RE_space_model", "_output.csv")
filename_runs <- paste0(run_path, date, "_base_RE_space_model", "_data.RData")

# save output
write.csv(out, file = filename_outputs)

# save model selection and metadata to folder
model_info <- model_output[c('jags_out', 'dic', 'metadata')]
save(model_info, file = filename_runs)

##### Archive ####
