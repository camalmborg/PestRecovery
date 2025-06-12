# recovery rate state space model script
# this script contains the code for running a jags model for estimating recovery rates

# load environment if needed:
load("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Environments/2025_06_12_environment.RData")

# libraries:
librarian::shelf(rjags, coda)

# time series data:
tcg <- tcg_ts  # harvard forest time series data from making plots for hotspot means - will need to update code for gee data

# make the disturbance magnitude version:
start <- grep("^2013", names(tcg))
end <- grep("^2015", names(tcg))
steady_state <- apply(tcg[,start:end], 1, mean)
dist <- grep("^2017", names(tcg))
post_dist <- tcg[,(dist + 1):ncol(tcg)]
dm_post_dist <- steady_state - post_dist

recov_state_space <- "model {
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
  x[s,1] ~ dnorm(x_ic, t_ic)
}

### Random Effects: (add later)
for (s in sites){
  asite[s] ~ dnorm(0, tausite)
}

for (t in 1:(nt-1)){
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

# data for model:
recov_sample <- as.matrix(dm_post_dist[sort(sample(1:nrow(tcg), 10)),])
# time series length:
time = 1:ncol(recov_sample)
sites = 1:nrow(recov_sample)
# first x:
x1 <- mean(tcg[,grep("^2017",names(tcg))])

# make list object
model_data <- list(y = recov_sample,
                   nt = length(time),
                   #time = rep(time, length(sites)),
                   sites = sites, 
                   t_obs = 0.001, a_obs = 0.001,
                   t_add = 0.001, a_add = 0.001,
                   r_ic = 1, r_prec = 0.001,
                   x_ic = x1, t_ic = 0.01)

# attempt model run:
jags_model <- jags.model(file = textConnection(recov_state_space),
                         data = model_data, 
                         n.chains = 3)
run_model <- coda.samples(jags_model,
                          variable.names = c("x[2,1]", "R[2,1]", "atime", "asite",
                                             "tau_obs", "tau_add", "tausite", "tautime",
                                             "r0"),
                          n.iter = 5000)

# let's see if it worked:

vars <- varnames(run_model)[c(1:3, 61:63, 100:103)]
params <- run_model[,vars]
plot(params)

plot(run_model)
out <- as.matrix(run_model)
test_site <- out[,grep("^x\\[1,", colnames(out))]
