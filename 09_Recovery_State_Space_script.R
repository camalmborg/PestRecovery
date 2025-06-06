# recovery rate state space model script
# this script contains the code for running a jags model for estimating recovery rates

# libraries:
librarian::shelf(rjags, coda)

# time series data:
tcg <- tcg_ts  # harvard forest time series data from making plots for hotspot means - will need to update code for gee data
# time series length:
time = 1:ncol(tcg)
nt = length(time)

recov_state_space <- "
model {
for (s in sites){

### Data Model:
  for (t in 1:n){
    y[t,s] ~ dnorm(x[t,s], tau_obs
  }

### Process Model:
for (t in 2:n){
    R[t,s] <- r0 #+ betar[z]  # here is where covariates of R would go? + alpha site + alpha t
    mu[t,s] <- R[t,s] * x[t-1,s]  # logit here?
    x[t,s] ~ dnorm(mu[t,s], tau_add)
  }
  
}

### Random Effects: (add later)
# for (s in sites){
#   asite ~ dnorm(0, tausite)
# }
# 
# for (t in 1:n){
#   atime ~ dnorm(0, tautime)
# }

### Priors:
x[1] ~ dnorm(x_ic, t_ic)
r0 ~ dnorm(r_ic, rprec)  # initial condition r
tau_obs ~ dgamma(t_obs, a_obs)
tau_add ~ dgamma(t_add, a_add)
#tausite ~ dgamma(0.001, 0.001)
#tautime ~ dgamma(0.001, 0.001)
"
