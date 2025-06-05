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
### Data Model:
for (t in 1:n){
  y[t] ~ dnorm(x[t], tau_obs
}

### Process Model:
for (t in 2:n){
  mu[t] <- R[s] * x[t-1]
  x[t] ~ dnorm(mu[t], tau_add)
}

for (s in 1:sites){
  R[s] <- 
}

