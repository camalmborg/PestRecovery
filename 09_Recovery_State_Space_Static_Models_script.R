# recovery rate state space model script for base and random effect models
# this script contains the code for running a jags model for estimating recovery rates
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)

# load environment if needed:
#load("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Environments/2025_06_12_environment.RData")
#load("Environments/2025_06_29_environment.RData")
load("Environments/2025_07_03_environment.RData")

# libraries:
librarian::shelf(rjags, coda)

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
                   dmags_tcg_y1 = dist_hist_tcg$`2016-05-01`,
                   dmags_tcg_y2 = dist_hist_tcg$`2017-05-01`,
                   dmag_tcg_sum = dist_hist_tcg$dist_2_yrs,
                   land_cov = as.factor(nlcd$category),
                   pct_tree_cov = nlcd$percent_tree_cover,
                   precip_2014 = static_daym$prcp[which(static_daym$year == 2014)],
                   precip_2015 = static_daym$prcp[which(static_daym$year == 2015)],
                   temp_mx_2014 = static_daym$tmax[which(static_daym$year == 2014)],
                   temp_mx_2015 = static_daym$tmax[which(static_daym$year == 2015)],
                   temp_min_2014 = static_daym$tmin[which(static_daym$year == 2014)],
                   temp_min_2015 = static_daym$tmin[which(static_daym$year == 2015)],
                   vpd_2014 = static_daym$vp[which(static_daym$year == 2014)],
                   vpd_2015 = static_daym$vp[which(static_daym$year == 2015)]) %>%
  mutate(precip_mean = rowMeans(select(., precip_2014, precip_2015))) %>%
  mutate(temp_mx_mean = rowMeans(select(., temp_mx_2014, temp_mx_2015))) %>%
  mutate(temp_min_mean = rowMeans(select(., temp_min_2014, temp_min_2015))) %>%
  mutate(vpd_mean = rowMeans(select(., vpd_2014, vpd_2015)))

# make list object
model_data <- list(y = recov_data,
                   nt = length(time),
                   sites = sites, 
                   t_obs = 0.001, a_obs = 0.001,
                   t_add = 0.001, a_add = 0.001,
                   r_ic = 1, r_prec = 0.001,
                   x_ic = x1, t_ic = 0.01)

recov_state_space_uni_static <- "model {
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