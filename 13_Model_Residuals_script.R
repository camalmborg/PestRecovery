### Script for testing spatial and temporal autocorrelation in best model residuals ###

## Load libraries
library(dplyr)
library(tidyverse)
library(rjags)
library(spatial)  
library(maps) 
library(sp)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

## Load Best Model
# all model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]
# best model:
best <- 38  # also important: 42, 29, 7
best_model <- models[best] 
# load model:
load(paste0(dir, "model_runs/", best_model))

## Get Model Information and Data
# get parameter values:
best_model_params <- model_params[best,]
# model name:
model_name <- best_model
if (grepl("cov", name) == TRUE){
  model_name <- str_match(name, "cov_(.*?)_data")[,2]
} else {
  model_name <- str_match(name, "_(.*?)_model")[,2]
}
# load model inputs:
model_inputs <- model_info$metadata$model_data

## Extract Model Parameters:
# model output:
jags_out <- model_info$jags_out
# model parameter names:
jags_vars <- varnames(jags_out)
params <- jags_out[,grep("r|R|^tau|^b|^at", jags_vars)]
# remove burn in:
burn_in = 25000
params_burn <- window(params, start = burn_in)
# convert output to dataframe:
params_out <- as.data.frame(as.matrix(params_burn))
# separate out specific params:
x_params <- grep("^x", colnames(out))
beta_params <- grep("^b",colnames(out))
taus <- grep("tau", colnames(out))
r_ints <- grep("r0", colnames(out))
r_rates <- grep("R", colnames(out))
# extract:
model_x <- out[,x_params]
x_ci <- apply(model_x, 2, quantile, c(0.05, 0.5, 0.95))
x_cols <- colnames(x_ci)
x_means <- matrix(NA, nrow = 5000, ncol = 7)
for (i in 1:nrow(x_means)){
  # get correct columns for each row:
  cols <- which(!is.na(str_match(x_cols, paste0(paste0("x\\[", i, ",")))))
  # get values:
  col_values <- x_ci[2, cols]
  # fill in table:
  x_means[i, 1:7] <- col_values[1:7]
}
# beta parameters:
best_params <- model_params[best,]
beta_one <- best_params$`beta[1]`
beta_two <- best_params$`beta[2]`
# tau parameters (convert sd = 1/precision):
tau_obs <- 1/best_params$tau_obs
tau_add <- 1/best_params$tau_add
# time random effect parameters:
a_times <- c(best_params$`atime[1]`, best_params$`atime[2]`, best_params$`atime[3]`,
            best_params$`atime[4]`, best_params$`atime[5]`, best_params$`atime[6]`)
# recovery rate intercept:
r_int <- best_params$r0
# starting TCG:
x_init <- tcg %>%
  select(`2017-05-01`:`2023-05-01`)
# covariates:
cov_one <- model_inputs$cov_one
cov_two <- model_inputs$cov_two

## Prepare to run forward in time:
y_pred <- matrix(NA, nrow = 5000, ncol = 6)
for (i in 1:nrow(y_pred)){
  for(j in 1:ncol(y_pred)){
    R[i,j] <- r_int + a_times[j] + (beta_one*cov_one[i, j]) + (beta_two*cov_two[i])
    mu[i,j] <- R[i,j] * x_init[i,j]
    y_pred[i,j] <- mu[i,j]
  }
}

## Residual calculation
# get the model input for y and covariates:
y <- as.data.frame(model_inputs$y) %>% select(-c(`2024-05-01`))
resid <- as.data.frame(y_pred) - y
colnames(resid) <- colnames(y)


### Spatial Autocorrelation in Model Residuals ###
# Prepare a spatial data set for variogram
# add coords:
resid$lat <- coords$lat
resid$lon <- coords$lon
resid_spatial <- sp::SpatialPointsDataFrame(coords, data = resid)

# Make maps for each year, color by resid values

# test plot:
plot(resid$lon, resid$lat, pch = 16)
map("state",add=TRUE)

## Make variograms for each year
# make trend surface:
surf <- surf.ls(4, resid$lon, resid$lat, na.omit(resid[,1]))  # 4 had lowest AIC
# project a matrix over the region:
tr <- trmat(surf, 
            min(resid_spatial$lon) , max(resid_spatial$lon),
            min(resid_spatial$lat), max(resid_spatial$lat), 
            50) # 50x50m matrix
image(tr, asp = 3/5) 

vg <- spatial::variogram(surf, 100)
cg <- spatial::correlogram(surf, 1000, xlim = c(0,1))

# repeat for each year


### Temporal Autocorrelation in model residuals ###
# have a matrix of residuals by time
# create 1 vector that is starting point in time and another vector that is ending point in time
# vector that is timepoint 1 across sites, timepoint 2 across sites, through k-1 )(cols 1:5)
# vector that is cols 2:6
# ask the correlation between these
# lag 2 is 1:4/3:6
# lag 3 is 1:3/4:6

