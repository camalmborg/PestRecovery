### Script for testing spatial and temporal autocorrelation in best model residuals ###

## Load libraries
library(dplyr)
library(tidyverse)
library(rjags)
library(sf)
library(spatial)  
library(ggplot2)
library(tigris)

# load environment if needed:
load("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Environments/2025_07_07_environment.RData")

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)
# load model performance information:
dic_sort <- read.csv("2025_11_30_all_recov_models_dics.csv", row.names = 1)

## Load Best Model
# all model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]
# best model:
best <- dic_sort$model_number[1]
best_model <- models[best] 
# load model:
load(paste0(dir, "model_runs/", best_model))
# load model_params:
model_params <- read.csv(file = "2025_11_30_all_recov_models_param_means.csv")

## Get Model Information and Data
# get parameter values:
best_params <- model_params[best,]
name = best_model
# model name:
model_name <- best_model
if (grepl("cov", name) == TRUE){
  model_name <- str_match(name, "cov_(.*?)_data")[,2]
} else {
  model_name <- str_match(name, "_(.*?)_model")[,2]
}
# load model inputs:
model_inputs <- model_info$metadata$model_data

# get baselines:
tcg_base <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  # rename:
  rename_with(~ str_replace_all(.x, c("^\\s*X" = "", "\\." = "-"))) %>%
  # get baseline for anomolies:
  mutate(baseline = rowMeans(select(., `2010-05-01`:`2015-05-01`), na.rm = TRUE), .before = 1) %>%
  # create anomalies from baseline:
  mutate(across(!baseline, ~ baseline - .x))

tcg <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  # rename:
  rename_with(~ str_replace_all(.x, c("^\\s*X" = "", "\\." = "-"))) #%>%

  
## Load forecast result
forecast <- read.csv("Recovery_Forecasts/2025-12-03_ens_1500_model_1_start_year_2017_reforecast_result.csv")
# prepare for residual calculation:
y_pred <- forecast %>%
  # add baseline to predictions:
  mutate(across(!site, ~ tcg_base$baseline - .x)) %>%
  # ensemble means:
  mutate(site = as.factor(site)) %>%
  group_by(site) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  select(-c(site, V1))
  

## Residual calculation
# get the model input for y and covariates:
y <- tcg %>% select(`2018-05-01`:`2023-05-01`)
resid <- y_pred - y
colnames(resid) <- paste0("year_", seq(1, ncol(resid), length.out = ncol(resid)))


### Spatial Autocorrelation in Model Residuals ###
# Prepare a spatial data set for variogram
# add coords:
resid$lat <- coords$lat
resid$lon <- coords$lon
resid_spatial <- sp::SpatialPointsDataFrame(coords, data = resid)

# get states for mapping:
ma_ct_ri <- tigris::states(cb = TRUE) %>%
  filter(NAME %in% c("Massachusetts", "Connecticut", "Rhode Island"))
# convert spatial data points to points:
resid_sf <- st_as_sf(resid_spatial)
# set crs to be same as states:
resid_sf <- st_set_crs(resid_sf, st_crs(ma_ct_ri))

## Make maps for each year
mapping_residuals <- function(resid_col, resid_sf){
  # select chosen column in residuals data frame:
  column <- paste0("year_", resid_col)
  map_data <- c(column, "lon", "lat", "geometry")
  # grab map data:
  resid_map <- resid_sf[,map_data]
  colnames(resid_map) <- c("resids", map_data[2:4])
  resid_map <- resid_map %>%
    filter(!is.na(resids))
  
  # make the map:
  residual_map <- ggplot(resid_map) +
    geom_sf(data = ma_ct_ri, fill = "grey", color = "black", size = 0.5) +
    # add the points:
    geom_sf(aes(fill = resids), size = 1.5,
            color = "black", shape = 21, stroke = 0.1) +
    scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red", midpoint = 0) +
    theme_bw() +
    labs(title = paste0("Residuals: Predicted - Observed TCG, Year ", as.character(resid_col)),
      fill = "Residual") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  
  #save the map:
  # plot save location:
  save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/Residuals/"
  # Save the plot to a PNG file:
  ggsave(paste0(save_dir, Sys.Date(), "_", "Year_", as.character(resid_col), "_residuals_map.png"),
         plot = residual_map,
         width = 10, height = 6,
         dpi = 600)
}


## Make variograms for each year
spatial_ac <- function(resid, resid_col){
  # make trend surface:
  surf <- surf.ls(0, resid$lon, resid$lat, na.omit(resid[,resid_col]))
  # project matrix over region:
  # tr <- trmat(surf,
  #             min(resid$lon) , max(resid$lon),
  #             min(resid$lat), max(resid$lat),
  #             50) # 50x50m matrix
  vg <- spatial::variogram(surf, 100, plotit = F)
  cg <- spatial::correlogram(surf, 1000, xlim = c(0,1), plotit = F)
  
  # plot variogram:
  vg_data <- data.frame(xp = vg$x, yp = vg$y)
  vg_plot <- ggplot(data = vg_data, aes(x = xp, y = yp)) +
    geom_point() +
    labs(title = paste0("Variogram: Year ", as.character(resid_col))) +
    theme_bw() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  #vg_plot
  # plot save location:
  save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/Residuals/"
  # Save the plot to a PNG file:
  ggsave(paste0(save_dir, Sys.Date(), "_", "Year_", as.character(resid_col), "_variogram_plot.png"),
         plot = vg_plot,
         width = 10, height = 6,
         dpi = 600)
  
  # plot correlogram:
  cg_data <- data.frame(xp = cg$x, yp = cg$y)
  cg_plot <- ggplot(data = cg_data, aes(x = xp, y = yp)) +
    geom_point() +
    labs(title = paste0("Correlogram: Year ", as.character(resid_col))) +
    theme_bw()
  #cg_plot
  # plot save location:
  save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/Residuals/"
  # Save the plot to a PNG file:
  ggsave(paste0(save_dir, Sys.Date(), "_", "Year_", as.character(resid_col), "_correlogram_plot.png"),
         plot = cg_plot,
         width = 10, height = 6,
         dpi = 600)
}

# Run for each year:
for (i in 1:6){
  mapping_residuals(resid_col = i, resid_sf = resid_sf)
  spatial_ac(resid = resid, resid_col = i)
}


### Temporal Autocorrelation in model residuals ###
# resid = matrix of residuals by time
# create two vectors for covariance in repeated measures:
# vector that is timepoint 1 across sites, timepoint 2 across sites, through k-1 )(cols 1:5):
t_one <- as.vector(c(resid[,1], resid[,2], resid[,3], resid[,4], resid[,5]))
# vector that is cols 2:6:
t_minus_one <- as.vector(c(resid[,2], resid[,3], resid[,4], resid[,5], resid[,6]))
# ask the correlation between these:
lag_one_cor <- round(cor(t_one, t_minus_one, use = "complete.obs"), 4)
# do stats::cor.test for significance:
cor_test <- stats::cor.test(t_one, t_minus_one, method = "pearson")
cor_lag_one_test <- cor_test$p.value
plot(t_one, t_minus_one)

# do lag 2:
# vector that is cols 1:4:
t_one <- as.vector(c(resid[,1], resid[,2], resid[,3], resid[,4]))
# vector that is cols 3:6:
t_minus_one <- as.vector(c(resid[,3], resid[,4], resid[,5], resid[,6]))
# ask the correlation between these:
lag_two_cor <- round(cor(t_one, t_minus_one, use = "complete.obs"), 4)
# do stats::cor.test for significance:
cor_test <- stats::cor.test(t_one, t_minus_one, method = "pearson")
cor_lag_two_test <- cor_test$p.value
plot(t_one, t_minus_one)

# do lag 3:
# vector that is cols 1:3:
t_one <- as.vector(c(resid[,1], resid[,2], resid[,3]))
# vector that is cols 4:6:
t_minus_one <- as.vector(c(resid[,4], resid[,5], resid[,6]))
# ask the correlation between these:
lag_three_cor <- round(cor(t_one, t_minus_one, use = "complete.obs"), 4)
# do stats::cor.test for significance:
cor_test <- stats::cor.test(t_one, t_minus_one, method = "pearson")
cor_lag_three_test <- cor_test$p.value
plot(t_one, t_minus_one)

# do lag 4:
# vector that is cols 1:3:
t_one <- as.vector(c(resid[,1], resid[,2]))
# vector that is cols 4:6:
t_minus_one <- as.vector(c(resid[,5], resid[,6]))
# ask the correlation between these:
lag_four_cor <- round(cor(t_one, t_minus_one, use = "complete.obs"), 4)
# do stats::cor.test for significance:
cor_test <- stats::cor.test(t_one, t_minus_one, method = "pearson")
cor_lag_four_test <- cor_test$p.value
plot(t_one, t_minus_one)

# do lag 5:
# vector that is cols 1:3:
t_one <- as.vector(c(resid[,1]))
# vector that is cols 4:6:
t_minus_one <- as.vector(c(resid[,6]))
# ask the correlation between these:
lag_five_cor <- round(cor(t_one, t_minus_one, use = "complete.obs"), 4)
# do stats::cor.test for significance:
cor_test <- stats::cor.test(t_one, t_minus_one, method = "pearson")
cor_lag_five_test <- cor_test$p.value
plot(t_one, t_minus_one)

# make a table:
cor_lags <- data.frame(lag = c(1:5),
                       cor = c(lag_one_cor, lag_two_cor, lag_three_cor, lag_four_cor, lag_five_cor),
                       p_val = c(cor_lag_one_test, cor_lag_two_test, cor_lag_three_test, cor_lag_four_test, cor_lag_five_test))

plot(cor_lags$lag, cor_lags$cor)

### ARCHIVE ###

# # make trend surface:
# surf <- surf.ls(0, resid$lon, resid$lat, na.omit(resid[,1]))  # 4 had lowest AIC
# summary(surf)
# # project a matrix over the region:
# tr <- trmat(surf,
#             min(resid_spatial$lon) , max(resid_spatial$lon),
#             min(resid_spatial$lat), max(resid_spatial$lat),
#             50) # 50x50m matrix
# image(tr, asp = 5/5)
# 
# vg <- spatial::variogram(surf, 100)
# cg <- spatial::correlogram(surf, 1000, xlim = c(0,1))

# repeat for each year

# # extract:
# model_x <- params_out[,x_params]
# x_ci <- apply(model_x, 2, quantile, c(0.05, 0.5, 0.95))
# x_cols <- colnames(x_ci)
# x_means <- matrix(NA, nrow = 5000, ncol = 6)
# for (i in 1:nrow(x_means)){
#   # get correct columns for each row:
#   cols <- which(!is.na(str_match(x_cols, paste0(paste0("x\\[", i, ",")))))
#   # get values:
#   col_values <- x_ci[2, cols]
#   # fill in table:
#   x_means[i, 1:6] <- col_values[1:6]
# }

# ## Extract Model Parameters:
# # model output:
# jags_out <- model_info$jags_out
# # model parameter names:
# jags_vars <- varnames(jags_out)
# params <- jags_out[,grep("x|r|R|^tau|^b|^at", jags_vars)]
# # remove burn in:
# burn_in = 50000
# params_burn <- window(params, start = burn_in)
# # convert output to dataframe:
# params_out <- as.data.frame(as.matrix(params_burn))
# # separate out specific params:
# x_params <- grep("^x", colnames(params_out))
# beta_params <- grep("^b",colnames(params_out))
# taus <- grep("tau", colnames(params_out))
# r_ints <- grep("r0", colnames(params_out))
# r_rates <- grep("R", colnames(params_out))
# # beta parameters:
# beta_one <- best_params$beta.1.
# beta_two <- best_params$beta.2.
# beta_three <- best_params$beta.3.
# # tau parameters (convert sd = 1/precision):
# tau_obs <- 1/best_params$tau_obs
# tau_add <- 1/best_params$tau_add
# # time random effect parameters:
# a_times <- c(best_params$atime.1., best_params$atime.2., best_params$atime.3.,
#              best_params$atime.4., best_params$atime.5., best_params$atime.6.)
# # recovery rate intercept:
# r_int <- best_params$r0
# # starting TCG:
# x_init <- model_info$jags_out
# 
# # covariates:
# cov_one <- as.matrix(model_inputs$cov_one)
# cov_two <- as.matrix(model_inputs$cov_two)
# 
# ## Prepare to run forward in time:
# y_pred <- matrix(NA, nrow = 5000, ncol = 6)
# for (i in 1:nrow(y_pred)){
#   for(j in 1:ncol(y_pred)){
#     R <- r_int + a_times[j] + (beta_one*cov_one[i, j]) + (beta_two*cov_two[i])
#     mu <- R * x_init[i,j]
#     y_pred[i,j] <- mu
#   }
# }
