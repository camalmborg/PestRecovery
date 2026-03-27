### Script for sensitivities to drivers analyses ###

## Load libraries
library(dplyr)
library(tidyverse)
library(rjags)
library(coda)
library(ggplot2)
library(patchwork)
library(sf)
library(tigris)

## Loading TCG and TCG baselines:
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


## Loading parameter values from best models:
# load model performance information:
dic_sort <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/2025_11_30_all_recov_models_dics.csv", row.names = 1)

## Load Best Model
# all model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]
# best model:
best <- dic_sort$model_number[1]
best_model <- models[best] 
name = best_model
# model name:
model_name <- best_model
if (grepl("cov", name) == TRUE){
  model_name <- str_match(name, "cov_(.*?)_data")[,2]
} else {
  model_name <- str_match(name, "_(.*?)_model")[,2]
}
# load model:
load(paste0(dir, "model_runs/", best_model))
# load model_params:
model_params <- read.csv(file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/2025_11_30_all_recov_models_param_means.csv")

# collect site and time random effect model params:
jags_out <- model_info$jags_out
vars <- varnames(jags_out)
site_params <- jags_out[,grep("^as", vars)]
time_params <- jags_out[,grep("^at", vars)]
# remove burn in:
burn_in = 25000
site_params_burn <- window(site_params, start = burn_in)
time_params_burn <- window(time_params, start = burn_in)
out_site <- as.matrix(site_params_burn)
out_time <- as.matrix(time_params_burn)
# space random effects:
space_re <- apply(out_site, 2, mean, na.rm = T)
# make into spatial layer for map later:
coords <- data.frame(lon = tcg$longitude,
                     lat = tcg$latitude)
space_re_map_data <- sp::SpatialPointsDataFrame(coords, data = as.data.frame(space_re))


## Get Model Information and Data:
# get parameter values:
best_params <- model_params[best,]

# beta parameters:
betas <- as.numeric(best_params[grep("beta", names(best_params))])
betas <- betas[-which(is.na(betas))]
# r0 intercept:
r0 <- as.numeric(best_params[grep("r0", names(best_params))])
# taus:
tau_add <- sqrt(1/as.numeric(best_params[grep("tau_add", names(best_params))]))

# get model inputs:
model_in <- model_info$metadata$model_data
# get X starting:
start <- mean(tcg_base$`2017-05-01`, na.rm = T)
# standard deviation from start:
start_sd <- sd(tcg_base$`2017-05-01`, na.rm = T)
# mean baseline value:
base_mean <- mean(tcg_base$baseline, na.rm = T)

## Run Sensitivity Analyses for Parameters:
# sd's from mean:
sd_from_mean <- c(-2, -1, 0, 1, 2)

# sensitivity analyses list for covariates:
sens_list <- list()
for (i in 1:3){
  # disturbance magnitude correction case:
  if(i == 3){
    dm = 1
  } else {
    dm = 0
  }
  # set which Xo == sd from mean:
  Xo <- c(rep(0,3))
  
  # matrix for collecting one run:
  sens_mx <- matrix(NA, nrow = 5, ncol = 7)
  for(s in 1:5){
    # set up for filling in Xo:
    Xo[i] <- sd_from_mean[s]
    # set up X for filling in time:
    X <- c(rep(0, 7))
    X[1] <- start + dm*sd_from_mean[s]*start_sd
    for(t in 2:7){
      R <- r0 + betas[1]*Xo[1] + betas[2]*Xo[2] + betas[3]*(Xo[3]) #+ time_re[t-1]*Xo[4] + space_re[s]*Xo[5] + tau_add*Xo[6]
      X[t] <- R*X[t-1]
    }
    sens_mx[s,] <- base_mean - X
  }
  sens_list[[i]] <- as.data.frame(sens_mx)
}

## Run sensitivity analyses for random effects and error terms (with MCMC)
# loop:
for (i in 4:6){
  # set which Xo == sd from mean:
  Xo <- c(rep(0,6))
  Xo[i] = 1
  # matrix for collecting one run:
  sens_mx <- matrix(NA, nrow = 5000, ncol = 7)
  
    # run MCMC
    for(s in 1:5000){ 
      # set up X for filling in time:
      X <- c(rep(0, 7))
      X[1] <- start
      # sample time re params:
      alpha_t = out_time[sample.int(nrow(out_time),1),]
      
      for(t in 2:7){
        epsilon = rnorm(1 , 0, tau_add)
        R <- r0 + betas[1]*Xo[1] + betas[2]*Xo[2] + betas[3]*(Xo[3]) + alpha_t[t-1]*Xo[4] + space_re[s]*Xo[5] + epsilon*Xo[6]
        X[t] <- R*X[t-1]
      }
      sens_mx[s,] <- base_mean - X
    }
  ci <- apply(sens_mx, 2, quantile, pnorm(sd_from_mean))
  sens_list[[i]] <- as.data.frame(ci)
}


## Making Plots:
covs <- c("VPD", "Mean Max Temp 2015", "Disturbance Magnitude", 
          "Year Random Effect", "Spatial Random Effect", "Process Error")
years <- c(2017:2023)

#'@param sens_list = sensitivity analyses from loop
#'@param param_num = numeric; list member number
sens_plot_fx <- function(sens_list, param_num){
  # make plot data:
  plot_data <- sens_list[[param_num]] |>
    # rename columns with years:
    setNames(as.character(years)) |>
    # add column for covariate:
    mutate(sfm = sd_from_mean, .before = 1) |>
    # pivot for making plot data:
    pivot_longer(cols = -sfm,
                 names_to = "year",
                 values_to = "value") |>
    mutate(year = as.numeric(year))
  
  # make plot:
  sens_plot <- ggplot(data = plot_data, aes(x = year, y = value, group = sfm, color = sfm)) +
    geom_line() +
    scale_color_gradient(low = "blue", high = "red",
                         guide = guide_colorbar(direction = "vertical", title.position = "top")) +
    labs(title = covs[[param_num]],
         x = "Year",
         y = "TCG",
         color = "SD from\nMean") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.position = c(0.9, 0.25))#,
          # legend.title = element_text(size = 12, margin = margin(b = 13))) 
  return(sens_plot)
}

## Combining plots into single figure
# covariate plots:
VPD_plot <- sens_plot_fx(sens_list, 1) + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
temp_plot <- sens_plot_fx(sens_list, 2) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
distmag_plot <- sens_plot_fx(sens_list, 3) + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
year_re_plot <- sens_plot_fx(sens_list, 4) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
space_re_plot <- sens_plot_fx(sens_list, 5)
tau_add_plot <- sens_plot_fx(sens_list, 6) + theme(axis.title.y = element_blank())

# combine plots:
combined <- (VPD_plot + temp_plot + distmag_plot + year_re_plot + space_re_plot + tau_add_plot) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "right")
combined

# save:
# plot save location:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/"
png(filename = paste0(save_dir, Sys.Date(), "_sensitivities_plots_combined.png"),
    width = 10, height = 12, units = "in", res = 600)
combined
dev.off()


### Map of spatial random effects ###
# get states for mapping:
ma_ct_ri <- tigris::states(cb = TRUE) %>%
  filter(NAME %in% c("Massachusetts", "Connecticut", "Rhode Island"))

# convert spatial data points to points:
space_re_map <- st_as_sf(space_re_map_data)
# set crs to be same as states:
space_re_map <- st_set_crs(space_re_map, st_crs(ma_ct_ri))

# make the map:
spatial_re_map <- ggplot(space_re_map) +
  geom_sf(data = ma_ct_ri, fill = "grey", color = "black", size = 0.5) +
  # add the points:
  geom_sf(aes(fill = space_re), size = 1.5,
          color = "black", shape = 21, stroke = 0.15) +
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Magnitude of Spatial Random Effect",
       fill = "Spatial RE\nValue") +
  theme_bw() +
  theme(panel.grid = element_line(linewidth = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
#spatial_re_map

#save the map:
# Save the plot to a PNG file:
ggsave(paste0(save_dir, Sys.Date(), "_spatial_random_effect_map.png"),
       plot = spatial_re_map,
       width = 10, height = 6,
       dpi = 600)


### Time series of temporal random effect ###
# time parameters:
time_re <- apply(out_time, 2, mean, na.rm = T)
time_re_upper <- apply(out_time, 2, quantile, c(0.025))
time_re_lower <- apply(out_time, 2, quantile, c(0.925))
# make data for time series plot later:
time_re_time_series <- data.frame(year = as.character(c(2018:2023)),
                                  time_re = time_re, 
                                  upper = time_re_upper,
                                  lower = time_re_lower)

time_re_plot <- ggplot(data = time_re_time_series, mapping = aes(x = year, y = time_re, group = 1))+
  geom_point(size = 2.5, color = "navy") +
  geom_line(color = "navy", linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "navy", alpha = 0.15) +
  labs(x = "Year", y = "Temporal Random Effect") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
time_re_plot
# save the plot:
ggsave(paste0(save_dir, Sys.Date(), "_temporal_random_effect_plot.png"),
       plot = time_re_plot,
       width = 7, height = 3,
       dpi = 600)




### Archive ###

# # quickly attach to lat/lon for making a map later:
# space_re_map_data <- data.frame(lat = tcg$latitude, 
#                                 lon = tcg$longitude,
#                                 space_re = space_re)


# test <- sens_list[[3]]
# plot(test[1,], type = "l", ylim = c(min(test - 0.5), max(test)))
# for(i in 2:nrow(test)){
#   lines(test[i,], col = i)
# }

# test <- sens_list[[1]] |>
#   # rename columns with years:
#   setNames(as.character(years)) |>
#   # add column for covariate:
#   mutate(sfm = sd_from_mean, .before = 1) |>
#   # pivot for making plot data:
#   pivot_longer(cols = -sfm,
#                names_to = "year",
#                values_to = "value")

