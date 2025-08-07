### Sample Time Series plots ###

## Load libraries and necessary environments
librarian::shelf(dplyr, tidyverse, rjags, coda)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

## pull in model output files if not in the environment already
model_params <- read.csv("2025_07_31_all_base_uni_recov_models_param_means.csv")
load("2025_07_31_recov_models_outputs_list.RData")  # object is called model_outputs
# model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]

## Make sample time series for selected models:
# base = models[1]
# base with time random effect = models[4]
# uni_vp = models[8]
# uni_dmags_tcg_y2 = models[11]

# choose model:
m_num <- 8  # change this when changing models
model_pick <- models[m_num]
# load model_info object
load(paste0(dir, "model_runs/", model_pick))
model_inputs <- model_info$metadata$model_data
# model parameters:
out <- as.matrix(model_info$jags_out)
# separate out specific params:
x_params <- grep("^x", colnames(out))
beta_params <- grep("^b",colnames(out))
taus <- grep("tau", colnames(out))
r <- grep("r0", colnames(out))

# making time series:
sample <- 161
y <- as.matrix(model_inputs$y)
obs <- y[sample,]  # example sites: 1, 6, 9, 52, 56, 77**, 87, 98, 104, 124, 125, 144, 156*, 161**, 1550, 1990
# get confidence intervals of x from model outputs:
xs <- out[,x_params]
x_samp <- xs[,grep(paste0("x\\[", as.character(sample),","), colnames(xs))]
x_ci <- apply(x_samp, 2, quantile, c(0.05, 0.5, 0.95))

# test plot:
#plot(as.Date(names(obs)), obs)
#plot(as.Date(names(obs)), x_ci[2,])

## Making pretty time series plots
# load libraries for plots:
librarian::shelf(ggplot2, hrbrthemes)

# prepare model data:
plot_data <- data.frame(date = as.Date(names(obs)),
                        obs = obs,
                        x_low = x_ci[1,],
                        x_med = x_ci[2,],
                        x_high = x_ci[3,])

# base_model_data <- data.frame(date = as.Date(names(obs)),
#                               obs = obs,
#                               x_low = x_ci[1,],
#                               x_med = x_ci[2,],
#                               x_high = x_ci[3,])

plot_name <- sub(".*uni_(.*?)_data.*", "\\1", model_pick)

# make the plot layering observation and model preds:
time_series <- ggplot(data = plot_data) +
  # time series for observations:
  geom_point(aes(x = date, y = obs), color = "black", size = 2) +
  geom_line(aes(x = date, y = obs), color = "black", linetype = "solid") +
  # time series for model:
  geom_point(aes(x = date, y = x_med),
             color = "red", size = 2) +
  geom_line(aes(x = date, y = x_med),
            color = "red", linetype = "dashed") +
  # add confidence intervals:
  geom_ribbon(aes(x = date, ymin = x_low, ymax = x_high),
              fill = "red", alpha = 0.25) +
  # add base model for compare:
  geom_point(data = base_model_data, aes(x = date, y = x_med),
             color = "navy", size = 1) +
  geom_line(data = base_model_data, aes(x = date, y = x_med),
            color = "navy", linetype = "dashed", linewidth = 0.5) +
  geom_ribbon(data = base_model_data, aes(x = date, ymin = x_low, ymax = x_high),
              fill = "navy", alpha = 0.15) +
  # set the axis limits:
  #ylim(c(-0.095, 0.0025)) +
  # plot labels:
  labs(title = "Sample Time Series",
       y = "Tasseled Cap Greenness", 
       x = "Year") +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  theme_bw() +
  theme(legend.position = "right")

time_series

# save them:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/"
setwd(save_dir)
# all vars:
# Save the plot to a PNG file
ggsave("2025_08_06_sample_time_series_with_base_ESA.png", 
       plot = time_series,
       dpi = 600)


