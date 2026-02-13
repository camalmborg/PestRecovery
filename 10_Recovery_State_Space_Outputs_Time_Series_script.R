### Sample Time Series plots ###

## Load libraries and necessary environments
librarian::shelf(dplyr, tidyverse, rjags, coda)
# load libraries for plots:
librarian::shelf(ggplot2, hrbrthemes)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

# load environment if needed:
load("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Environments/2025_07_07_environment.RData")

## pull in model output files if not in the environment already
# model files:
models <- list.files(paste0(dir, "model_runs/"))[grep("RData", list.files(paste0(dir, "model_runs/")))]
# best model:
dic_sort <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/2025_11_30_all_recov_models_dics.csv")

## Load Best Model
# all model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]
# best model:
best <- as.numeric(dic_sort$model_number[1])
best_model <- models[best] 
# load model:
load(paste0(dir, "model_runs/", best_model))
# load model_params:
#model_params <- read.csv(file = "2025_11_30_all_recov_models_param_means.csv")

# choose model:
#m_num <-  # change this when changing models
model_pick <- models[m_num]
# load model_info object
load(paste0(dir, "model_runs/", model_pick))
model_inputs <- model_info$metadata$model_data
# model parameters:
jags_out <- model_info$jags_out
# remove burn in:
burn_in = 50000
jags_burn <- window(jags_out, start = burn_in) 
# get parameters from model output:
out <- as.matrix(jags_burn)
# separate out specific params:
x_params <- grep("^x", colnames(out))

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

# forecast time period:
years <- 2017:2023
# observations:
tcg <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  # rename:
  rename_with(~ str_replace_all(.x, c("^\\s*X" = "", "\\." = "-"))) %>%
  # select years:
  select(`2017-05-01`:`2023-05-01`) %>%
  # rename with just years:
  rename_with(~ sub("-.*", "", .x))

## Load forecast result
forecast <- read.csv("Recovery_Forecasts/2025-12-03_ens_1500_model_1_start_year_2017_reforecast_result.csv")
# prepare for residual calculation:
y_pred <- forecast %>%
  # add baseline to predictions:
  mutate(across(-site, ~ tcg_base$baseline[site] - .x)) %>%
  # ensemble means:
  mutate(site = as.factor(site)) %>%
  group_by(site) %>%
  summarise(across(where(is.numeric), ~ quantile(.x, c(0.50), na.rm = TRUE)))
y_upper <- forecast %>%
  mutate(across(-site, ~ tcg_base$baseline[site] - .x)) %>%
  mutate(site = as.factor(site)) %>%
  group_by(site) %>%
  summarise(across(where(is.numeric), ~ quantile(.x, c(0.90), na.rm = TRUE))) 
y_lower <- forecast %>%
  mutate(across(-site, ~ tcg_base$baseline[site] - .x)) %>%
  mutate(site = as.factor(site)) %>%
  group_by(site) %>%
  summarise(across(where(is.numeric), ~ quantile(.x, c(0.10), na.rm = TRUE)))
y_up <- forecast %>%
  mutate(across(-site, ~ tcg_base$baseline[site] - .x)) %>%
  mutate(site = as.factor(site)) %>%
  group_by(site) %>%
  summarise(across(where(is.numeric), ~ quantile(.x, c(0.75), na.rm = TRUE))) 
y_low <- forecast %>%
  mutate(across(-site, ~ tcg_base$baseline[site] - .x)) %>%
  mutate(site = as.factor(site)) %>%
  group_by(site) %>%
  summarise(across(where(is.numeric), ~ quantile(.x, c(0.25), na.rm = TRUE)))
  


# making time series:
#a <- sample(3000:5000, 1)
#sample <- a
sample <- 3372

# from model outputs:
xs <- out[,x_params]
x_samp <- xs[,grep(paste0("x\\[", as.character(sample),","), colnames(xs))]
y_ci <- apply(x_samp, 2, quantile, c(0.10, 0.5, 0.90))
y_ci <- tcg_base$baseline[sample] - y_ci

# from forecast outputs:
f_ci <- rbind(y_upper[sample,], y_pred[sample,], y_lower[sample,], y_up[sample,], y_low[sample,])

# observation:
obs <- tcg[sample,]

## Making Time Series
# prepare model data:
plot_data <- data.frame(date = as.numeric(names(obs)),
                        obs = as.numeric(obs),
                        y_low = as.numeric(y_ci[1,]),
                        y_med = as.numeric(y_ci[2,]),
                        y_high = as.numeric(y_ci[3,]),
                        x_low = as.numeric(f_ci[3, -1]),
                        x_med = as.numeric(f_ci[2, -1]),
                        x_high = as.numeric(f_ci[1, -1]),
                        x_low2 = as.numeric(f_ci[4, -1]),
                        x_high2 = as.numeric(f_ci[5, -1]))

plot_name <- sub(".*multi_(.*?)_data.*", "\\1", model_pick)

# make the plot layering observation and model preds:
time_series <- ggplot(data = plot_data) +
  # time series for observations:
  geom_point(aes(x = date, y = obs, color = "Observations"), size = 2.5) +
  geom_line(aes(x = date, y = obs, color = "Observations"), linetype = "solid") +
  # time series for model:
  geom_point(aes(x = date, y = y_med,
             color = "Model"), size = 2) +
  geom_line(aes(x = date, y = y_med,
            color = "Model"), linetype = "dashed") +
  # add confidence intervals:
  geom_ribbon(aes(x = date, ymin = y_low, ymax = y_high,
              fill = "Model 90% Interval"), alpha = 0.25) +
  # add base model for compare:
  geom_point(aes(x = date, y = x_med,
             color = "Forecast"), size = 2) +
  geom_line(aes(x = date, y = x_med,
            color = "Forecast"), linetype = "dotdash", linewidth = 0.5) +
  geom_ribbon(aes(x = date, ymin = x_low, ymax = x_high,
              fill = "Forecast 90% Interval"), alpha = 0.10) +
  geom_ribbon(aes(x = date, ymin = x_low2, ymax = x_high2,
              fill = "Forecast 75% Interval"), alpha = 0.20) +
  scale_x_continuous(breaks = sort(unique(plot_data$date))) +
  # colors:
  scale_color_manual(name = "Lines",
                     breaks = c("Observations", "Model", "Forecast"),
                     values = c("Observations" = "black", 
                                "Model" = "red", 
                                "Forecast" = "navy")) +
  scale_fill_manual(name = "Confidence Intervals",
                    breaks = c("Model 90% Interval",
                               "Forecast 90% Interval",
                               "Forecast 75% Interval"),
                    values = c("Model 90% Interval" = "red",
                               "Forecast 90% Interval" = "navy",
                               "Forecast 75% Interval" = "navy")) +
  guides(color = guide_legend(order = 1),
         fill  = guide_legend(order = 2)) +
  # set the axis limits:
  #ylim(c(-1, 1)) +
  # seeing plot closer to time series:
  #coord_cartesian(ylim = c(min(obs, na.rm = T) - 0.005, max(obs, na.rm = T) + 0.01)) +
  # plot labels:
  labs(title = "Sample Time Series (single Landsat pixel)",
       y = "Tasseled Cap Greenness", 
       x = "Year") +
  #scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

time_series
#sample

# save it:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/"
setwd(save_dir)
# save:
png(filename = paste0(save_dir, "2026_02_09_time_series_sample_point.png"),
    height = 6,
    width = 10,
    units = "in",
    res = 600)
time_series
dev.off()
# ggsave("2026_02_09_time_series_sample_point.png", 
#        plot = time_series, 
#        width = 14, 
#        height = 8, 
#        dpi = 600)


# # playing around:
# sites <- sample(1:5000, size = 10, replace = FALSE)
# testing <- tcg[sites,] %>%
#   pivot_longer(cols = everything(), names_to = "year", values_to = "obs")
# 
# 
# test_plot <- ggplot(data = testing) +

# # test plot:
# plot(1:7, obs, ylim = c(min(y_ci[,-1]), max(y_ci[,-1])))
# lines(1:7, y_ci[2,-1])
# lines(1:7, y_ci[1,-1], col = "red")
# lines(1:7, y_ci[3,-1], col = "blue")
# #plot(as.Date(names(obs)), x_ci[2,])
  


# save them:
# save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/"
# setwd(save_dir)
# all vars:
# # Save the plot to a PNG file
# ggsave("2025_08_06_sample_time_series_with_base_ESA.png", 
#        plot = time_series,
#        height = 7,
#        width = 8,
#        units = "in",
#        dpi = 600)


