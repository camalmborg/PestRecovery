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

## Make sample time series for selected models:
best = dic_sort[which(dic_sort$perform == 1),]
m_num = best$model_number
best = models[m_num]

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
beta_params <- grep("^b",colnames(out))
taus <- grep("tau", colnames(out))
a_time <- grep("at", colnames(out))
a_site <- grep("as", colnames(out))
r <- grep("r0", colnames(out))

# ## Load forecasts
# # set working directory:
# dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/"
# setwd(dir)
# 
# # load model forecast files:
# files <- list.files(pattern = "result.csv$")
# #model_num = as.numeric(Sys.getenv("SGE_TASK_ID"))
# model_num = 1
# # get years of analysis:
# start_year <- as.numeric(stringr::str_extract(files[model_num], "(?<=start_year_)\\d+"))
# years <- start_year:2023
# # loading predicted forecast values file:
# model_out <- read.csv(files[model_num])
# pred <- model_out %>%
#   # rename columns with years:
#   rename_with(~ as.character(years)[seq_along(.)], .cols = -1)


years <- 2017:2023
# load observation data:
tcg_bs <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  select(starts_with("X")) %>%
  rename_with(~ str_replace_all(., c("X|_tcg_mean" = "", "\\." = "-"))) %>%
  # get baseline value 2010-2015:
  mutate(baseline = rowMeans(select(., `2010-05-01`:`2015-05-01`), na.rm = TRUE), .before = 1) %>%
  # create anomalies from baseline:
  mutate(across(!baseline, ~ baseline - .x)) %>%
  # add site and site and lat lon:
  mutate(site = 1:nrow(tcg), .before = 1) %>% 
  mutate(longitude = coords$lon, .before = 2) %>%
  mutate(latitude = coords$lat, .before = 3)

tcg_obs <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  select(starts_with("X")) %>%
  rename_with(~ str_replace_all(., c("X|_tcg_mean" = "", "\\." = "-"))) %>%
  rename_with(~ str_replace_all(., "-.*", "")) %>%
  # select 2017-2023:
  select(c(`2017`:`2023`))
  

# making time series:
a <- sample(1:5000, 1)
sample <- a

# from model outputs:
xs <- out[,x_params]
x_samp <- xs[,grep(paste0("x\\[", as.character(sample),","), colnames(xs))]
y_ci <- apply(x_samp, 2, quantile, c(0.05, 0.5, 0.95))

# baseline of sample:
obs_base <- tcg$baseline[sample]
y_ci_samp <- obs_base - y_ci

# observation:
obs <- tcg_obs[sample,]

# test plot:
plot(1:7, obs)
lines(1:7, y_ci_samp[2,])
lines(1:7, y_ci_samp[1,], col = "red")
lines(1:7, y_ci_samp[3,], col = "blue")
#plot(as.Date(names(obs)), x_ci[2,])

## Making Time Series
# prepare model data:
plot_data <- data.frame(date = as.numeric(names(obs)),
                        obs = as.numeric(obs),
                        y_low = as.numeric(y_ci[1,]),
                        y_med = as.numeric(y_ci[2,]),
                        y_high = as.numeric(y_ci[3,]))

plot_name <- sub(".*multi_(.*?)_data.*", "\\1", model_pick)

# make the plot layering observation and model preds:
time_series <- ggplot(data = plot_data) +
  # time series for observations:
  geom_point(aes(x = date, y = obs), color = "black", size = 2) +
  geom_line(aes(x = date, y = obs), color = "black", linetype = "solid") +
  # time series for model:
  geom_point(aes(x = date, y = y_med),
             color = "red", size = 2) +
  geom_line(aes(x = date, y = y_med),
            color = "red", linetype = "dashed") +
  # add confidence intervals:
  geom_ribbon(aes(x = date, ymin = y_low, ymax = y_high),
              fill = "red", alpha = 0.25) +
  # # add base model for compare:
  # geom_point(data = base_model_data, aes(x = date, y = x_med),
  #            color = "navy", size = 1) +
  # geom_line(data = base_model_data, aes(x = date, y = x_med),
  #           color = "navy", linetype = "dashed", linewidth = 0.5) +
  # geom_ribbon(data = base_model_data, aes(x = date, ymin = x_low, ymax = x_high),
  #             fill = "navy", alpha = 0.15) +
  # set the axis limits:
  #ylim(c(-1, 1)) +
  # seeing plot closer to time series:
  #coord_cartesian(ylim = c(min(obs, na.rm = T) - 0.005, max(obs, na.rm = T) + 0.01)) +
  # plot labels:
  labs(title = "Sample Time Series",
       y = "Tasseled Cap Greenness", 
       x = "Year") +
  #scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  theme_bw() +
  theme(legend.position = "right")

time_series


# # playing around:
# sites <- sample(1:5000, size = 10, replace = FALSE)
# testing <- tcg[sites,] %>%
#   pivot_longer(cols = everything(), names_to = "year", values_to = "obs")
# 
# 
# test_plot <- ggplot(data = testing) +
  


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


