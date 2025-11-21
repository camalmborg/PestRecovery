### Sample Time Series plots ###

## Load libraries and necessary environments
librarian::shelf(dplyr, tidyverse, rjags, coda)
# load libraries for plots:
librarian::shelf(ggplot2, hrbrthemes)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

## pull in model output files if not in the environment already
#model_params <- read.csv("2025_07_31_all_base_uni_recov_models_param_means.csv")
#load("2025_07_31_recov_models_outputs_list.RData")  # object is called model_outputs
# model files:
models <- list.files(paste0(dir, "model_runs/"))[grep("RData", list.files(paste0(dir, "model_runs/")))]
# best model:
dic_sort <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/2025_11_17_multi_recov_models_dics.csv")

## Make sample time series for selected models:
# base = dic_sort[which(dic_sort$covariate == "base"),]
# base = base$model_number
# base = models[base]
# base with time random effect = models[4]
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
out <- as.matrix(model_info$jags_out)
# separate out specific params:
x_params <- grep("^x", colnames(out))
beta_params <- grep("^b",colnames(out))
taus <- grep("tau", colnames(out))
r <- grep("r0", colnames(out))

## Load forecasts
# set working directory:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/"
setwd(dir)

# load model forecast files:
files <- list.files(pattern = "result.csv$")
#model_num = as.numeric(Sys.getenv("SGE_TASK_ID"))
model_num = 1
# get years of analysis:
start_year <- as.numeric(stringr::str_extract(files[model_num], "(?<=start_year_)\\d+"))
years <- start_year:2023
# loading predicted forecast values file:
model_out <- read.csv(files[model_num])
pred <- model_out %>%
  # rename columns with years:
  rename_with(~ as.character(years)[seq_along(.)], .cols = -1)

# load observation data:
tcg <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  # select columns with observations for 2017-2023:
  select(matches(as.character(years))) %>%
  # rename columns for years:
  rename_with(~ as.character(years)[seq_along(.)])


# making time series:
a <- sample(1:5000, 1)
sample <- a

# from model outputs:
xs <- out[,x_params]
x_samp <- xs[,grep(paste0("x\\[", as.character(sample),","), colnames(xs))]
y_ci <- apply(x_samp, 2, quantile, c(0.05, 0.5, 0.95))

# # sampling:
# test <- pred[pred$`2017` < 0.15,]

# a <- sample(unique(test$site), 1)

# #y <- as.matrix(model_inputs$y)
# y <- pred %>%
#   # select sample site:
#   filter(site == sample) %>%
#   # get only tcg values:
#   select(-c(site))
# y_ci <- apply(y, 2, quantile, c(0.05, 0.5, 0.95))

# observations:
# example sites: 77**, 346*, 4706*, 4911*, 4859**, 
# 4584*, 4466**, 4178**, 4367, 3754**, 3223**, 4473, 3818
obs <- tcg[sample,]
# test plot:
plot(as.Date(names(obs), format = "%Y"), obs)
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


