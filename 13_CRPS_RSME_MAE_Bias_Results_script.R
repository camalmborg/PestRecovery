### CRPS and Bias Results Script ###

## Load libraries
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(patchwork)
library(gt)

## Loading CRPS and Bias run results
# working directories:
crps_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/CRPS/"
rmse_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/RMSE_Bias/"
# files:
crps_files <- list.files(crps_dir)[grep("across_site", list.files(crps_dir))]
rmse_files <- list.files(rmse_dir)[grep("rmse_mae_bias", list.files(rmse_dir))]

## CRPS
# loading all as list:
crps_load <- lapply(paste0(crps_dir, crps_files), read.csv)
# identify the maximum col num:
max_col <- max(sapply(crps_load, ncol))
# colnames: 
columns <- c("Model", as.character(2017:2023))

# how to make it into one dataframe:
crps_all <- matrix(data = NA, nrow = length(crps_load), ncol = max_col)
model_names <- c()
colnames(crps_all) <- columns
# loop to fill in matrix:
for (i in 1:length(crps_load)){
  # select name:
  name <- crps_files[i]
  name_extract <- str_extract(name, "(?<=_model_).*?(?=_crps)")
  # add name to column
  model_names[i] <- name_extract
  # load member:
  crps <- crps_load[[i]] %>%
    rename(model = 1) %>%
    select(-model) %>%
    rename_with(~ str_replace_all(., c("X" = "")))
  # put it in the matrix:
  for (j in colnames(crps)){
    crps_all[i, j] <- as.numeric(crps[,j])
  }
}
crps_all <- as.data.frame(crps_all)
crps_all$Model <- model_names

## RMSE, MAE, and Bias
# load files:
bias_load <- lapply(paste0(rmse_dir, rmse_files), read.csv)
# identify the maximum col num:
max_col <- max(sapply(bias_load, ncol))
# colnames: 
columns <- c("Model", "Metric", as.character(2017:2023))
# metric:
metric <- c("RMSE", "MAE", "Bias")

# how to make it into one dataframe:
rmse_all <- matrix(data = NA, nrow = length(bias_load), ncol = max_col + 1)
mae_all <- matrix(data = NA, nrow = length(bias_load), ncol = max_col + 1)
bias_all <- matrix(data = NA, nrow = length(bias_load), ncol = max_col + 1)
colnames(rmse_all) <- columns
colnames(mae_all) <- columns
colnames(bias_all) <- columns
model_names <- c()

# loading one:
#i = 1
# loop to fill in matrices:
for (i in 1:length(bias_load)){
  # select name:
  name <- rmse_files[i]
  name_extract <- str_extract(name, "(?<=_model_).*?(?=_rmse)")
  # add name to column
  model_names[i] <- name_extract
  # separate all bias data:
  all_bias <- bias_load[[i]] %>%
    rename(Metric = 1) %>%
    rename_with(~ str_replace_all(., c("X" = ""))) %>%
    pivot_longer(
      cols = starts_with("2"),
      names_to = "Year",
      values_to = "Score")
  # select each group:
  rmse <- all_bias %>%
    filter(Metric == "RMSE") %>%
    pivot_wider(names_from = Year,
                values_from = Score) %>%
    select(-Metric)
  mae <- all_bias %>%
    filter(Metric == "MAE") %>%
    pivot_wider(names_from = Year,
                values_from = Score) %>%
    select(-Metric)
  bias <- all_bias %>%
    filter(Metric == "Bias") %>%
    pivot_wider(names_from = Year,
                values_from = Score) %>%
    select(-Metric)
  # put it in the matrix:
  for (j in colnames(rmse)){
    rmse_all[i, j] <- suppressWarnings(as.numeric(rmse[,j]))
    mae_all[i, j] <- suppressWarnings(as.numeric(mae[,j]))
    bias_all[i, j] <- suppressWarnings(as.numeric(bias[,j]))
  }
}
rmse_all <- as.data.frame(rmse_all)
mae_all <- as.data.frame(mae_all)
bias_all <- as.data.frame(bias_all)
rmse_all$Model <- model_names
rmse_all$Metric <- "RMSE"
mae_all$Model <- model_names
mae_all$Metric <- "MAE"
bias_all$Model <- model_names
bias_all$Metric <- "Bias"


## Preparing data for score vs lead time and score vs 1-year lag plots
# function for getting results:
#'@param test = crps/rmse/mae/bias data
#'@param starts = character object with years being forecast 
get_plot_data <- function(test, starts, name){
  result <- test %>%
    # add column for model number for indexing:
    mutate(model_num = if_else(str_detect(Model, "^base_"), "base", str_extract(Model, "\\d+")), .after = 1)
  # make empty matrix to fill with diagonals:
  plots <- matrix(NA, nrow = nrow(result), ncol = length(starts))
  for (i in 1:nrow(plots)){
    model_num <- result$model_num[i]
    get_years <- result[i, grep("2", colnames(result))]
    nas <- which(is.na(get_years)) 
    if (length(nas) > 0){
      casts <- result[result$model_num == model_num, starts[-nas]]
    } else {
      casts <- result[result$model_num == model_num, starts]
    }
    diag <- diag(as.matrix(casts))
    plots[i,1:length(diag)] <- diag
  }
  colnames(plots) <- starts
  # get diagonal means and one-year lags:
  plot_data <- as.data.frame(plots) %>%
    mutate(model_num = result$model_num, .before = 1) %>%
    mutate(diag_mean = rowMeans(select(., -1), na.rm = TRUE))
  # which have all years:
  nums <- which(complete.cases(plot_data))
  plot_data$yr_one_lag <- c(plots[nums[1],], plots[nums[2],], plots[nums[3],], plots[nums[4],])
  # years for each:
  n_years <- rep(1:length(starts), 4)
  cast_years <- rep(starts, 4)
  plot_data$year <- n_years
  plot_data$cast_year <- cast_years
  # make plot data:
  plot_data <- plot_data[,-grep("2", colnames(plot_data))]
  # add column for test type:
  plot_data$Test <- name
  return(plot_data)
}

# years of forecasts:
starts <- as.character(2018:2023)
# use function:
crps_plot_data <- get_plot_data(crps_all, starts, "CRPS")
rmse_plot_data <- get_plot_data(rmse_all, starts, "RMSE")
mae_plot_data <- get_plot_data(mae_all, starts, "MAE")
bias_plot_data <- get_plot_data(bias_all, starts, "Bias")

## Making Plots
# diagonal means (metric vs lead time):
crps_plot <- ggplot(crps_plot_data, aes(x = year, y = diag_mean, color = as.factor(model_num), group = as.factor(model_num))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:6) +
  labs(title = "CRPS",
       x = NULL,
       y = "CRPS Score",
       color = "Forecast Model") +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        axis.text.x = element_blank())

rmse_plot <- ggplot(rmse_plot_data, aes(x = year, y = diag_mean, color = as.factor(model_num), group = as.factor(model_num))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:6) +
  labs(title = "RMSE",
       x = NULL,
       y = "Residual TCG",
       color = "Forecast Model") +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        axis.text.x = element_blank())

mae_plot <- ggplot(mae_plot_data, aes(x = year, y = diag_mean, color = as.factor(model_num), group = as.factor(model_num))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:6) +
  labs(title = "MAE",
       x = NULL,
       y = "Residual TCG",
       color = "Forecast Model") +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        axis.text.x = element_blank())

bias_plot <- ggplot(bias_plot_data, aes(x = year, y = diag_mean, color = as.factor(model_num), group = as.factor(model_num))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:6) +
  labs(title = "Bias",
       x = "Lead Time (Forecast)",
       y = "Residual TCG",
       color = "Forecast Model") +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))

# metric vs year logs:
crps_yr_lag_plot <- ggplot(crps_plot_data, aes(x = cast_year, y = yr_one_lag, color = as.factor(model_num), group = as.factor(model_num))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_discrete(breaks = 2018:2023) +
  labs(title = "Year-Lag CRPS",
       x = NULL,
       y = NULL,
       color = "Forecast Model") +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        axis.text.x = element_blank())

rmse_yr_lag_plot <- ggplot(rmse_plot_data, aes(x = cast_year, y = yr_one_lag, color = as.factor(model_num), group = as.factor(model_num))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_discrete(breaks = 2018:2023) +
  labs(title = "Year-Lag RMSE",
       x = NULL,
       y = NULL,
       color = "Forecast Model") +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        axis.text.x = element_blank())

mae_yr_lag_plot <- ggplot(mae_plot_data, aes(x = cast_year, y = yr_one_lag, color = as.factor(model_num), group = as.factor(model_num))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_discrete(breaks = 2018:2023) +
  labs(title = "Year-Lag MAE",
       x = NULL,
       y = NULL,
       color = "Forecast Model") +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        axis.text.x = element_blank())

bias_yr_lag_plot <- ggplot(bias_plot_data, aes(x = cast_year, y = yr_one_lag, color = as.factor(model_num), group = as.factor(model_num))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_x_discrete(breaks = 2018:2023) +
  labs(title = "Year-Lag Bias",
       x = "Year",
       y = NULL,
       color = "Forecast Model") +
  theme_bw() +
  theme(plot.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12))

# plot as a group:
combine_plots <- ((crps_plot / rmse_plot / mae_plot / bias_plot) | (crps_yr_lag_plot / rmse_yr_lag_plot / mae_yr_lag_plot / bias_yr_lag_plot)) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "Forecast Performance Metrics",
                  theme = theme(plot.title = element_text(size = 16))) &
  scale_color_manual(name = "Forecast Model:",
                     values = c("1" = "#5778a4", 
                                "2" = "#e49444",
                                "3" = "#d1615d",
                                "base" = "#85b6b2"),
                     labels = c("1" = "Model 1",
                                "2" = "Model 2",
                                "3" = "Model 3",
                                "base" = "Base")) &
  theme(legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
combine_plots

# save them
# set up directory path:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/"
# save combine:
png(filename = paste0(save_dir, "2026_02_11_crps_rmse_mae_bias_plots.png"),
    height = 10,
    width = 9,
    units = "in",
    res = 600)
combine_plots
#dev.off()
# ggsave(combine_plots,
#        filename = paste0(save_dir, "2026_02_11_crps_rmse_mae_bias_plots.png"),
#        height = 15,
#        width = 14,
#        dpi = 600)

## Making table of crps/rmse/mae/bias results
all_forecast_perform <- rbind(crps_plot_data, rmse_plot_data, mae_plot_data, bias_plot_data) |>
  select(model_num, year, diag_mean, cast_year, yr_one_lag, Test) |>
  rename(Model = model_num) |>
  rename(Score = diag_mean) |>
  rename(`Year-Lag Score` = yr_one_lag) |>
  rename(`Lead Time Year` = year) |>
  rename(`Forecast Start Year` = cast_year) |>
  # rounding:
  mutate(across(c(Score, `Year-Lag Score`), function(x) round(x, 3)))
afp_table <- gt(all_forecast_perform, groupname_col = "Test") |>
  # adding breaks for column labels:
  cols_label(`Lead Time Year` = html("Lead Time<br>Year"),
             `Year-Lag Score` = html("Year-Lag<br>Score"),
             `Forecast Start Year` = html("Forecast<br>Start Year")) |>
  # make columns bold:
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) |>
  # add row spanners for tests:
  tab_style(style = cell_text(style = "italic"),
            locations = cells_row_groups()) |>
  tab_style(style = cell_fill(color = "grey86"),
            locations = cells_row_groups()) |>
  # add lines to separate two analyses:
  tab_style(style = cell_borders(sides = "right",
                                 color = "grey86",
                                 weight = px(3),
                                 style = "solid"),
            locations = list(cells_body(columns = c(Model, Score)),
                             cells_column_labels(c(Model, Score)))) |>
  # add thinner lines between other columns:
  tab_style(style = cell_borders(sides = "right",
                                 color = "grey86",
                                 weight = px(1),
                                 style = "solid"),
            locations = list(cells_body(columns = c(`Lead Time Year`, `Forecast Start Year`)),
                             cells_column_labels(c(`Lead Time Year`, `Forecast Start Year`)))) |>
  tab_options(table.font.size = 16) 
afp_table

# set up directory path:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/"
# save table:
gtsave(afp_table,
       file = paste0(save_dir, "2026_02_11_forecast_crps_all_table.rtf"))
# as csv
write.csv(all_forecast_perform, "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/Residuals/crps_rmse_mae_bias_table_all.csv")

### ARCHIVE ###

## Plotting CRPS
# x <- colnames(crps_all)[-1]
# y <- log10(crps_all[1, 2:ncol(crps_all)])
# plot(x, y, type = "l")
# for (i in 2:nrow(crps_all)){
#   lines(x, log10(crps_all[i, 2:ncol(crps_all)]))
# }

## Plotting RSME
# x <- colnames(rmse_all)[-c(1,2)]
# y <- log10(rmse_all[1, 3:ncol(rmse_all)])
# plot(x, y, type = "l")
# for (i in 2:nrow(rmse_all)){
#   lines(x, log10(rmse_all[i, 3:ncol(rmse_all)]))
# }

## Tableau 10 colors:
# blue: #5778a4,
#   orange: #e49444,
#   red: #d1615d,
#   teal: #85b6b2,
