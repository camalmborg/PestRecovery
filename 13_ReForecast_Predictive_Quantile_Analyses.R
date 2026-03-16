### Predictive Quantiles analysis ###

# load libraries
librarian::shelf(dplyr, tidyverse, stringr, ggplot2)

## Getting TCG observations
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


## Comparing forecasts with observations:
# set wd:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

## Load forecast
forecasts <- grep("^2025", list.files("Recovery_Forecasts/"))
# setting task id for running array job on SCC:
task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
# loading forecast:
file <- list.files("Recovery_Forecasts/")[forecasts[task_id]]
# year start:
year <- as.numeric(str_extract(file, "(?<=year_)\\d+(?=_)"))
fyears <- year:2023
# load forecast:
f <- read.csv(paste0("Recovery_Forecasts/", file))
# add back baseline to get TCG predictions:
fcast <- f |>
  # add baseline to predictions:
  mutate(across(-site, ~ tcg_base$baseline[site] - .x))
# rename columns:
colnames(fcast) <- c("site", as.character(fyears))

# match tcg years:
tcg_f <- tcg[,c(names(fcast)[-1])]

## Get predictive quantiles
# collecting quantiles:
qsite <- matrix(data = NA, nrow = nrow(tcg_f), ncol = ncol(tcg_f))
colnames(qsite) <- (as.character(fyears))
## Compare each ensemble member for each site
for (i in 1:nrow(tcg_f)){
  for (j in 1:ncol(tcg_f)){
    fsite <- fcast[which(fcast$site == i), -1][,j]
    qsite[i,j] <-  sum(fsite < tcg_f[i,j])/1500 
  }
}


## Saving result
# set save directory:
save_dir <- "Recovery_Forecasts/Predictive_Quantiles/"
# save:
file_name <- gsub(".csv", "", file)
write.csv(qsite, file = paste0(save_dir, file_name, "_pred_quant.csv"))



### For Making Summary Histograms ###

# load some predictive quantiles
pdqs <- list.files("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/Recovery_Forecasts/Predictive_Quantiles/")
pdq <- read.csv(paste0("Recovery_Forecasts/Predictive_Quantiles/",pdqs[2]))[-1]

# make histogram data:
pdq_hist_data <- pdq |>
  # # remote first column:
  select(-1) |>
  # fix column names:
  rename_with(~ str_replace_all(., c("X" = "", "\\." = "-"))) |>
  # add site column:
  mutate(site = 1:nrow(pdq), .before = 1) |>
  # pivot:
  pivot_longer(-site, names_to = "year", values_to = "quants")
  
# make histogram:
pdq_hist <- ggplot(data = pdq_hist_data, aes(x = quants)) +
  geom_histogram(bins = 25, fill = "slategray2") +
  facet_wrap(~year, nrow = 1, ncol = 2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  labs(x = "Quantiles", y = "Count") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14))
pdq_hist

# save it:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/Predictive_Quantiles/"
setwd(save_dir)


# # save:
# png(filename = paste0(save_dir, "2026_03_15_predictive_histogram_forecast_vs_tcg.png"),
#     height = 6,
#     width = 8,
#     units = "in",
#     res = 600)
# 
# tf_df |>
#   pivot_longer(everything(), names_to = "column", values_to = "value") |>
#   count(column, value) |>
#   filter(!is.na(value)) |>
#   ggplot(aes(x = column, y = n, fill = factor(value, levels = c(0,1),
#                                               labels = c("Underestimated","Overestimated")))) +
#   scale_fill_manual(values = c("Underestimated" = "skyblue1", "Overestimated" = "indianred1")) +
#   geom_col(position = position_dodge(width = 0.6), width = 0.55) +
#   labs(x = "Year", y = "Count", fill = "Value") +
#   theme_bw() +
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         legend.title = element_blank())
# dev.off()