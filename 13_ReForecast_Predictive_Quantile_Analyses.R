### Predictive Quantiles analysis

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

# load forecast:
f <- read.csv("Recovery_Forecasts/2025-12-03_ens_1500_model_1_start_year_2017_reforecast_result.csv")
fcast <- f %>%
# add baseline to predictions:
mutate(across(-site, ~ tcg_base$baseline[site] - .x))




# where forecast > tcg values:
tf_df <- as.data.frame((y_pred[, -1] > tcg) * 1)
colnames(tf_df) <- c(as.character(2017:2023))
tf_sums <- colSums(tf_df, na.rm = T)

# save it:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/"
setwd(save_dir)
# save:
png(filename = paste0(save_dir, "2026_03_15_predictive_histogram_forecast_vs_tcg.png"),
    height = 6,
    width = 8,
    units = "in",
    res = 600)

tf_df |>
  pivot_longer(everything(), names_to = "column", values_to = "value") |>
  count(column, value) |>
  filter(!is.na(value)) |>
  ggplot(aes(x = column, y = n, fill = factor(value, levels = c(0,1),
                                              labels = c("Underestimated","Overestimated")))) +
  scale_fill_manual(values = c("Underestimated" = "skyblue1", "Overestimated" = "indianred1")) +
  geom_col(position = position_dodge(width = 0.6), width = 0.55) +
  labs(x = "Year", y = "Count", fill = "Value") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),      
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 14),
        legend.title = element_blank())
dev.off()