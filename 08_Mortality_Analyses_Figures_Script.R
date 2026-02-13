# Script for making figures for mortality analyses

# load libraries
#install.packages("ggplot2")
#install.packages("hrbrthemes")
#install.packages("ggridges")
#install.packages("gt")
#install.packages("gtExtras)
librarian::shelf(tidyverse, dplyr, ggplot2, RColorBrewer, hrbrthemes, 
                 ggridges, gt, gtExtras, webshot2)

# # navigate to folder:
# dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
# setwd(dir)
#load("Environments/2025_05_08_environment.RData")
#load("Environments/2025_04_07_environment.RData")
load("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Environments/2025_05_15_environment.RData")

## plotting percent dead basal area observed and predicted
# make the data set:
plot_data <- data.frame(
  x = 1:nrow(data_sort),
  y_obs = (data_sort$pdba * 100),
  y_pred = (apply(ypred, 2, mean) * 100),
  hot = data_sort$hotspot,
  lat = data_sort$lat,
  lon = data_sort$lon,
  plot = data_sort$plot
)

# making plot:
pdba_pred_obs <- ggplot(plot_data, aes(x = x, y = y_obs)) +
  geom_ribbon(aes(ymin = (ci[1,] * 100), ymax = (ci[3,] * 100)), 
              fill = "dodgerblue2", alpha = 0.25) +
  geom_point() +
  geom_point(aes(x = x, y = y_pred),
             col = "dodgerblue2") +
  ggtitle("Mortality Observed (black) and Predicted (blue)") +
  xlab("Plots (sorted)") +
  ylab("Percent Dead Basal Area in Plot") +
  theme_bw() +
  theme(plot.title = element_text(size = 15),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
pdba_pred_obs

# predicted vs observed plot
pred_obs_compare <- ggplot(plot_data, aes(x = y_pred, y = y_obs)) +
  geom_point(size = 1.5) +
  ggtitle("Mortality Predicted vs. Observed") +
  xlab("Predicted Percent Dead Basal Area in Plot") +
  ylab("Observed Percent Dead Basal Area in Plot") +
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline (intercept = 0, slope = 1, linewidth = 1, linetype = "dashed", color = "firebrick1") +
  theme_bw() +
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
pred_obs_compare

### Histograms and density plots
# mortality across all plots
hist_mort_obs <- ggplot(data_sort, aes(x = (pdba * 100))) +
  geom_histogram(binwidth = 1, bins = 25,
                 fill = "dodgerblue2") +
  xlim(0, 100) +
  ylim(0, 6.5) +
  ggtitle("Percent Dead Basal Area Density") +
  xlab("Percent Dead Basal Area across Plots") +
  ylab("Count") +
  theme_bw() +
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
hist_mort_obs

density <- ggplot(data_sort, aes(x = (pdba * 100), y = (after_stat(density)))) +
  geom_density(color = "dodgerblue", fill = "dodgerblue", alpha = 0.35) +
  ggtitle("Percent Dead Basal Area Density across Plots") +
  xlab("Percent Dead Basal Area") +
  ylab("Density") +
  xlim(0, 100) +
  theme_bw() +
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
density

# mortality by hotspot ridge plot:
ridge_across_hotspots <- ggplot(data_sort, aes(x = pdba, y = fct_rev(as.factor(hotspot)), 
                                              fill = as.factor(hotspot))) +
  geom_density_ridges(scale = 2) +
  scale_fill_manual(values=c("khaki1", "deeppink", "forestgreen", "darkorchid4", "sienna1", "cornflowerblue"),
                    name = "Hotspot") +
  scale_y_discrete(expand = expand_scale(add = c(.5, 2.5))) +
  xlab("Percent dead basal area") +
  ylab("Hotspot") +
  ggtitle("Dead basal area by hotspot") +
  xlim(0, 1) +
  theme_bw() +
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
ridge_across_hotspots

# mortality by hotspot density plot:
density_across_hotspots <- ggplot(data_sort, aes(x = pdba, fill = as.factor(hotspot))) +
  geom_density(alpha = 0.35) +
  scale_fill_manual(values=c("sienna1", "deeppink", "forestgreen", "darkorchid4", "khaki1", "cornflowerblue"),
                    name = "Hotspot") +
  xlab("Percent dead basal area") +
  ggtitle("Dead basal area by hotspot") +
  xlim(0, 1) +
  theme_bw()
density_across_hotspots

### Saving the plots
# set up directory path:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures"

# predicted and observed with confidence interval:
png(filename = paste0(save_dir, "/2026_02_09_pred_obs_ci_plot_square.png"),
    width = 6,
    height = 6, 
    units = "in",
    res = 600)
pdba_pred_obs
dev.off()
# ggsave(plot = pdba_pred_obs,
#        filename= paste0(save_dir, "/2026_02_09_pred_obs_ci_plot_square.png"),
#        width=6,
#        height=6)

# predicted vs observed plot:
png(filename = paste0(save_dir, "/2026_02_09_pred_obs_compare_plot.png"),
    width = 6,
    height = 6, 
    units = "in",
    res = 600)
pred_obs_compare
dev.off()
# ggsave(plot = pred_obs_compare,
#        filename= paste0(save_dir, "/2026_02_09_pred_obs_compare_plot.png"),
#        width=6,
#        height=6)

# histogram of percent dead basal area across plots:
png(filename = paste0(save_dir, "/2026_02_09_hist_pdba_across_plots.png"),
    width = 6,
    height = 6, 
    units = "in",
    res = 600)
hist_mort_obs
dev.off()
# ggsave(plot = hist_mort_obs,
#        filename = paste0(save_dir, "/2026_02_09_hist_pdba_across_plots.png"),
#        width = 8,
#        height = 6)

# ridgeline plot of mortality across hotspots:
png(filename = paste0(save_dir, "/2026_02_09_hist_pdba_ridges_hotspots.png"),
    width = 6,
    height = 8, 
    units = "in",
    res = 600)
ridge_across_hotspots
dev.off()
# ggsave(plot = ridge_across_hotspots,
#        filename = paste0(save_dir, "/2026_02_09_hist_pdba_ridges_hotspots.png"),
#        width = 8,
#        height = 6)

# densities of mortality across hotspots:
png(filename = paste0(save_dir, "/2026_02_11_hist_pdba_density.png"),
    width = 6,
    height = 6, 
    units = "in",
    res = 600)
density
dev.off()
# ggsave(plot = density,
#        filename = paste0(save_dir, "/2026_01_09_hist_pdba_density.png"),
#        width = 8,
#        height = 6)


### Results Table and Heatmap
# make a data frame with top model information:
top_models <-data.frame(Rank = 1:5, 
                        Variables = c("PBA oak, DMag Y1", "PBA oak, DMag CS Y1 + Y2", "PBA oak, DMag Y1 + Y2", "PBA oak, DMag CS Y1", "PBA oak"),
                        "ΔDIC" = round(all_dics$delDIC[1:5], 3))
best_models <- multi_results$results[c(11, 1, 10, 12),]
best_uni <- uni_results$results[3,]
best_results <- rbind(best_models, best_uni)
#betas <- grep("^beta", colnames(multi_results$results))
top_models$int <-round(as.numeric(best_results[,"beta_int"]), 3)
top_models$b1 <- round(as.numeric(best_results[, "beta_1"]), 3)
top_models$b2 <- round(as.numeric(best_results[, "beta_2"]), 3)

# make the table:
results_table <- gt(top_models) %>%
  tab_spanner(
    label = "Model Performance",
    columns = c(Rank, Variables, "ΔDIC")) %>%
  tab_spanner(
    label = "Parameter Estimates",
    columns = c(int, b1, b2)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners()) %>%
  tab_options(table.font.size = 16) #%>%
  #gt_theme_538()
results_table
# save table
gtsave(results_table,
       file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/2026_02_09_mortality_results_table.html")


### Making univariate and multivariate mortality results tables:
# univariate names:
mort_models <- data.frame(covariates = c(1:13),
                          names = c("tcg_2016", "tcg_2017", "dmag_tcg_2016", "dmag_tcg_2017",
                                    "cs_2016", "cs_2017", "dmag_cs_2016", "dmag_cs_2017",
                                    "tba_oak", "pba_oak", "base", "dmag_tcg_2016_2017", "dmag_cs_2016_2017"))
unilookup <- setNames(mort_models$names, mort_models$covariates)

# univariate table:
uni_mort_results_no_log <- as.data.frame(uni_results$results[grep("nolog", uni_results$results[,2]),]) |>
  # match names:
  mutate(Model = gsub("_nolog", "", Model)) |>
  mutate(Model = as.numeric(Model)) |>
  arrange(Model) |>
  mutate(covariate = mort_models$names[-c(11:13)], .after = Model) |>
  # add log/no log:
  mutate(Logit = "no logit", .after = covariate) |>
  # select columns:
  select(-c(beta_2, beta_3)) |>
  # rename columns:
  rename(Intercept = beta_int) |>
  rename(Beta = beta_1)
uni_mort_results_log <- as.data.frame(uni_results$results[-grep("nolog", uni_results$results[,2]),]) |>
  # match names:
  mutate(Model = gsub("_log", "", Model)) |>
  mutate(Model = as.numeric(Model)) |>
  arrange(Model) |>
  mutate(covariate = mort_models$names, .after = Model) |>
  # add log/no log:
  mutate(Logit = "logit", .after = covariate) |>
  # select columns:
  select(-c(beta_2, beta_3)) |>
  # rename columns:
  rename(Intercept = beta_int) |>
  rename(Beta = beta_1)
# combine:
all_uni_mort_models <- rbind(uni_mort_results_no_log, uni_mort_results_log) |>
  mutate(Model = 1:nrow(all_uni_mort_models)) |>
  # select only what we want:
  select(-Run) |>
  # convert columns to numeric
  mutate(across(-c(Model, covariate, Logit), as.numeric)) |>
  mutate(across(where(is.numeric), ~ round(.x, digits = 4)))

uni_mort_final_table <- gt(all_uni_mort_models)
gtsave(uni_mort_final_table,
       file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/2026_02_12_uni_mortality_params_table.rtf")
write.csv(all_uni_mort_models, file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Mortality_Model_Runs/all_uni_mort_models.csv")
  
# multivariate table:
# multivariate names:
mort_multi_models <- data.frame(covariates = c(1:6),
                                names = c("pba_oak", "tba_oak", "dmag_cs_2016_2017", "dmag_tcg_2016_2017", "dmag_tcg_2016", "dmag_cs_2016"))
# make a look up table:
lookup <- setNames(mort_multi_models$names, mort_multi_models$covariates)

# multivariate table:
multi_mort_results <- as.data.frame(multi_results$results) |>
  # matching to covariates:
  mutate(Model = gsub(".*covs_", "", Model)) |>
  separate(Model, into = c("covariate 1", "covariate 2", "covariate 3"), sep = "_") |>
  # fill in with lookup table:
  mutate(across(c(`covariate 1`, `covariate 2`, `covariate 3`), ~ lookup[.x])) |>
  # remove the ones with overlapping covariates:
  slice(-c(4:9, 13:15)) |>
  # select columns:
  select(-c(Run, `covariate 3`, beta_3)) |>
  # rename columns:
  rename(Intercept = beta_int) |>
  mutate(Model = 1:nrow(multi_mort_results), .before = 1) |>
  # convert columns to numeric
  mutate(across(-c(`covariate 1`, `covariate 2`), as.numeric)) |>
  mutate(across(where(is.numeric), ~ round(.x, digits = 4)))

multi_mort_final_table <- gt(multi_mort_results)
gtsave(multi_mort_final_table,
       file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/2026_02_12_multi_mortality_params_table.rtf")
write.csv(multi_mort_results, file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Mortality_Model_Runs/all_multi_mort_models.csv")


## Fixing mortality dic table:
all_dics_uni_log <- all_dics |>
  filter(type == "uni", mod_info == "log") |>
  # fill in with lookup table:
  mutate(across(c(model), ~ unilookup[.x])) |>
  rename(`covariate 1` = model) |>
  # make logit column
  mutate(logit = "logit", .after = `covariate 1`) |>
  mutate(`covariate 2` = NA, .after = `covariate 1`) |>
  # select columns we want:
  select(-c(type, mod_info)) |>
  # round DICs and delDIC:
  mutate(across(where(is.numeric), ~ round(.x, digits = 4)))

all_dics_uni_nolog <- all_dics |>
  filter(type == "uni", mod_info == "nolog") |>
  # fill in with lookup table:
  mutate(across(c(model), ~ unilookup[.x])) |>
  rename(`covariate 1` = model) |>
  mutate(logit = "no logit", .after = `covariate 1`) |>
  mutate(`covariate 2` = NA, .after = `covariate 1`) |>
  # select columns we want:
  select(-c(type, mod_info)) |>
  # round DICs and delDIC:
  mutate(across(where(is.numeric), ~ round(.x, digits = 4)))

all_dics_multi <- all_dics |>
  filter(type != "uni") |>
  separate(mod_info, into = c("covariate 1", "covariate 2", "covariate 3"), sep = " ") |>
  # fill in with lookup table:
  mutate(across(c(`covariate 1`, `covariate 2`, `covariate 3`), ~ lookup[.x])) |>
  # remove the ones with overlapping covariates:
  slice(-c(5:7, 11:12, 14, 15:17)) |>
  # select columns we want:
  select(-c(model, type, `covariate 3`)) |>
  mutate(logit = "logit", .after = `covariate 2`) |>
  # round DICs and delDIC:
  mutate(across(where(is.numeric), ~ round(.x, digits = 4)))

# rebind them and sort:
all_mort_dics_table <- rbind(all_dics_uni_log, all_dics_uni_nolog, all_dics_multi) |>
  arrange(delDIC) |>
  mutate(perform = 1:nrow(all_mort_dics_table))

final_mort_dic_table <- gt(all_mort_dics_table)
gtsave(final_mort_dic_table,
       file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/2026_02_12_all_mort_dic_table.rtf")
write.csv(all_mort_dics_table, file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Mortality_Model_Runs/all_mort_dics.csv")

