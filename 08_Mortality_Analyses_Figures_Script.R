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
ggsave(plot = pdba_pred_obs,
       filename= paste0(save_dir, "/2026_02_09_pred_obs_ci_plot_square.png"),
       width=6,
       height=6)

# predicted vs observed plot:
ggsave(plot = pred_obs_compare,
       filename= paste0(save_dir, "/2026_02_09_pred_obs_compare_plot.png"),
       width=6,
       height=6)

# histogram of percent dead basal area across plots:
ggsave(plot = hist_mort_obs,
       filename = paste0(save_dir, "/2026_02_09_hist_pdba_across_plots.png"),
       width = 8,
       height = 6)

# ridgeline plot of mortality across hotspots:
ggsave(plot = ridge_across_hotspots,
       filename = paste0(save_dir, "/2026_02_09_hist_pdba_ridges_hotspots.png"),
       width = 8,
       height = 6)

# densities of mortality across hotspots:
ggsave(plot = density,
       filename = paste0(save_dir, "/2026_01_09_hist_pdba_density.png"),
       width = 8,
       height = 6)


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

