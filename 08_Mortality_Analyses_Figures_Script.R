# Script for making figures for mortality analyses

# load libraries
#install.packages("ggplot2")
#install.packages("hrbrthemes")
#install.packages("ggridges")
librarian::shelf(tidyverse, dplyr, ggplot2, RColorBrewer, hrbrthemes, ggridges)

# # navigate to folder:
# dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
# setwd(dir)
#load("Environments/2025_05_08_environment.RData")
#load("Environments/2025_04_07_environment.RData")
load("Envrionments/2025_05_15_environment.RData")

## plotting percent dead basal area observed and predicted
# make the data set:
plot_data <- data.frame(
  x = 1:nrow(data_sort),
  y_obs = data_sort$pdba,
  y_pred = apply(ypred, 2, mean),
  hot = data_sort$hotspot,
  lat = data_sort$lat,
  lon = data_sort$lon,
  plot = data_sort$plot
)
# making plot:
pdba_pred_obs <- ggplot(plot_data, aes(x = x, y = y_obs)) +
  geom_ribbon(aes(ymin = ci[1,], ymax = ci[3,]), 
              fill = "dodgerblue2", alpha = 0.25) +
  geom_point() +
  geom_point(aes(x = x, y = y_pred),
             col = "dodgerblue2") +
  ggtitle("Mortality Observed (black) and Predicted (blue)") +
  xlab("Plots (sorted)") +
  ylab("Percent dead basal area in plot") +
  theme_classic()
pdba_pred_obs

# predicted vs observed plot
pred_obs_compare <- ggplot(plot_data, aes(x = y_pred, y = y_obs)) +
  geom_point() +
  xlab("Predicted percent dead basal area in plot") +
  ylab("Observed percent dead basal area in plot") +
  xlim(0,0.87) +
  ylim(0,0.87) +
  geom_abline (intercept = 0, slope=1, linetype = "dashed", color = "firebrick1") +
  theme_bw()
pred_obs_compare

### Histograms and density plots
# mortality across all plots
hist_mort_obs <- ggplot(data_sort, aes(x = pdba)) +
  geom_histogram(binwidth = 0.01, bins = 50,
                 color = "dodgerblue4", fill = "dodgerblue2") +
  xlim(0, 0.9) +
  ylim(0, 6.5) +
  xlab("Percent dead basal area across plots") +
  theme_bw()
hist_mort_obs

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
  theme_bw()
ridge_across_hotspots

# mortality by hotspot density plot:
density_across_hotspots <- ggplot(data_sort, aes(x = pdba,fill = as.factor(hotspot))) +
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
       filename= paste0(save_dir, "/2025_05_15_pred_obs_ci_plot.png"),
       width=8,
       height=6)

# predicted vs observed plot:
ggsave(plot = pred_obs_compare,
       filename= paste0(save_dir, "/2025_05_15_pred_obs_compare_plot.png"),
       width=6,
       height=6)

# histogram of percent dead basal area across plots:
ggsave(plot = hist_mort_obs,
       filename = paste0(save_dir, "/2025_05_15_hist_pdba_across_plots.png"),
       width = 8,
       height = 6)

# ridgeline plot of mortality across hotspots:
ggsave(plot = ridge_across_hotspots,
       filename = paste0(save_dir, "/2025_05_15_hist_pdba_ridges_hotspots.png"),
       width = 8,
       height = 6)

# densities of mortality across hotspots:
ggsave(plot = density_across_hotspots,
       filename = paste0(save_dir, "/2025_05_15_hist_pdba_density_hotspots.png"),
       width = 8,
       height = 6)


### Results Table and Heatmap