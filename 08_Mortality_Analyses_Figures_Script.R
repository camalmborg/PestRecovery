# Script for making figures for mortality analyses

# load libraries
#install.packages("ggplot2")
#install.packages("hrbrthemes")
librarian::shelf(tidyverse, dplyr, ggplot2, RColorBrewer, hrbrthemes)

# # navigate to folder:
# dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
# setwd(dir)

## plotting percent dead basal area observed and predicted
# make the data set:
plot_data <- data.frame(
  x = 1:nrow(data_sort),
  y_obs = data_sort$pdba,
  y_pred = apply(ypred, 2, mean),
  hot = data_sort$hotspot
)
# making plot:
pdba_pred_obs <- ggplot(plot_data, aes(x = x, y = y_obs)) +
  geom_ribbon(aes(ymin = ci[1,], ymax = ci[3,]), 
              fill = "dodgerblue2", alpha = 0.50) +
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

### Histograms
