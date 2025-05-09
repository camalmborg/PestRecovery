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
  geom_point() +
  geom_point(aes(x = x, y = y_pred),
             col = "red") +
  geom_line(aes(x = x, y = ci[1,])) +
  geom_line(aes(x = x, y = ci[3,]))
pdba_pred_obs

pred_obs_compare <- ggplot(plot_data, aes(x = y_pred, y = y_obs)) +
  geom_point()
pred_obs_compare
