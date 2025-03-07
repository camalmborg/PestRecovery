### This is the script for running mortality analyses
### It is for the first set of analyses, using harvard forest field plot data
### Assessing drivers of plot-level mortality, using remote sensing and field observations

#### ----- Libraries ----- ####
#install.packages("librarian")
librarian::shelf(tidyverse, dplyr, ggplot2, mgcv)

#### ----- Load Data (if not in environment) ----- ####
# disturbance magnitude data
# harvard forest plot data

#### ----- prepping data to run models ----- ####
# preparing plot-level information from tree data:
tree_to_plot <- trees %>%
  select(hotspot, point, plot, tot_tree, tot_dead, pdead, tot_dba, pdba) %>%
  group_by(plot) %>% summarise(across(everything(), first, .groups = "keep"))

notree <- setdiff(plots$plot, tree_to_plot$plot)

#### ----- running models ----- ####
# using beta regression in mgcv package
# function example: gam(y ~ s(x), family = betar(link = "logit), data = data)

