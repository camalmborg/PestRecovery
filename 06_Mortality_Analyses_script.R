### This is the script for running mortality analyses
### It is for the first set of analyses, using harvard forest field plot data
### Assessing drivers of plot-level mortality, using remote sensing and field observations

#### ----- Libraries ----- ####
#install.packages("librarian")
librarian::shelf(tidyverse, dplyr, googledrive, ggplot2)

#### ----- Load Data (if not in environment) ----- ####
