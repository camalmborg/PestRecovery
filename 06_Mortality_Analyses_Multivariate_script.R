### This is the script for running multivariate mortality analyses
### It is for the first set of analyses, using harvard forest field plot data
### Assessing drivers of plot-level mortality, using remote sensing and field observations

#### ----- Libraries ----- ####
#install.packages("librarian")
#install.packages("rjags")
librarian::shelf(tidyverse, dplyr, rjags, ggplot2) 