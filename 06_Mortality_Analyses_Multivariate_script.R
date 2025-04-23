### This is the script for running multivariate mortality analyses
### It is for the first set of analyses, using harvard forest field plot data
### Assessing drivers of plot-level mortality, using remote sensing and field observations

#### ----- Libraries ----- ####
#install.packages("librarian")
#install.packages("rjags")
librarian::shelf(tidyverse, dplyr, rjags, ggplot2, combinat) 

#### ----- Load Data (if not in environment) ----- ####
# set working directory:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)
load("Environments/2025_04_23_environment.RData")

#### ----- Preparing data for models ----- ####
## make predictor data set:
pred_m <- pred %>%
  # select the top performers:
  select(pba_oak, tba_oak, dmag_y1, dmag_y2, dmag_cs_y1, dmag_cs_y2, cs_y1) %>%
  # add the sums of dmag columns (tcg dmag, cs dmag)
  mutate(dmag_y1_y2 = dmag_y1 + dmag_y2, .after = dmag_y2) %>%
  mutate(dmag_cs_y1_y2 = dmag_cs_y1 + dmag_cs_y2, .after = dmag_cs_y2)

## make matrix of models using column numbers:
# making data frame for up to 3-variable models to start
# two variable models:
two_var_models <- as.data.frame(t(cbind(combn(c(1, 3:9), 2),
                                        combn(c(2, 3:9), 2),
                                        combn(c(3, 4:9), 2)))) %>%  # combinations of variables
  # filter all with 1, 2, 3 as 'starting' variables
  filter(V1 %in% 1:3) %>%
  # select only distinct rows
  distinct(V1, V2, .keep_all = TRUE) %>%
  # add a placeholder for 3rd variable to join together
  mutate(V3 = NA)
  
# three variable models:
three_var_models <- as.data.frame(t(cbind(combn(c(1, 3:9), 3),
                                          combn(c(2, 3:9), 3)))) %>% # combinations of variables
  # filter all with 1, 2 starting and variable 2 = 3
  filter(V1 %in% 1:2 & V2 == 3) %>%
  # select only distinct rows
  distinct(V1, V2, V3, .keep_all = TRUE)

# make one data frame with all models listed:
mult_models <- rbind(two_var_models, three_var_models)
