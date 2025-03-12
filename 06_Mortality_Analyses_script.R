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
## preparing plot-level information from tree data:
tree_to_plot <- trees %>%
  select(hotspot, point, plot, tot_tree, tot_dead, pdead, tot_dba, pdba) %>%
  group_by(plot) %>% summarise(across(everything(), first))

notree <- setdiff(plots$plot, tree_to_plot$plot)

## making the matrix for predictors variables:
# disturbance magnitude and tcg/scores during disturbance:
cols <- c(grep("^20", names(hf_dmag)), 
          grep("^dmag", names(hf_dmag)))
dmags <- cbind(hf_dmag[cols], hf_dmag_cs[cols])[-which(plots$plot == notree),]  # remote no-tree plot
colnames(dmags) <- c("tcg_y1", "tcg_y2", "dmag_y1", "dmag_y2",
                     "cs_y1", "cs_y2", "dmag_cs_y1", "dmag_cs_y2")

## daymet data:
# making seasonal values:
spring = c(3,4,5)
summer = c(6,7,8)
seas_daym <- hf_daym %>%
  group_by(site, year) %>% filter(month %in% spring) #%>%
  group_by(year) %>% aggregate()

# combine data:

#### ----- running models ----- ####
# using beta regression in mgcv package
# function example: gam(y ~ s(x), family = betar(link = "logit), data = data)

