### This is the script for running mortality analyses
### It is for the first set of analyses, using harvard forest field plot data
### Assessing drivers of plot-level mortality, using remote sensing and field observations

#### ----- Libraries ----- ####
#install.packages("librarian")
#install.packages("censReg")
librarian::shelf(tidyverse, dplyr, ggplot2, mgcv, censReg)

#### ----- Load Data (if not in environment) ----- ####
# disturbance magnitude data
# harvard forest plot data

#### ----- prepping data to run models ----- ####
### preparing plot-level information from tree data:
tree_to_plot <- trees %>%
  select(hotspot, point, plot, tot_tree, tot_dead, pdead, tot_dba, pdba) %>%
  group_by(plot) %>% summarise(across(everything(), first))

notree <- setdiff(plots$plot, tree_to_plot$plot)
#plot_out <- which(plots$plot == notree)

### making the matrix for predictors variables:
## disturbance magnitude and tcg/scores during disturbance:
cols <- c(grep("^20", names(hf_dmag)), 
          grep("^dmag", names(hf_dmag)))
dmags <- cbind(hf_dmag[cols], hf_dmag_cs[cols])[-which(plots$plot == notree),]  # remote no-tree plot
colnames(dmags) <- c("tcg_y1", "tcg_y2", "dmag_y1", "dmag_y2",
                     "cs_y1", "cs_y2", "dmag_cs_y1", "dmag_cs_y2")

## oak dead basal area:
oak_to_plot <- oaktrees %>%
  select(hotspot, point, plot, tot_oak, tba_oak, pba) %>%
  group_by(plot) %>% summarise(across(everything(), first))
oak_to_plot <- oak_to_plot[-which(plots$plot == notree),]

oak_dat <- tree_to_plot %>%
  select(plot, hotspot, point) %>%
  left_join(., oak_to_plot, by = 'plot', na_matches = "na") %>%
  replace(is.na(.), 0)

## daymet data:
# making seasonal values:
spring = c(3,4,5)
summer = c(6,7,8)
# spring seasonal aggregates:
spr_daym <- hf_daym %>%
  rename(spr_precip = prcp..mm.day., spr_tmax = tmax..deg.c., spr_tmin = tmin..deg.c., spr_vpd = vp..Pa.) %>%
  group_by(site, year) %>% filter(month %in% spring) %>%
  summarise(across(c(spr_precip, spr_tmax, spr_tmin, spr_vpd), mean))
# summer seasonal aggregates:
summ_daym <- hf_daym %>%
  rename(summ_precip = prcp..mm.day., summ_tmax = tmax..deg.c., summ_tmin = tmin..deg.c., summ_vpd = vp..Pa.) %>%
  group_by(site, year) %>% filter(month %in% summer) %>%
  summarise(across(c(summ_precip, summ_tmax, summ_tmin, summ_vpd), mean))
# spring/summer monthly means for 2016:
daym16 <- hf_daym %>%
  filter(., month %in% c(spring, summer) & year == 2016) %>%
  rename(precip = prcp..mm.day., tmax = tmax..deg.c., tmin = tmin..deg.c., vpd = vp..Pa.) %>%
  mutate(date = month.abb[month]) %>% 
  pivot_wider(names_from = date, 
              values_from = c(precip, tmax, tmin, vpd)) %>%
  group_by(site) %>% summarise(across(everything(), sum, na.rm = T)) %>%
  ungroup()%>%
  rename_at(vars(-c(site, year, month)), ~paste0("2016_", .))
  
# combine data for seasonal 2014-2015:  
daym <- cbind(spr_daym[which(spr_daym$year == 2014), grep("^spr", names(spr_daym))], 
              summ_daym[which(summ_daym$year == 2014), grep("^summ", names(summ_daym))]) %>%
                rename_with(~paste0("2014_", .x)) %>%
  bind_cols(spr_daym[which(spr_daym$year == 2015), grep("^spr", names(spr_daym))],
            summ_daym[which(summ_daym$year == 2015), grep("^summ", names(summ_daym))]) %>%
  rename_at(vars(starts_with('s')), ~paste0('2015_', .)) %>%
  bind_cols(spr_daym[which(spr_daym$year == 2016), grep("^spr", names(spr_daym))],
            summ_daym[which(summ_daym$year == 2016), grep("summ", names(summ_daym))]) %>%
  rename_at(vars(starts_with('s')), ~paste0('2016_', .))

# combine daymet data:
dayms <- cbind(daym, daym16[,grep("^2", names(daym16))])[-which(plots$plot == notree),]
# clean up:
rm(daym, daym16, spr_daym, summ_daym)

## make predictor variable data set:
pred <- cbind(dmags,                           # disturbance magnitude and disturbance year tcg/scores data
              oak_dat$tba_oak, oak_dat$pba,    # oak tree density data
              dayms)                           # daymet data

## make response variable data set:
resp <- cbind.data.frame(plot = tree_to_plot$plot, 
                         hotspot = tree_to_plot$hotspot,
                         dead = tree_to_plot$tot_dead > 0,
                         pdead = tree_to_plot$pdead, 
                         pdba = tree_to_plot$pdba) %>%
  mutate(dead_bin = ifelse(dead == TRUE, 1, 0), .after = 2)

#### ----- running models ----- ####
# using beta regression in mgcv package
# function example: gam(y ~ s(x), family = betar(link = "logit"), data = data)

test_data <- cbind.data.frame(y = resp$pdead, x = pred$tcg_y1)
test <- gam(y ~ s(x), family = betar(link = "logit"), data = test_data)

