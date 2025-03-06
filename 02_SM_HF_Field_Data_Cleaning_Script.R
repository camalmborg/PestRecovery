# This is the script for cleaning Harvard Forest field data from Summer 2022
# Clean up: data organized, demographic calculations for running mortality analyses

#### ----- Load libraries: ----- ####
#install.packages("librarian")
librarian::shelf(tidyverse, googledrive)

#### ----- Load data: ----- ####
trees_raw <- read.csv("Data/hf_field_data/Tree_data.csv")
seedlings_raw <- read.csv("Data/hf_field_data/Oak_regen_data.csv")
understory_raw <- read.csv("Data/hf_field_data/Understory_data.csv")
plots_raw <- read.csv("Data/hf_field_data/2022_hf_plot_data.csv")

#### ----- Clean and prep data: ----- ####
## Plot data:
# select plots that are from 2022 data collection:
plots <- plots_raw[which(!is.na(plots_cor$latitude_2022)),]  # grab just the 2022 field sites
# clean up for analyses:
plots <- plots %>%
  mutate(plot = paste(hotspot, point, sep = "-"), .after = point) %>%  # add column for plot number
  mutate(latitude = as.numeric(str_replace(latitude_2022, "_2022", "")), .after = plot) %>%  # clean lat/long columns
  mutate(longitude = as.numeric(str_replace(longitude_2022, "_2022", "")), .after = latitude) %>%
  select(-c(latitude_2017, latitude_2022, longitude_2017, longitude_2022)) %>% # remove old lat/long columns
  mutate(invasives = ifelse(invasives_2022 == "no", 0, 1)) %>%  # change invasives from "y/n" to "1/0"
  mutate(recent_timber_harvest = ifelse(recent_harv_plot == "no", 0, 1)) %>% #change timber harvest from "y/n" to "1/0"
  mutate(harv_since_2017 = ifelse(harv_since_2017_GIS == "no", 0, 1)) %>%
  select(-c(meas_2017, date_2017, observ_2017, recent_harv_plot, harv_since_2017_GIS))

## Trees data:
# select plots from 2022:
trees <- trees_raw[which(trees_raw$year == 2022),]
BAF = 10  # basal area factor from variable radius plots
trees <- trees %>%
  mutate(plot = paste(hotspot, point, sep = "-"), .after = point) %>%  # add column for plot number
  select(-c(ref_clss_2017, conf_2017)) %>%
  mutate(ba_m2 = ((dbh_cm/100)^2) * 3.14, .after = dbh_cm) %>% # basal area calculation (m^2) from dbh in cm
  group_by(hotspot,plot) %>% mutate(treeper  = n()) %>%  # number of trees per plot
  group_by(hotspot) %>% mutate(hotcount = length(unique(plot))) %>%  # number of plots per hotspot
  group_by(plot) %>% mutate(tba_m2 = sum(ba_m2)) %>%
  mutate(exp_fac = BAF / ba_m2 / hotcount) %>% # expansion factor, BAF/BA/number of plots
  group_by(hotspot) %>% mutate(TPA = sum(exp_fac))

## Oak-specific data:
oaks <- unique(trees$genusp)[grep("^QUER", unique(trees$genusp))]  # all oak species

percent <- function(x, count){  # percentage function for calculating ba and dead tree percentages
  (x / count)*100
}

# make prep oak tree data
oaktrees_raw <- subset(trees, genusp %in% oaks)
oaktrees <- oaktrees_raw %>% 
  group_by(plot) %>% mutate(tot_oak = n(),  # number of oaks
                            dead = ifelse(cond_2022 == "L", 0, 1),
                            tba_oak = sum(ba_m2),  # total pak basal area per plot
                            pba = percent(tba_oak,tba_m2),  # percent oak basal area in plot
                            deadoak = length(which(cond_2022 == "D")),  # total dead oaks
                            pdead = percent(deadoak,treeper),  # percent dead oak of all trees in plot
                            pdo = percent(deadoak,tot_oak)) %>%  # percent dead oak of all oaks in plot
  group_by(plot,cond_2022) %>% mutate(dba = sum(ba_m2)*dead,
                                 pdba = percent(dba,tba_m2)) %>%
  group_by(hotspot) %>% mutate(tba_oak_hot = sum(ba_m2),  # total oak basal area in hotspot
                               pba_hot = percent(tba_oak,tba_m2),  # percent oak basal area in hotspot
                               deadoak_hot = sum(deadoak),  # sum dead oak in hotspot
                               pdo_hot = percent(deadoak,tot_oak),  # percent dead oak in hotspot
                               dba_hot = sum(dba),  # total dead basal area in hot spot
                               pdba_hot = percent(dba_hot,sum(tba_m2))) # percent dead basal area in hotspot


## Understory data prep:
# get 2022 sites:
understory <- understory_raw[which(understory_raw$year == 2022),] # grab just the 2022 field sites
# make plot column
understory <- understory %>%
  mutate(plot = paste(hotspot, point, sep = "-"), .after = point)  # add column for plot numbe

undground <- understory[understory$type == "g",]
undund <- understory[understory$type == "u",]

## Seedlings data prep:
oak_codes <- c("BO", "RO", "WO", "CO") # oak species for subsetting
# adding columns for number of seedlings per transect and plot:
seedlings <- seedlings_raw %>%
  mutate(plot = paste(hotspot, point, sep = "-"), .after = point) %>%  # add column for plot number
  mutate(all_seed = rowSums(seedlings_raw[,oak_codes])) %>%
  group_by(plot,loc) %>% mutate(plot_seed = sum(all_seed))

## clean up environment:
rm(list=setdiff(ls(), c("plots", "seedlings", "trees", "understory", "oaktrees", "oaks")))


# Update timeline:
# 2024-12-20 seedling and understory begin cleaning up
# 2024-12-19 cleaned and prepped plot, tree, and oak tree data
# 2025_02_14 re-writing for updated and corrected HF data from Audrey
