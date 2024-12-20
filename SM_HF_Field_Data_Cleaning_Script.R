# This is the script for cleaning Harvard Forest field data from Summer 2022
# Clean up: data organized, demographic calculations for running mortality analyses

#### ----- Load libraries: ----- ####
#install.packages("librarian")
librarian::shelf(tidyverse, googledrive)

#### ----- Load data: ----- ####
trees_raw <- read.csv("Data/HF_Field_Data/Tree_data.csv")
seedlings_raw <- read.csv("Data/HF_Field_Data/Seedlings_long.csv")
understory_raw <- read.csv("Data/HF_Field_Data/Und_ground_survey.csv")
plots_raw <- read.csv("Data/HF_Field_Data/Plot_data.csv")

#### ----- Clean and prep data: ----- ####
## Plot data:
plots <- plots_raw %>%
  mutate(plot = str_replace(plot, " ", "-")) %>%  # remove space between plot number
  mutate(hotspot = as.numeric(substr(plot, 1, 1))) %>%  # add hot spot identifier column
  mutate(latitude = str_replace(latitude, " N", "")) %>%  # clean lat/long columns
  mutate(longitude = str_replace(longitude, " W", "")) %>%
  mutate_at(c('latitude','longitude'), as.numeric) %>%
  mutate(longitude = longitude * -1) %>%  # correct longitude now that numeric to make W
  mutate(invasives = ifelse(invasives == "no", 0, 1)) %>%  # change invasives from "y/n" to "0/1"
  mutate(recent_timber_harvest = ifelse(recent_timber_harvest == "no", 0, 1))  #change timber harvest from "y/n" to "0/1"

## Trees data:
BAF = 10  # basal area factor from variable radius plots
trees <- trees_raw %>%
  mutate(plot = str_replace(plot, " ", "-")) %>%  # clean plot names
  mutate(hotspot = as.numeric(substr(plot, 1, 1))) %>%  # add hotspot column
  mutate(ba = (dbh^2) * 0.005454) %>% # basal area calculation from dbh
  group_by(hotspot,plot) %>% mutate(treeper  = n()) %>%  # number of trees per plot
  group_by(hotspot) %>% mutate(hotcount = length(unique(plot))) %>%  # number of plots per hotspot
  group_by(plot) %>% mutate(tba = sum(ba)) %>%
  mutate(exp_fac = BAF / ba / hotcount) %>% # expansion factor, BAF/BA/number in plot
  group_by(hotspot) %>% mutate(TPA = sum(exp_fac))

## Oak-specific data:
oaks <- c("BO", "RO", "WO", "CO") # oak species for subsetting

percent <- function(x, count){  # percentage function for calculating ba and dead tree percentages
  (x / count)*100
}

# make prep oak tree data
oaktrees <- subset(trees, spp %in% oaks)
oaktrees <- oaktrees %>% 
  group_by(plot) %>% mutate(tot_oak = n(),  # number of oaks
                            dead = ifelse(Cond == "L", 0, 1),
                            tba_oak = sum(ba),  # total pak basal area per plot
                            pba = percent(tba_oak,tba),  # percent oak basal area in plot
                            deadoak = length(which(Cond == "D")),  # total dead oaks
                            pdead = percent(deadoak,treeper),  # percent dead oak of all trees in plot
                            pdo = percent(deadoak,tot_oak)) %>%  # percent dead oak of all oaks in plot
  group_by(plot,Cond) %>% mutate(dba = sum(ba)*dead,
                                 pdba = percent(dba,tba)) %>%
  group_by(hotspot) %>% mutate(tba_oak_hot = sum(ba),  # total oak basal area in hotspot
                               pba_hot = percent(tba_oak,tba),  # percent oak basal area in hotspot
                               deadoak_hot = sum(deadoak),  # sum dead oak in hotspot
                               pdo_hot = percent(deadoak,tot_oak),  # percent dead oak in hotspot
                               dba_hot = sum(dba),  # total dead basal area in hot spot
                               pdba_hot = percent(dba_hot,sum(tba))) # percent dead basal area in hotspot


## Understory data prep:
understory <- understory_raw %>%
  mutate(plot = str_replace(plot, " ", "-")) %>% 
  mutate(hotspot = as.numeric(substr(plot, 1, 1)))

undground <- understory[understory$type == "g",]
undund <- understory[understory$type == "u",]

## Seedlings data prep:
seedlings <- seedlings_raw %>%
  mutate(plot = str_replace(plot, " ", "-")) %>% 
  mutate(hotspot = as.numeric(substr(plot, 1, 1))) %>%
  mutate(all_seed = rowSums(seedlings[,oaks]))

# Update timeline:
# 2024-12-20 seedling and understory begin cleaning up
# 2024-12-19 cleaned and prepped plot, tree, and oak tree data
