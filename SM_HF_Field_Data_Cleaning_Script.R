# This is the script for cleaning Harvard Forest field data from Summer 2022
# Clean up: data organized, demographic calculations for running mortality analyses

### ----- Load libraries: ----- ###
#install.packages("librarian")
librarian::shelf(tidyverse, googledrive)

### ----- Load data: ----- ###
trees <- read.csv("Data/HF_Field_Data/Tree_data.csv")
seedlings <- read.csv("Data/HF_Field_Data/Seedlings_long.csv")
understory <- read.csv("Data/HF_Field_Data/Und_ground_survey.csv")
plots <- read.csv("Data/HF_Field_Data/Plot_data.csv")

### ----- Clean and prep data: ----- ###
## Plot data:
plots <- plots %>%
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
trees <- trees %>%
  mutate(plot = str_replace(plot, " ", "-")) %>% 
  mutate(hotspot = as.numeric(substr(plot, 1, 1))) %>%
  mutate(ba = (dbh^2) * 0.005454) %>% # basal area calculation from dbh  #%>%
  group_by(hotspot,plot) %>% mutate(treeper  = n()) %>%  # number of trees per plot
  group_by(hotspot) %>% mutate(hotcount = length(unique(plot))) %>%  # number of plots per hotspot
  mutate(exp_fac = BAF / ba / hotcount) %>% # expansion factor, BAF/BA/number in plot
  group_by(hotspot) %>% mutate(TPA = sum(exp_fac))

# Update timeline:
# 2022-12-19