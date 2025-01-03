# This is the script for HF field data sites from Earth Engine product
# forest condition score and tasseled cap greenness data 

#### ----- Libraries ----- ####
#install.packages("librarian")
librarian::shelf(tidyverse, googledrive, ggplot2)

#### ----- Load Data ----- ####
loader <- function(path, name){
  read.csv(paste0(path, name))
}
file <- "Data/GEE_Data/HF_Field_Sites/"
hf_lat_lon_raw <- loader(file, "hf_lat_lon_corrected.csv")
scores_raw <- loader(file, "2025_01_02_HF_sample_score_mean.csv")
tcg_raw <- loader(file, "2025_01_02_HF_sample_score_mean.csv")

#### ----- Cleaning and Prepping ----- ####
# clearning up data
# lat/lon data
hf_lat_lon <- hf_lat_lon_raw %>%
  mutate(plot = str_replace(hotplot, "_", "-")) %>%  #make matching "plot" column for matching tree, plots, seedling, and understory data
  mutate(hotspot = as.numeric(substr(plot, 1, 1)))  #make hotspot identifier column

# condition score data from GEE
scores <- scores_raw %>%
  select(starts_with("X")) %>%  #select columns with score data
  rename_with(~ str_replace_all(., c("X|_score_mean" = "", 
                                     "\\." = "-")))
# tasseled cap greenness data from GEE
tcg <- tcg_raw %>%
  select(starts_with("X")) %>%
  rename_with(~ str_replace_all(., c("X|_score_mean" = "",
                                     "\\." = "-")))
# objects with column means from cs and tcg data:
score_means <- apply(scores, 2, mean, na.rm = T)
tcg_means <- apply(tcg, 2, mean, na.rm = T)
# dates object
years <- names(scores) %>%
  str_replace_all(.,"-.*", "")  #remove everything after year 
  #str_replace_all(., "(?<=6).*", "")  #replaces everything after 6
dates <- paste0("June ", years)

#### ----- Plots ----- ####
# time series plots
time_series <- ggplot()

# Update timeline
# 2024-12-20 created; adding GEE data
