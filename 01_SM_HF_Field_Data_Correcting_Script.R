# Script for making asset to sample condition scores and tasseled cap greenness values
# from Earth Engine asset
# created Feb 4, 2025
# includes updated lat/long points corrected from HF REU 2022 data

# load libraries:
library(dplyr)
library(tidyverse)

# open plot file:
file <- "Data/hf_field_data/2022_hf_plot_data.csv"
plots_cor <- read.csv(file)

# select plots that are from 2022 data collection:
plots_2022 <- plots_cor[which(!is.na(plots_cor$latitude_2022)),]
#length(complete.cases(plots_2022$latitude_2022)) == nrow(plots_2022) # check: if TRUE, all lat/lon included
sites <- plots_2022[,c("point", "latitude_2022", "longitude_2022")] # make sites table for GEE
colnames(sites) <- c("point", "latitude", "longitude")
sites$site <- 1:nrow(plots_2022) # add site number column

# save sites table as csv for GEE asset upload
write.csv(sites, "Data/hf_field_data/hf_2022_sites.csv", row.names = FALSE)
