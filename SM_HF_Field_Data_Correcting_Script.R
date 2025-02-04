# Script for making asset to sample condition scores and tasseled cap greenness values
# from Earth Engine asset
# created Feb 4, 2025
# includes updated lat/long points corrected from HF REU 2022 data

# load libraries:
library(dplyr)
library(tidyverse)

# open plot file:
file <- "Data/HF_Field_Data_fixed/plots20172022.csv"
plots_cor <- read.csv(file)

# select plots that are from 2022 data collection:
plots_2022 <- plots_cor[which(!is.na(plots_cor$latitude.2022)),]
#length(complete.cases(plots_2022$latitude.2022)) == nrow(plots_2022) # check: if TRUE, all lat/lon included
sites <- plots_2022[,c("point", "latitude.2022", "longitude.2022")] # make sites table for GEE
colnames(sites) <- c("point", "latitude", "longitude")
sites$site <- 1:nrow(plots_2022) # add site number column

# save sites table as csv for GEE asset upload
write.csv(sites, "Data/HF_Field_Data_fixed/hf_2022_sites.csv")
