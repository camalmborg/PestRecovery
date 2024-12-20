# This is the script for HF field data sites from Earth Engine product
# forest condition score and tasseled cap greenness data 

#### ----- Libraries ----- ####
#install.packages("librarian")
librarian::shelf(tidyverse, googledrive)

#### ----- Load Data ----- ####
loader <- function(path, name){
  read.csv(paste0(path, name))
}
file <- "Data/GEE_Data/HF_Field_Sites/"
hf_lat_long <- loader(file, "hf_lat_lon.csv")
scores <- loader(file, "HF_sample_score_mean.csv")
tcg <- loader(file, "HF_sample_score_mean.csv")



# Update timeline
# 2024-12-20 created; adding GEE data
