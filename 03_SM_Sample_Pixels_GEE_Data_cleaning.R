# This is the script for Earth Engine product data
# forest condition score and tasseled cap greenness data 

dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)

#### ----- Libraries ----- ####
librarian::shelf(tidyverse, googledrive, ggplot2, RColorBrewer)

#### ----- Load Data ----- ####
# function for loading things from different places:
loader <- function(path, folder, name){
  read.csv(paste0(path, folder, name))
}

# data folders
home <- "Data/"  # where the data lives
field <- "hf_field_data/"  # 2022 field data
gee <- "gee_data/"         # Google Earth Engine data


# load Google Earth Engine samples 5000 points:
#scores_raw <- loader(home, gee, "2025_02_18_5000_points_sample_score_mean.csv")
#tcg_raw <- loader(home, gee, "2025_02_18_5000_points_sample_tcg_mean.csv")
scores_raw <- loader(home, gee, "2025_02_28_growing_season_sample_score_mean.csv")
tcg_raw <- loader(home, gee, "2025_02_28_growing_season_sample_tcg_mean.csv")

# extract the coordinates from the .geo column:
get_lat_lon <- function(data){
  geo <- as.data.frame(data[,".geo"])
  coords<-matrix(nrow=nrow(geo),ncol=2)
  for (i in 1:nrow(geo)){
    #longitudes:
    lon<-str_extract(geo[i,], "\\d+\\.*\\d*")
    coords[i,1]<-as.numeric(lon)*-1
    
    #latitudes:
    extlon<-sub(lon,"",geo[i,])
    coords[i,2]<-as.numeric(str_extract(extlon, "\\d+\\.*\\d*"))
  }
  coords <- as.data.frame(coords)
  colnames(coords)<-c("lon","lat")
  return(coords)
}

coords <- get_lat_lon(scores_raw)

#### ----- Cleaning and Prepping ----- ####
# cleaning up the data!
# condition score data from GEE
coords <- get_lat_lon(scores_raw)
scores <- scores_raw %>%
  select(starts_with("X")) %>%  #select columns with score data
  rename_with(~ str_replace_all(., c("X|_score_mean" = "", "\\." = "-"))) %>%  #name columns with just dates
  mutate(site = 1:nrow(scores_raw), .before = 1) %>%
  mutate(longitude = coords$lon, .before = 2) %>%  # add coordinates
  mutate(latitude = coords$lat, .before = 3)

coords <- get_lat_lon(tcg_raw)
tcg <- tcg_raw %>%
  select(starts_with("X")) %>%
  rename_with(~ str_replace_all(., c("X|_tcg_mean" = "", "\\." = "-"))) %>%
  mutate(site = 1:nrow(tcg_raw), .before = 1) %>%
  mutate(longitude = coords$lon, .before = 2) %>%
  mutate(latitude = coords$lat, .before = 3)

# clean environment:
rm("scores_raw", "tcg_raw")

# time series plots
# just get the time series from the scores and tcg:
scores_ts <- scores[, c(grep("^2", names(scores)))]
tcg_ts <- tcg[,c(grep("^2", names(scores)))]
