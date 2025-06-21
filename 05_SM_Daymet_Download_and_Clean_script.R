# Daymet script for getting temperature, precipitation, and vapor pressure defecit data

#### ----- Libraries ----- ####
#install.packages("librarian")
#install.packages("daymetr")
librarian::shelf(tidyverse, dplyr, googledrive, daymetr, lubridate, data.table)

# set wd:
wd <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(wd)

####---- Load Data ----####
# load clean data if not in environment:
coords <- read.csv("Data/sites_5k_samp_coords.csv")
# add sites IDs to coords:
coords <- coords %>%
  select(c(lat, lon)) %>%
  mutate(site = 1:nrow(coords), .before = 1)

plots <- read.csv("Data/plots_clean.csv")
hf_coords <- data.frame(site = 1:nrow(plots), 
                        lat = plots$latitude, 
                        lon = plots$longitude)

# write tables for bulk daymet downloads:
# 5k sample:
write.table(coords, "coords_5k_samp.csv",
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)
# harvard forest sites:
write.table(hf_coords, "coords_hf_samp.csv",
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)


####---- Daymet Download function ----####
#'@param start = starting year of data you want
#'@param end = end year of data you want
#'@param points = csv file where lat/lons are for bulk download
#'@param data = daymet data from download (list)
#'@param vars = daymet variables you want (character vector)
daymet_grab_n_agg <- function(start, end, points, vars){
  # set up conditions for download, years of data and csv file with lat/lon points:
  start = start
  end = end
  file = paste0(getwd(), points)
  
  # download:
  data <- download_daymet_batch(
    file_location = file,
    start = start,
    end = end
  )
  
  # aggregate into monthly means:
  agg_daymet <- function(data, vars){
    agg <- list()
    for (i in 1:length(data)){
      dm <- data[[i]]$data
      dm <- dm %>%
        select(c(year, yday, vars)) %>%  # remote unnecessary data
        mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"), .before = 3) %>%  # get dates for month
        mutate(month = month(date), .before = 4) %>%  # add month 
        group_by(year,month) %>% summarise(across(-c(yday, date), mean)) %>% # compute monthly means for each variable
        mutate(site = i, .before = 1)  # add site ID
      agg[[i]] <- dm
    }
    return(agg)
  }
  
  # identify variables and run:
  vars <- vars # variables I want
  agg <- agg_daymet(data, vars)
  
  # unlist daymet variables for matching to sites:
  dm_data <- as.data.frame(do.call(rbind, agg))
  return(dm_data)
}

vars <- c("prcp..mm.day.", "tmax..deg.c.", "tmin..deg.c.","vp..Pa.")
hf_daym <- daymet_grab_n_agg(2014, 2024, "/coords_hf_samp.csv", vars)
samp_daym <- daymet_grab_n_agg(2014, 2024, "/coords_5k_samp.csv", vars)

# clean up environment:
#rm(coords, hf_coords)

### Update timeline:
# 2025-03-01 created script, wrote functions, tested functions
# 2025-03-01 saved HF daymet data
# 2025-03-06 turned into one function to run automatically

####---- ARCHIVE ----####
# # start with a 10 point sample:
# samp <- coords[1:10,c("site", "lat", "lon")]
# write.table(samp, "dm_samp.csv",
#             sep = ",",
#             col.names = TRUE,
#             row.names = FALSE)

# # daymetr package bulk download call:
# # set start and end dates:
# start = 2014  # time period for data grab
# end = 2017
# #file = paste0(getwd(),"/coords_hf_samp.csv")  # site coordinates
# file = paste0(getwd(), "/coords_5k_samp.csv")

#' # aggregate into monthly means:
#' #'@param data = daymet data from download (list)
#' #'@param vars = daymet variables you want (character vector)
#' agg_daymet <- function(data, vars){
#'   agg <- list()
#'   for (i in 1:length(data)){
#'     dm <- data[[i]]$data
#'     dm <- dm %>%
#'       select(c(year, yday, vars)) %>%  # remote unnecessary data
#'       mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"), .before = 3) %>%  # get dates for month
#'       mutate(month = month(date), .before = 4) %>%  # add month 
#'       group_by(year,month) %>% summarise(across(-c(yday, date), mean)) %>% # compute monthly means for each variable
#'       mutate(site = i, .before = 1)  # add site ID
#'     agg[[i]] <- dm
#'   }
#'   return(agg)
#' }
#' 
#' # identify variables and run:
#' vars <- c("prcp..mm.day.", "tmax..deg.c.", "tmin..deg.c.","vp..Pa.")  # variables I want
#' agg <- agg_daymet(data, vars)
#' 
#' # unlist daymet variables for matching to sites:
#' dm_data <- as.data.frame(do.call(rbind, agg))

