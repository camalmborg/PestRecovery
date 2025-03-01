# Daymet script for getting temperature, precipitation, and vapor pressure defecit data

#### ----- Libraries ----- ####
#install.packages("librarian")
#install.packages("daymetr")
librarian::shelf(tidyverse, dplyr, googledrive, daymetr, lubridate)

####---- Load Data ----####
# load clean data if not in environment:
# coords <- read.csv("sites_5k_samp_coords.csv)
# plots <- read.csv("plots_clean.csv)

# add sites IDs to coords:
coords <- coords %>%
  mutate(site = 1:nrow(coords), .before = 1)

####---- Daymet Download ----####
# start with a 5 point sample:
samp <- coords[1:10,c("site", "lat", "lon")]
write.table(samp, "dm_samp.csv",
            sep = ",",
            col.names = TRUE,
            row.names = FALSE)

# daymetr package bulk download call:
test <- download_daymet_batch(
  file_location = paste0(getwd(),"/dm_samp.csv"),
  start = 2010,
  end = 2024
)
# ^ this works

# aggregate into monthly means:
agg_daymet <- function(data, vars){
  agg <- list()
  for (i in 1:length(data)){
    dm <- data[[i]]$data
    dm <- dm %>%
      select(c(year, yday, vars)) %>%  # remote unnecessary data
      mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"), .before = 3) %>%  # get dates for month
      mutate(month = month(date), .before = 4) %>%  # add month 
      group_by(month) %>% summarise(across(-c(year, yday, date), mean))  # compute monthly means for each variable
    agg[[i]] <- dm
  }
  return(agg)
}

# testing...
vars <- c("prcp..mm.day.", "tmax..deg.c.", "tmin..deg.c.","vp..Pa.")  # variables I want
agg <- agg_daymet(test, vars)

# unlist the daymet variables for matching to sites:
unlist_daymet <- function(agg){
  for (i in length(agg)){
    dm <- agg[[i]]
  }
}

### Update timeline:
# 2025-03-01 created script