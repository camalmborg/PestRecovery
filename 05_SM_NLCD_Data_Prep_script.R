# National Land Cover Database data preparation:

### load libraries
librarian::shelf(dplyr, terra)

### load raw data
# directory:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/nlcd"
setwd(dir)
# load data:
nlcd <- read.csv("sample_nlcd.csv") %>%
  select(lat, lon, landcover, impervious, percent_tree_cover)

# legend landcover types:
legend <- read.csv("NLCD_landcover_legend.csv", header = T, na.strings=c("","NA")) %>%
  filter(!is.na(Legend))




# coordinates:
# coords object loaded from 03_SM_Sample_Pixels_GEE_Data_cleaning

# # load tiffs:
# tiffs <- list.files(dir)[grep("tiff", list.files(dir))][1:3]
# # note:
# # tiffs[1] = fractional cover
# # tiffs[2] = 
# # tiffs[3] = land cover
# 
# ### Processing tiff files
# # load tiffs:
# nlcd_frac_cover <- terra::rast(tiffs[1])
# nlcd_land_conf <- terra::rast(tiffs[2])
# nlcd_land_cov <- terra::rast(tiffs[3])
# # change projections:
# nlcd_frac_cover <- terra::project(nlcd_frac_cover, "epsg:4326")
# nlcd_land_conf <- terra::project(nlcd_land_conf, "epsg:4326")
# nlcd_land_cov <- terra::project(nlcd_land_cov, "epsg:4326")
# # convert sites to vector:
# sites <- terra::vect(coords[,2:3])
# # extract
# nlcd_extract_frac_cov <- terra::extract(nlcd_frac_cover, coords[,2:3])
# nlcd_extract_land_conf <- terra::extract(nlcd_land_conf, coords[,2:3])
# nlcd_extract_land_cov <- terra::extract(nlcd_land_cov, coords[,2:3])
# # compile:
# nlcd_data <- data.frame(lat = coords$lat, lon = coords$lon)
# nlcd_data["frac_cover"] <- nlcd_extract_frac_cov$Annual_NLCD_FctImp_2016
# nlcd_data["land_cov"] <- nlcd_extract_land_cov$Annual_NLCD_LndCov_2016
