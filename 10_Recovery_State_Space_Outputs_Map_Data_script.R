### Recovery Rate Across Sites script

## Load libraries and necessary environments
librarian::shelf(dplyr, tidyverse, rjags, coda)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

## pull in model output files if not in the environment already
model_params <- read.csv("2025_07_31_all_base_uni_recov_models_param_means.csv")
load("2025_07_31_recov_models_outputs_list.RData")  # object is called model_outputs
# model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]

# choose model:
m_num <- 8  # change this when changing models
model_pick <- models[m_num]
# load model_info object
load(paste0(dir, "model_runs/", model_pick))
model_inputs <- model_info$metadata$model_data
# model parameters:
out <- as.matrix(model_info$jags_out)
# separate out specific params:
x_params <- grep("^x", colnames(out))
beta_params <- grep("^b",colnames(out))
taus <- grep("tau", colnames(out))
r <- grep("r0", colnames(out))
big_r <- grep("R", colnames(out))
# from model inputs:
y <- as.matrix(model_inputs$y)

## Getting recovery rates around region
# get R's from model results:
recov_rates <- out[,big_r]
# get means:
rr <- as.data.frame(apply(recov_rates, 2, mean))
# organize by rowname:
recov_rates <- matrix(NA, nrow = 5000, ncol = 5)
colnames(recov_rates) <- colnames(y)[-1]
for (i in 2:6){
  rrs <- grep(paste0(",", i, "\\]"), rownames(rr))
  recov_rates[,i-1] <- rr[rrs,]
}
rm(rrs)
# get mean over time:
r_means <- apply(recov_rates, 1, mean)

# join with coordinates:
reg_recov <- data.frame(lon = coords$lon,
                        lat = coords$lat,
                        recov = r_means)

## Make a map
# load terra library:
librarian::shelf(terra, sf, tigris, ggplot2)

# make it into a point vector:
recov_vec <- vect(reg_recov, geom = c("lon", "lat"), crs = "EPSG:4326")

# get states:
states <- states(cb = TRUE) %>%
  filter(NAME %in% c("Massachusetts", "Connecticut", "Rhode Island"))
states <- vect(states)
states <- project(states, recov_vec)

# test plots:
plot(recov_vec)
plot(states, add = T)

# convert both to sf for final plots:
recov_vec <- st_as_sf(recov_vec)
states <- st_as_sf(states)

# make a nice ggmap:
recov_map <- ggplot(recov_vec) +
  # add the points:
  geom_sf(aes(fill = r_means), size = 2,
          color = "black", shape = 21, stroke = 0.1) +
  # change colors:
  scale_fill_gradient(low = "orange", high = "limegreen") +
  # add the state outlines:
  geom_sf(data = states, fill = NA, color = "black", size = 0.5) +
  # add labels
  labs(title = "Mean Annual Recovery Rates",
       fill = "Predicted Recovery Rates") +
  theme_bw() +
  theme(panel.grid = element_blank())

recov_map


## Make an inset map:
usa <- tigris::states(cb = TRUE, resolution = "20m") %>%
  filter(!STUSPS %in% c("AK", "HI", "PR", "GU", "VI", "AS", "MP"))

# bounding box of data:
bbox <- st_bbox(recov_vec)

# Create an inset map with a rectangle showing the extent of your main map
inset_map <- ggplot() +
  geom_sf(data = usa) +
  annotate("rect",
           xmin = bbox$xmin, xmax = bbox$xmax,
           ymin = bbox$ymin, ymax = bbox$ymax,
           fill = NA, color = "red", linewidth = 0.7) +
  theme_void()

inset_map
