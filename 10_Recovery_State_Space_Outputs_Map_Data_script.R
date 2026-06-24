### Recovery Rate Across Sites script

## Load libraries and necessary environments
librarian::shelf(dplyr, tidyverse, rjags, coda)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)
# load model performance information:
dic_sort <- read.csv("2025_11_30_all_recov_models_dics.csv", row.names = 1)

## Load Best Model
# all model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]
# best model:
best <- dic_sort$model_number[1]
best_model <- models[best] 
# load model:
load(paste0(dir, "model_runs/", best_model))
# load model_params:
model_params <- read.csv(file = "2025_11_30_all_recov_models_param_means.csv")

## Get Model Information and Data
# get parameter values:
best_params <- model_params[best,]
name = best_model
# model name:
model_name <- best_model
if (grepl("cov", name) == TRUE){
  model_name <- str_match(name, "cov_(.*?)_data")[,2]
} else {
  model_name <- str_match(name, "_(.*?)_model")[,2]
}
# load model inputs:
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
rr <- apply(recov_rates, 2, mean)
# organize by rowname:
recov_rates_mx <- matrix(NA, nrow = 5000, ncol = 6)
colnames(recov_rates_mx) <- colnames(y)[-1]
for (i in 2:7){
  rrs <- grep(paste0(",", i, "\\]"), names(rr))
  recov_rates_mx[,i-1] <- rr[rrs]
}
rm(rrs, rr)

# get mean over time:
r_means <- apply(recov_rates_mx, 1, mean)
# recovery anomaly:
recov_anom <- r_means - mean(r_means)

# forecast recovery rates:
## Load forecast result
forecast <- read.csv("Recovery_Forecasts/2025-12-03_ens_1500_model_1_start_year_2017_reforecast_result.csv")
# prepare for residual calculation:
y_pred <- forecast %>%
  # add baseline to predictions:
  mutate(across(!site, ~ tcg_base$baseline - .x)) %>%
  # ensemble means:
  mutate(site = as.factor(site)) %>%
  group_by(site) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  select(-c(site))
# recovery slopes:
recov_slopes <- c()
time <- 1:7
for (i in 1:nrow(y_pred)){
  recov_y <- unlist(y_pred[i,])
  slope <- lm(recov_y ~ time)
  recov_slopes[i] <- slope$coefficients[2]
}
# recovery anomaly:
recov_anom <- recov_slopes - mean(recov_slopes, na.rm = T)


# join with coordinates:
reg_recov <- data.frame(lon = coords$lon,
                        lat = coords$lat,
                        #recov = r_means
                        recov_anom)[-which(is.na(recov_anom)),]

# attempting to make recovery rate observed at each point 06/24/2026
# forecast time period:
years <- 2017:2023
# observations:
tcg <- read.csv("/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Data/tcg_5ksamp_clean.csv")[-1] %>%
  # rename:
  rename_with(~ str_replace_all(.x, c("^\\s*X" = "", "\\." = "-"))) %>%
  # select years:
  select(`2017-05-01`:`2023-05-01`) %>%
  # rename with just years:
  rename_with(~ sub("-.*", "", .x))

# find recovery rate for each
# make slopes dataframe with lat lon:
slopes <- as.matrix(data.frame(site = 1:nrow(tcg), 
                               lon = coords$lon, 
                               lat = coords$lat,
                               slopes = NA))
# find line between 2017-2023:
for(i in 1:nrow(tcg)){
  row <- tcg[i,] |>
    pivot_longer(cols = everything(),
                 names_to = "year",
                 values_to = "value") |>
    mutate(year = as.numeric(year))
  slope <- coef(lm(value ~ year, data = row))[["year"]]
  slopes[i,4] <- slope
}
# make data frame again:
slopes <- as.data.frame(slopes)

## Make a map
# load terra library:
librarian::shelf(terra, sf, tigris, ggplot2)

# make it into a point vector:
recov_vec <- vect(slopes, geom = c("lon", "lat"), crs = "EPSG:4326")

# get states:
states <- states(cb = TRUE) %>%
  filter(NAME %in% c("Massachusetts", "Connecticut", "Rhode Island"))
states <- vect(states)
states <- project(states, recov_vec)

# state bounding box:
bbox <- st_bbox(states)

# all US states for background
us_states <- states(cb = TRUE) |>
  st_as_sf() |>
  st_transform(st_crs(recov_vec))

# test plots:
plot(recov_vec)
plot(states, add = T)

# convert both to sf for final plots:
recov_vec <- st_as_sf(recov_vec)
states <- st_as_sf(states)

# make a nice ggmap:
recov_map <- ggplot(recov_vec) +
  # add the US background
  geom_sf(data = us_states, fill = "grey95", color = "grey40", linewidth = 0.2) +
  # add the state outlines:
  geom_sf(data = states, fill = "grey70") +
  # add the points:
  geom_sf(aes(fill = recov_anom), size = 1.75,
          color = "black", shape = 21, stroke = 0.1) +
  # change colors:
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red", midpoint = 0) +
  # add the state outlines:
  geom_sf(data = states, fill = NA, color = "black", size = 0.5) +
  # crop to MA/CT/RI:
  coord_sf(xlim = c(bbox["xmin"] - 0.2, bbox["xmax"] + 0.2),
           ylim = c(bbox["ymin"] - 0.2, bbox["ymax"] + 0.2),
           expand = FALSE) +
  # add labels
  labs(fill = "Observed Recovery\nRate Slope (TCG/year)\nAnomaly from Mean") +
  theme_bw() +
  theme(panel.grid = element_line(linetype = "dashed"),
        panel.background = element_rect(fill = "lightsteelblue1"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.justification = c(0.5, 0))

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
           fill = NA, color = "red", linewidth = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank())

inset_map

library(cowplot)

# Convert inset map to a grob object
inset_grob <- ggplotGrob(inset_map)

# Combine the maps
final_map <- ggdraw(recov_map) +
  draw_plot(inset_grob, x = 0.68, y = 0.58, 
            width = 0.3, height = 0.3) # Adjust x, y, width, height for placement and size

print(final_map)


# save them:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/"
setwd(save_dir)
# Save the plot to a PNG file
png(filename = "2026_06_24_recov_rate_anom_map_revised.png",
    width = 10, height = 8, units = "in",
    res = 600)
final_map
dev.off()

# ggsave("2026_03_08_recov_rate_map_ESA.png", 
#        plot = final_map,
#        dpi = 600)



# join with coordinates:
disturb <- tcg[,1]
disturb_anom <- tcg[,1] - mean(tcg[,1], na.rm = T)

reg_disturb <- data.frame(lon = coords$lon,
                        lat = coords$lat,
                        disturb = disturb,
                        disturb_anom = disturb_anom)
# make it into a point vector:
disturb_vec <- vect(reg_disturb, geom = c("lon", "lat"), crs = "EPSG:4326")
# convert to sf for final plots:
disturb_vec <- st_as_sf(disturb_vec)


# make a nice ggmap:
disturb_2017_map <- ggplot(disturb_vec) +
  # add the state outlines:
  geom_sf(data = states, fill = "grey70") +
  # add the points:
  geom_sf(aes(fill = disturb), size = 1.75,
          color = "black", shape = 21, stroke = 0.1) +
  # change colors:
  scale_fill_gradient2(low = "dodgerblue", mid = "white", high = "red", 
                       midpoint = median(disturb_vec$disturb, na.rm = T)) +
  # add the state outlines:
  geom_sf(data = states, fill = NA, color = "black", size = 0.5) +
  # add labels
  labs(fill = "TCG 2017\nAnomaly from Mean") +
  theme_bw() +
  theme(panel.grid = element_line(linetype = "dashed"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.justification = c(0.5, 0))

disturb_2017_map
