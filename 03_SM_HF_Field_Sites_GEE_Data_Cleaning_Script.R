# This is the script for HF field data sites from Earth Engine product
# forest condition score and tasseled cap greenness data 

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

# load clean plot data if not in environment:
plots <- read.csv(paste0(home,"plots_clean.csv"))

# load Google Earth Engine samples from field sites:
scores_raw <- loader(home, gee, "2025_02_18_hf_sites_sample_score_mean.csv")
tcg_raw <- loader(home, gee, "2025_02_18_hf_sites_sample_tcg_mean.csv")

#### ----- Cleaning and Prepping ----- ####
# cleaning up the data!
# condition score data from GEE
scores <- scores_raw %>%
  arrange(site) %>%  # puts them in site order
  select(starts_with("X")) %>%  #select columns with score data
  rename_with(~ str_replace_all(., c("X|_score_mean" = "", "\\." = "-"))) %>%
  mutate(hotspot = plots$hotspot, .before = 1) %>%  # add hotspot column
  mutate(point = plots$point, .before = 2) %>%  # add point (plot number) column
  mutate(plot = plots$plot, .before = 3) %>%  # add plot name column
  mutate(latitude = plots$latitude) %>%  # add lat and lon
  mutate(longitude = plots$longitude)

# tasseled cap greenness data from GEE
# same as above with tcg data
tcg <- tcg_raw %>%
  arrange(site) %>%
  select(starts_with("X")) %>%
  rename_with(~ str_replace_all(., c("X|_tcg_mean" = "","\\." = "-"))) %>%
  mutate(hotspot = plots$hotspot, .before = 1) %>%
  mutate(point = plots$point, .before = 2) %>%
  mutate(plot = plots$plot, .before = 3) %>%
  mutate(latitude = plots$latitude) %>%
  mutate(longitude = plots$longitude)

# make these their own things to define from 5000 sample points:
hf_scores <- scores
hf_tcg <- tcg

#### ----- Plots ----- ####
# time series plots
# just get the time series from the scores and tcg:
scores_ts <- scores[, c(grep("^2", names(scores)))]
tcg_ts <- tcg[,c(grep("^2", names(scores)))]

# get means from hotspots:
hotspot_means <- scores %>%
  pivot_longer(cols = c(grep("^2", names(scores)))) %>%
  group_by(name, hotspot) %>% summarise(mean = mean(value, na.rm = T), 
                                        lower = quantile(value, 0.05, na.rm = T),
                                        upper = quantile(value, 0.95, na.rm = T),
                                        .groups = "drop") %>%  #this one works
  mutate(date = as.Date(name)) %>%
  mutate(hotspot = as.factor(hotspot)) %>%
  arrange(hotspot) %>%
  select(hotspot, mean, lower, upper, date)
  # pivot_wider(
  #   names_from = name,
  #   values_from = c(mean, lower, upper)
# )

# alternative for including all data long format:
hotspot_means_full <- scores %>%
  pivot_longer(cols = c(grep("^2", names(scores)))) %>%
  group_by(name, hotspot) %>% mutate(mean = mean(value, na.rm = T),
                                     lower = quantile(value, 0.05, na.rm = T),
                                     upper = quantile(value, 0.95, na.rm = T))


### Time series plot for mean condition scores by hotspot ###
# make plot themes
# color palette:
colors <- brewer.pal(n = 6, name = "Set2")
# theme:
ts_theme <- theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_text(size = 13),
        axis.line = element_line(colour = "black"))

# make the plot:
time_series <- ggplot(data = hotspot_means, aes(x = date, y = mean, color = hotspot)) +
  geom_line(aes(group = hotspot), linewidth = 0.5, show.legend = FALSE,
            position = position_dodge(width = 0.2)) +
  scale_color_manual(values = colors) 
                      #+
                        # geom_line(aes(group = hotspot), linewidth = 0.5, show.legend = F,
                        #           position = position_dodge(width = 0.2)))


# Update timeline
# 2025-02-20 started making time series plots of hot spot mean values (set up for cs or tcg)
# 2025-02-20 starting to clean and prepare HF GEE data - updating file grabbing and csvs
# 2024-12-20 created; adding GEE data

#### ----- Archive ----- ####
# # load lat/lon data:
# hf_lat_lon_raw <- loader(home, field, "hf_2022_sites.csv")
# # load plot data:
# plots_raw <- loader(home, field, "2022_hf_plot_data.csv")
# plots <- plots_raw[which(!is.na(plots_cor$latitude_2022)),]  # grab just the 2022 field sites
# # lat/lon data
# hf_lat_lon <- hf_lat_lon_raw %>%
#   mutate(plot = str_replace(hotplot, "_", "-")) %>%  #make matching "plot" column for matching tree, plots, seedling, and understory data
#   mutate(hotspot = as.numeric(substr(plot, 1, 1)))  #make hotspot identifier column

# # objects with column means from cs and tcg data:
# score_means <- apply(scores, 2, mean, na.rm = T)
# tcg_means <- apply(tcg, 2, mean, na.rm = T)
# # dates object
# years <- names(scores) %>%
#   str_replace_all(.,"-.*", "")  #remove everything after year 
#   #str_replace_all(., "(?<=6).*", "")  #replaces everything after 6
# dates <- paste0("June ", years)