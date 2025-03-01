### This is the script for Harvard Forest 2022 field site disturbance magnitude and disturbance occurrance
### The script calculates a disturbance magnitude for the years of a disturbance, in this case 2016-2017 (into 2018)

#### ----- Libraries ----- ####
#install.packages("librarian")
librarian::shelf(tidyverse, dplyr, googledrive)

#### ----- Load and set up data ----- ####
# if not already in environment:
#hf_scores <- read.csv("Data/hf_plot_scores_clean.csv") # these come in with wrong variable names for the time series
#hf_tcg <- read.csv("Data/hf_plot_tcg_clean.csv")

# make the calculator function:
#'@param ts = time series data for detecting disturbance in canopy greenness
#'@param distyr = onset year of disturbance
dist_mag_calc <- function(ts, distyr, scrs){
  # separate columns with just the canopy observation data:
  series <- ts[,grep("^2", names(ts))]
  series$site <- 1:nrow(series)
  # get disturbance onset column:
  dist <- grep(as.character(distyr), colnames(series))
  # calculate a pre-disturbance steady state to compare disturbance condition:
  steady <- apply(series[,(dist-6):(dist-1)], 1, mean, na.rm = T)
  # calculate disturbance magnitudes:
  dmagy1 <- steady - series[,dist]      # onset year
  dmagy2 <- steady - series[,dist+1]    # following year
  # disturbance magnitudes being > 0 is defol:
  def1 <- dmagy1 > 0
  def2 <- dmagy2 > 0
  # combine dataset for output
  dmag_data <- cbind(series[,'site'], series[,dist], series[,dist+1], steady, dmagy1, dmagy2, def1, def2)
  colnames(dmag_data) <- c("site", names(series[dist]), names(series[dist+1]), "steady", "dmag1", "dmag2", "dm_gt_0_y1", "dm_gt_0_y2")
  return(dmag_data)
}

# Update timeline
# 2025-02-20 created script

####----ARCHIVE----####
# # calculating disturbance thresholds:
# sc_series <- scrs[,grep("^2", names(scrs))]  # use condition scores to identify sites @ threshold
# d <- quantile(sc_series[,(dist-6):(dist-1)], c(0.05), na.rm = T)  # 5% quantile in 5yr pre-onset
# d1 <- scrs[,dist] < d
# d2 <- scrs[,dist+1] < d