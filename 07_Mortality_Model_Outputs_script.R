# Script for visualizing, parsing, and plotting model outputs from mortality models

# load libraries
librarian::shelf(tidyverse, dplyr, rjags, coda)

# navigate to folder:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Mortality_Model_Runs/"
setwd(dir)

### Functions to get model outputs
# get model output:
mort_out <- function(dir, modelnum, log){
  num <- modelnum
  file <- list.files(paste0(dir, "model_outputs/"))
  if(log == TRUE){
    name <-  grep(paste0("modelrun_", as.character(modelnum), "_", "log"), file)
    model <- read.csv(paste0(dir, "model_outputs/", file[name]))
  } else {
    name <-  grep(paste0("modelrun_", as.character(modelnum), "_", "nolog"), file)
    model <- read.csv(paste0(dir, "model_outputs/", file[name]))
  }
  return(model)
}

# get model result:
mort_data <- function(dir, modelnum, log){
  num <- modelnum
  file <- list.files(paste0(dir, "model_runs"))
  if(log == TRUE){
    name <- grep(paste0("modelrun_", as.character(modelnum), "_", "log"), file)
  } else {
    name <- grep(paste0("modelrun_", as.character(modelnum), "_", "nolog"), file)
  }
  load(paste0(dir, "model_runs/", file[name]))
  return(model_info)
}

### Visualizations and parsing outputs:
dir = dir
modelnum = 1
log = TRUE
model <- mort_data(dir, modelnum, log)
output <- mort_out(dir, modelnum, log)

# out <- as.matrix(jags_out)
# pairs(out)
# cor(out)
# gelman.diag(jags_out)
# gelman.plot(jags_out)
# effectiveSize(jags_out)

# make y data for comparison plots
ymeans <- apply(out[,grep("y", colnames(out))], 2, mean)
plot(ymeans, resp$pdba)