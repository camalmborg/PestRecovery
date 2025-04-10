# Script for visualizing, parsing, and plotting model outputs from univariate mortality models

# load libraries
librarian::shelf(tidyverse, dplyr, rjags, coda)

# navigate to folder:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Mortality_Model_Runs/"
setwd(dir)

### Functions to get model outputs
# get model output:
#'@param dir = folder where the model_outputs and model_runs are located
#'@param modelnum = model number for univariate run
#'@param log = logical; TRUE = logit model run, FALSE = non-logit model run
mort_out <- function(dir, modelnum, log){
  num <- modelnum
  file <- list.files(paste0(dir, "model_outputs/"))
  if(log == TRUE){
    name <-  grep(paste0("modelrun_", as.character(modelnum), "_", "log"), file)
  } else {
    name <-  grep(paste0("modelrun_", as.character(modelnum), "_", "nolog"), file)
  }
  model <- read.csv(paste0(dir, "model_outputs/", file[name]))
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

# ### Visualizations and parsing outputs:
# # selecting model:
# dir = dir
# modelnum = 1
# log = TRUE
# # load model and output:
# model <- mort_data(dir, modelnum, log)
# output <- mort_out(dir, modelnum, log)
# 
# # visual inspections:
# jags_out <- model$jags_out
# vars <- varnames(jags_out)
# params <- jags_out[,grep("^alpha|b|q|tau", vars)]
# 
# gelman.diag(params)
# gelman.plot(params)
# effectiveSize(params)
# # visual inspect notes:
# # model 1-log: covergence good, burn in 20,000, effective sizes all >3900
# # model 2-log:
# 
# # making pred/obs plots:
# pred_obs_plot <- function(output, model){
#   predict <- apply(output[,grep("y", colnames(output))], 2, mean)
#   obs <- model$metadata$data$y
#   plot(predict, obs)
# }
# 
# # predict <- apply(output[,grep("y", colnames(output))], 2, mean)
# # obs <- model$metadata$data$y
# # plot(predict, obs)
# 
# # getting model dics:
# dic <- model$dic[[2]]

### loop for making pred-obs plots and recording dics
# number of models:
model_count <- length(list.files(paste0(dir, "model_outputs")))
model_numbers <- rep(1:10, 2)
log_types <- c(rep(TRUE, 10), rep(FALSE, 10))
# make empty list for filling in dics:
mort_dics <- matrix(NA, nrow = model_count, ncol = 3)
# the loop:
for (i in 1:model_count){
  # selecting model:
  dir = dir
  modelnum = model_numbers[i]
  log = log_types[i]
  # load model and output:
  model <- mort_data(dir, modelnum, log)
  output <- mort_out(dir, modelnum, log)
  # collect dic and model identity:
  mort_dics[i,1] <- modelnum
  mort_dics[i, 2] <- log
  mort_dics[i, 3] <- model$dic[[2]]
  # do a pred/obs plot:
  pred_obs_plot(output, model)
  print(i)
}

# make y data for comparison plots
# ymeans <- apply(out[,grep("y", colnames(out))], 2, mean)
# plot(ymeans, resp$pdba)

### Archive:
# out <- as.matrix(jags_out)
# pairs(out)
# cor(out)
# gelman.diag(jags_out)
# gelman.plot(jags_out)
# effectiveSize(jags_out)