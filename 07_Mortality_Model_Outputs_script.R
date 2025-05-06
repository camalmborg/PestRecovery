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

mort_data_load <- function(dir, model){
  load(paste0(dir, "model_runs/", model))
  return(model_info)
}


### loop for making pred-obs plots and recording dics
# number of models:
models <- list.files(paste0(dir, "model_runs"))
model_list <- models[grep("log", model_list)]
# make empty list for filling in dics:
mort_dics <- matrix(NA, nrow = length(model_list), ncol = 4)
#colnames(mort_dics) <- c("id", "model", "DIC")
colnames(mort_dics) <- c("id", "model", "log T/F", "DIC")
# the loop:
for (i in 1:length(model_list)){
  # selecting model:
  dir = dir
  # load model and output:
  model <- mort_data_load(dir, model_list[i])
  # collect dic and model identity:
  mort_dics[i, 1] <- i
  mort_dics[i, 2] <- model$metadata$run
  #mort_dics[i, 2] <- modelnum
  mort_dics[i, 3] <- model$metadata$mod
  #mort_dics[i, 3] <- log
  mort_dics[i, 4] <- model$dic[[2]]
  # track loop
  print(i)
}
mort_dic <- as.data.frame(mort_dics)
mort_dic$delDIC <- as.numeric(mort_dic$DIC) - min(as.numeric(mort_dic$DIC))
mort_dic <- mort_dic[order(mort_dic$delDIC),]
mort_dic$perform <- 1:nrow(mort_dic)
#write.csv(mort_dic, file = "2025_04_23_mortality_models_dic_table.csv")

### Making results tables:
#'@param dir = directory where model results are found
#'@param model_list = character vector with list of model result file names/RData objects
results_extract <- function(dir, model_list){
  # full list to fill:
  result_list <- list()
  # params sds, y's, mu's, and param means table
  param_sds_list <- list()
  y_list <- list()
  mu_list <- list()
  results <- matrix(NA, nrow = length(model_list), ncol = 14)
  colnames(results) <- c("Run", "Model", 
                         "alpha_1", "alpha_2", "alpha_3", "alpha_4", "alpha_5", "alpha_6",
                         "beta_int", "beta_1", "beta_2", "beta_3", "q", "tau")
  for(i in 1:length(model_list)){
    model <- mort_data_load(dir, model_list[i])
    mod_name <- str_extract(model_list[i], "(?<=modelrun_).*(?=_data)")
    jags_out <- model$jags_out
    vars <- varnames(jags_out)
    # split up outputs into model params and modeled y's
    # model params:
    params <- jags_out[,grep("^alpha|b|q|tau", vars)]
    param_out <- as.matrix(params)
    param_means <- apply(param_out, 2, mean, na.rm = TRUE)
    param_sds <- apply(param_out, 2, sd, na.rm = TRUE)
    param_sds_list[[i]] <- param_sds
    # model y's and mu's:
    mu <- jags_out[,grep("mu", vars)]
    y <- jags_out[,grep("y", vars)]
    mu_out <- as.matrix(mu)
    y_out <- as.matrix(y)
    mu_means <- apply(mu_out, 2, mean, na.rm = TRUE)
    y_means <- apply(y_out, 2, mean, na.rm = TRUE)
    mu_sds <- apply(mu_out, 2, sd, na.rm = TRUE)
    y_sds <- apply(y_out, 2, sd, na.rm = TRUE)
    y_res <- list(y_means = y_means, 
                  y_sds = y_sds)
    mu_res <- list(mu_means = mu_means, 
                   mu_sds = mu_sds)
    y_list[[i]] <- y_res
    mu_list[[i]] <- mu_list
    
    # fill in results table:
    results[i, 1] <- i
    results[i, 2] <- mod_name
    results[i, grep("^alpha", colnames(results))] <- param_means[grep("^alpha", names(param_means))]
    results[i, grep("^beta", colnames(results))][1:2] <- param_means[grep("^b", names(param_means))]
    results[i, "q"] <- param_means["q"]
    results[i, "tau"] <- param_means["tau"]
  }
  result_list <- list(results = results, 
                      param_sds = param_sds_list, 
                      model_ys = y_list, 
                      model_mus = mu_list)
  return(result_list)
}

uni_results <- results_extract(dir, model_list)
#save(uni_results, file = "2025_05_06_univar_model_results_list.RData")


### Visualizations and parsing outputs:
# selecting model:
dir = dir
modelnum = 12
log = TRUE
# # load model and output:
model <- mort_data(dir, modelnum, log)
output <- mort_out(dir, modelnum, log)
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
# # # making pred/obs plots:
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