# Script for visualizing, parsing, and plotting model outputs from multivariate mortality models

# load libraries
librarian::shelf(tidyverse, dplyr, rjags, coda)

# navigate to folder:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Mortality_Model_Runs/"
setwd(dir)

### Functions to get model outputs
# get model output:
#'@param dir = folder where the model_outputs and model_runs are located
#'@param modelnum = model number for multivariate run
mort_out <- function(dir, modelnum){
  num <- modelnum
  file <- list.files(paste0(dir, "model_outputs/"))
  multis <- file[grep("multi", file)]
  name <-  grep(paste0("modelrun_", as.character(modelnum), "_"), multis)
  model <- read.csv(paste0(dir, "model_outputs/", multis[name]))
  return(model)
}

# get model result object:
mort_data <- function(dir, modelnum){
  num <- modelnum
  file <- list.files(paste0(dir, "model_runs/"))
  multis <- file[grep("covs", file)]
  name <-  grep(paste0("modelrun_", as.character(modelnum), "_"), multis)
  load(paste0(dir, "model_runs/", multis[name]))
  return(model_info)
}

mort_data_load <- function(dir, model){
  load(paste0(dir, "model_runs/", model))
  return(model_info)
}

### loop for recording multivariate dics
# model_list:
models <- list.files(paste0(dir, "model_runs"))
model_list <- models[grep("covs", models)]
# make empty list for filling in dics:
mort_dics <- matrix(NA, nrow = length(model_list), ncol = 5)
#colnames(mort_dics) <- c("id", "model", "DIC")
colnames(mort_dics) <- c("id", "model", "ncovs" ,"covars", "DIC")
# the loop:
for (i in 1:length(model_list)){
  # selecting model:
  dir = dir
  # load model and output:
  model <- mort_data_load(dir, model_list[i])
  # collect dic and model identity:
  mort_dics[i, 1] <- i
  mort_dics[i, 2] <- model$metadata$run
  # number of covariates:
  mort_dics[i, 3] <- model$metadata$ncovs
  # covariates included (column numbers of multivariate preds object):
  mort_dics[i, 4] <- paste(model$metadata$covars, collapse = " ")
  # model dic:
  mort_dics[i, 5] <- model$dic[[2]]
  # tracking loop repeats:
  print(i)
}
mort_dic <- as.data.frame(mort_dics)
mort_dic$delDIC <- as.numeric(mort_dic$DIC) - min(as.numeric(mort_dic$DIC))
mort_dic <- mort_dic[order(mort_dic$delDIC),]
mort_dic$perform <- 1:nrow(mort_dic)
#write.csv(mort_dic, file = "2025_04_29_multi_mortality_models_dic_table.csv")

### Comparing with univariate models
# load univariate results:
uni_mort_dics <- read.csv("2025_04_23_mortality_models_dic_table.csv")[-1]
# load multivariate results:
multi_mort_dics <- read.csv("2025_04_29_multi_mortality_models_dic_table.csv")[-1]
# make a full dataframe with DICs from both uni and multi models:
uni <- cbind(model = uni_mort_dics$model, 
             type = "uni",
             mod_info = uni_mort_dics$log.T.F,
             DIC = uni_mort_dics$DIC)
multi <- cbind(model = multi_mort_dics$model,
               type = "multi",
               mod_info = multi_mort_dics$covars,
               DIC = multi_mort_dics$DIC)
# combine:
all_dic <- as.data.frame(rbind(uni, multi))
all_dic$delDIC <- as.numeric(all_dic$DIC) - min(as.numeric(all_dic$DIC))
all_dic <- all_dic[order(all_dic$delDIC),]
all_dic$perform <- 1:nrow(all_dic)
#write.csv(all_dic, file = "2024_04_29_all_mort_model_dics.csv")

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
    # extract jags results
    out <- model$jags_out
    # remove burn in:
    burnin <- 25000
    jags_out <- window(out, start = burnin)
    # get variable names:
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
    mu_list[[i]] <- mu_res
    
    # fill in results table:
    results[i, 1] <- i
    results[i, 2] <- mod_name
    results[i, grep("^alpha", colnames(results))] <- param_means[grep("^alpha", names(param_means))]
    if(length(grep("^b", names(param_means))) == 3) {
      results[i, grep("^b", colnames(results))[1:3]] <- param_means[grep("^b", names(param_means))]
      } else {
      results[i, grep("^b", colnames(results))] <- param_means[grep("^b", names(param_means))] 
      }
    results[i, "q"] <- param_means["q"]
    results[i, "tau"] <- param_means["tau"]
  }
  result_list <- list(results = results, 
                      param_sds = param_sds_list, 
                      model_ys = y_list, 
                      model_mus = mu_list)
  return(result_list)
}

multi_results <- results_extract(dir, model_list)
#save(multi_results, file = "2025_05_01_multivar_model_results_list.RData")
mort_results_table <- multi_results$results
#write.csv(mort_results_table, file = "2025_09_19_mort_results_table.csv")

### Visualizations and parsing outputs:
# selecting model:
dir = dir
modelnum = 1
# # load model and output:
model <- mort_data(dir, modelnum)
#output <- mort_out(dir, modelnum)

# visual inspections:
jags_out <- model$jags_out
vars <- varnames(jags_out)
params <- jags_out[,grep("^alpha|b|q|tau", vars)]

plot(params)
gelman.diag(params)
gelman.plot(params)
effectiveSize(params)
# visual inspect notes:
# model 3: converges very well, burn in ~25,000?, effective sizes all look great
