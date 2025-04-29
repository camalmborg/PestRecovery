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
uni_mort_dics <- read.csv("2025_04_23_mortality_models_dic_table.csv")
# load multivariate results:
multi_mort_dics <- read.csv("2025_04_29_multi_mortality_models_dic_table.csv")
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


### Visualizations and parsing outputs:
# selecting model:
dir = dir
modelnum = 4
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
