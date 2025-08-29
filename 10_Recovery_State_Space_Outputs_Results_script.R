# Checkin' out the recovery state space models

## Set up
# load libraries:
librarian::shelf(rjags, coda, dplyr, stringr)

# set correct working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

# model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]
# for multivariate:
#models <- models[which(grepl("multi", models) == T)]

### Collecting DICs
# extract from model metadata:
model_dics <- matrix(NA, nrow = length(models), ncol = 3)
colnames(model_dics) <- c("model_number", "covariate", "dic")
for (i in 1:length(models)){
  # load model:
  load(paste0(dir, "model_runs/", models[i]))
  # model name:
  name <- models[i]
  if (grepl("cov", name) == TRUE){
    model <- str_match(name, "cov_(.*?)_data")[,2]
  } else {
    model <- str_match(name, "_(.*?)_model")[,2]
  }
  model_dics[i,1] <- i
  model_dics[i,2] <- model
  # extract DIC:
  dic <- model_info$dic[[2]]
  # fill in table:
  model_dics[i,3] <- dic
  print(i)
}

# get deltaDICs:
dic_sort <- as.data.frame(model_dics[order(as.numeric(model_dics[,"dic"])),], decreasing = TRUE)
dic_sort$del_dic <- as.numeric(dic_sort$dic) - min(as.numeric(dic_sort$dic)) 
# save:
write.csv(dic_sort, "2025_08_29_all_recov_models_dics.csv")


### Getting beta parameters and calculating CIs
# # testing...
#models <- list.files(paste0(dir, "model_runs", "/archive"))#[grep("RData", list.files(paste0(dir, "model_runs")))]
#load(paste0(dir, "model_runs/archive/", models[8]))

#varnames(params)
model_params <- as.data.frame(matrix(NA, nrow = length(models), ncol = 12))
model_outputs <- list()
colnames(model_params) <- c("model_number", "covariate", 
                            "atime[1]", "atime[2]", "atime[3]", "atime[4]", "atime[5]", "atime[6]",
                            "beta", "r0", "tau_add", "tau_obs")
for (i in 1:length(models)){
  # progress bar:
  print(i)
  # load model:
  load(paste0(dir, "model_runs/", models[i]))
  # model name:
  name <- models[i]
  if (grepl("cov", name) == TRUE){
    model <- str_match(name, "cov_(.*?)_data")[,2]
  } else {
    model <- str_match(name, "_(.*?)_model")[,2]
  }
  # load model number and name:
  model_params[i,1] <- i
  model_params[i,2] <- model
  
  # collect model params:
  jags_out <- model_info$jags_out
  vars <- varnames(jags_out)
  params <- jags_out[,grep("r0|^tau|^b|^at", vars)]
  # remove burn in:
  burn_in = 25000
  params_burn <- window(params, start = burn_in)
  # convert output to dataframe:
  params_out <- as.data.frame(as.matrix(params_burn))
  # save model params to list:
  model_outputs[[i]] <- params_out
  names(model_outputs)[i] <- model
  # get param means:
  params_mean <- apply(params_out, 2, mean)
  # fill in table of parameter means:
  for (j in varnames(params)){
    model_params[i,j] <- round(params_mean[j], 3)
  }
}

# save:
write.csv(model_params, "2025_08_29_all_recov_models_param_means.csv")
save(model_outputs, file = "2025_08_29_all_recov_models_outputs_list.RData")

# remove things I don't need:
rm(model_info, params, params_burn, params_out, jags_out, vars)


# ### For model convergence checks and tests:
# # # load model
# load(paste0(dir, "model_runs/", models[1]))
# #
# # # load jags output:
# jags_out <- model_info$jags_out
# vars <- varnames(jags_out)
# params <- jags_out[,grep("r0|^tau", vars)]
# R_samp <- sample(vars[grep("R", vars)], 6)
# x_samp <- sample(vars[grep("x", vars)], 6)
# R_params <- jags_out[,R_samp]
# x_params <- jags_out[,x_samp]
# 
# # for the random effects model
# atime_params <- jags_out[,grep("^at", vars)]
# # asite_samp <- sample(vars[grep("^as", vars)], 6)
# # asite_params <- jags_out[,asite_samp]
# 
# # for betas (univariate):
# beta_params <- jags_out[,grep("^b", vars)]

# testing convergence with plots
# testing effective sizes
# gelman to determine burn in
