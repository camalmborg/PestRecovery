# Checkin' out the recovery state space models

## Set up
# load libraries:
librarian::shelf(rjags, coda, dplyr, stringr)

# set correct working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

# model files:
models <- list.files(paste0(dir, "model_runs"))[grep("RData", list.files(paste0(dir, "model_runs")))]


# # load model
# load(paste0(dir, "model_runs/", models[1]))
# 
# # load jags output:
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
# asite_samp <- sample(vars[grep("^as", vars)], 6)
# asite_params <- jags_out[,asite_samp]
# 
# # for betas (univariate):
# beta_params <- jags_out[,grep("^b", vars)]



### Collecting DICs
# extract from model metadata:
model_dics <- matrix(NA, nrow = length(models), ncol = 4)
colnames(model_dics) <- c("model_number", "covariate", "dic", "group")
for (i in 1:length(models)){
  # load model
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
  model_dics[i,4] <- "static group"
  print(i)
}

# get deltaDICs:
dic_sort <- as.data.frame(model_dics[order(as.numeric(model_dics[,"dic"])),], decreasing = TRUE)
dic_sort$del_dic <- min(as.numeric(dic_sort$dic)) - as.numeric(dic_sort$dic)
# save:
#write.csv(dic_sort, "2025_07_21_uni_static_recov_models_dics.csv")
                               