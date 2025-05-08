# This script is for checking out model results prior to making presentation
# and publication figures

# load libraries
librarian::shelf(tidyverse, dplyr, rjags, coda, boot, statip)

# navigate to folder:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Mortality_Model_Runs/"
setwd(dir)
# environment:
#load()

## model_list:
models <- list.files(paste0(dir, "model_runs"))
uni_model_list <- models[grep("log", models)]
multi_model_list <- models[grep("covs", models)]

## loading results:
load("2025_05_06_univar_model_results_list.RData")    # univariate results
load("2025_05_01_multivar_model_results_list.RData")  # multivariate results

## loading dic results:
all_dics <- read.csv("2024_04_29_all_mort_model_dics.csv")[-1]


### Visualizations and parsing outputs:
# selecting model:
dir = dir
num = 1
model <- all_dics$model[all_dics$perform == num]
# get the model name:
if(all_dics[num,]$type == "multi"){
  model_name <- multi_model_list[grep(paste0("modelrun_", model), multi_model_list)]
} else if (all_dics[num,]$type == "uni" & all_dics[num,]$mod_info == "log") {
  model_name <- uni_model_list[grep(paste0("modelrun_", model, "_log"), uni_model_list)]
} else if (all_dics[num,]$type == "uni" & all_dics[num,]$mod_info == "nolog") {
  model_name <- uni_model_list[grep(paste0("modelrun_", model, "_nolog"), uni_model_list)]
}
# load the model:
load(paste0(dir, "model_runs/", model_name))  # will load in environment as "model_info"

# visual inspections:
jags_out <- model_info$jags_out
vars <- varnames(jags_out)
params <- jags_out[,grep("^alpha|b|q|tau", vars)]

plot(params)
gelman.diag(params)
gelman.plot(params)
effectiveSize(params)


### Preparing for making figures and other outputs
# extracting data for predicted/obserbed comparisons:
model_data <- cbind(resp, 
                    pred, 
                    dmag_y1_y2 = pred_m$dmag_y1_y2, 
                    dmag_cs_y1_y2 = pred_m$dmag_cs_y1_y2)
# sort data by response pdba:
c <- vector()
for (i in 1:nrow(model_data)){
  if (model_data$pdba[i] == 0){
    c[i] <- "l"   
  } else if (model_data$pdba[i] > 0 & model_data$pdba[i] < 1) {
    c[i] <- "ld"
  } else if (model_data$pdba[i] == 1) {
    c[i] <- "d"
  }
}
# add to data frame:
model_data$c <- c
# sort by y:
data_sort <- model_data[order(model_data$pdba),]
rm(model_data)

## extracting data from model runs to compare with observations
# opening results for each model:
model_list_num <- if (model_name %in% uni_model_list){
  which(uni_model_list == model_name)
} else {
  which(multi_model_list == model_name)
}
# get output:
jags_out <- model_info$jags_out
model_output <- as.matrix(jags_out)



predicted <- apply(model_output, 2, mean)
betas <- grep("^b", names(predicted))
taus <- grep("^t", names(predicted))
mus <- grep("^mu", names(predicted))
# data from inputs:
covars <- model_info$metadata$data$x
# getting hot spot alpha term:
alphas <- grep("^a", names(predicted))
hot <- model_info$metadata$data$hot
# # model predictions for mu term:
# mu <- vector() 
# for (i in 1:model_info$metadata$data$sites){
#   prod <- as.matrix(covars) %*% as.matrix(predicted[betas])
#   mu[i] <- inv.logit(prod[i] + predicted[alphas[hot[i]]])
# }
# tau term:
tau <- 1/predicted[taus]

c = which(data_sort$c == "l")
d = which(data_sort$c == "ld")
e = which(data_sort$c == "d")
y <- vector()
theta <- vector()
for (i in  c){
  theta[i] <- pnorm(0, mu[i], tau)  
  y[i] <- rbern(1, prob = theta[i])
}

## between 0-1:
for (i in d){
  y[i] <- rnorm(1, mu[i], tau)
}

## right (1) censored:
for (i in e){
  theta[i] <- 1 - pnorm(1, mu[i], tau)
  y[i] <- rbern(1, prob = theta[i])
}

hotspots <- as.factor(data_sort$hotspot)
plot(1:156, data_sort$pdba, pch = 20)
points(1:156, y, pch = 20, col = hotspots)


mu <- predicted[mus]
m <- vector()
for (i in 1:156){
  if(mu[i] < 0){
    m[i] <- 0
  } else if (mu[i] > 1){
    m[i] <- 1
  } else {
    m[i] <- mu[i]
  }
}

ypred <- matrix(NA, nrow = 25000, ncol = 156)
for (i in 1:25000){
  ypred[i,] <- rnorm(1, m, tau)
}
ci <- apply(ypred, 2, quantile, c(0.025, 0.975))


# # true/false if model is multivar (T) or uni var (F)
# tf_multi <- model_name %in% multi_model_list
# # get model results
# if (tf_multi == TRUE){
#   param_means <- multi_results$results[model_list_num,]
# }

#model_input <- model_info$metadata$data
