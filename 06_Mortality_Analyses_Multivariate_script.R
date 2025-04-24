### This is the script for running multivariate mortality analyses
### It is for the first set of analyses, using harvard forest field plot data
### Assessing drivers of plot-level mortality, using remote sensing and field observations

#### ----- Libraries ----- ####
#install.packages("librarian")
#install.packages("rjags")
librarian::shelf(tidyverse, dplyr, rjags, coda, ggplot2, combinat) 

#### ----- Load Data (if not in environment) ----- ####
# set working directory:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)
load("Environments/2025_04_23_environment.RData")

#### ----- Preparing data for models ----- ####
## make predictor data set:
pred_m <- pred %>%
  # add the sums of dmag columns (tcg dmag, cs dmag)
  mutate(dmag_y1_y2 = dmag_y1 + dmag_y2, .after = dmag_y2) %>%
  mutate(dmag_cs_y1_y2 = dmag_cs_y1 + dmag_cs_y2, .after = dmag_cs_y2) %>%
  # select the top performers:
  select(pba_oak, tba_oak, dmag_cs_y1_y2, dmag_y1_y2, dmag_y1, cs_y1)

## make matrix of models using column numbers:
# making data frame for up to 3-variable models to start
# two variable models:
two_var_models <- as.data.frame(t(cbind(combn(c(1, 3:6), 2),
                                        combn(c(2, 3:6), 2),
                                        combn(c(3, 4:6), 2)))) %>%  # combinations of variables
  # filter all with 1, 2, 3 as 'starting' variables
  filter(V1 %in% 1:3) %>%
  # select only distinct rows
  distinct(V1, V2, .keep_all = TRUE) %>%
  # add a placeholder for 3rd variable to join together
  mutate(V3 = NA)
  
# three variable models:
three_var_models <- as.data.frame(t(cbind(combn(c(1, 3:6), 3),
                                          combn(c(2, 3:6), 3)))) %>% # combinations of variables
  # filter all with 1, 2 starting and variable 2 = 3
  filter(V1 %in% 1:2 & V2 == 3) %>%
  # select only distinct rows
  distinct(V1, V2, V3, .keep_all = TRUE)

# make one data frame with all models listed:
mult_models <- rbind(two_var_models, three_var_models)  # combine into one dataframe
rm(pred, two_var_models, three_var_models)  # cleaning up

# make model data set:
model_data <- cbind.data.frame(y = resp$pdba, 
                               pred_m,
                               hs = resp$hotspot)

#### ----- running multivariate models----- ####
# multivariate mortality model:
multmort_model <- read_file("Models/2025_04_22_multi_mort_model_with_logit.txt")

# multivariate mortality model run function:
#'@param model_data = pred/resp dataframe object
#'@param model_list = dataframe of multivariate models 
#'@param model = model character vector object from .txt file
#'@param niter = number of model iterations to run
#'@param diter = number of DIC iterations to run
#'@param run = loop number for model run for id
run_multi_mort_model <- function(model_data, model_list, model, niter, diter, run){
  # identifier:
  model_run = run
  # sort by 0-1 mortality percentages before running
  # sort data by y:
  c <- vector()
  for (i in 1:nrow(model_data)){
    if (model_data$y[i] == 0){
      c[i] <- "l"   
    } else if (model_data$y[i] > 0 & model_data$y[i] < 1) {
      c[i] <- "ld"
    } else if (model_data$y[i] == 1) {
      c[i] <- "d"
    }
  }
  # add to data frame:
  model_data$c <- c
  # sort by y:
  data_sort <- model_data[order(model_data$y),]
  
  # get the predictors:
  covs <- unlist(model_list[run,]) + 1  # add plus 1 to get past the y in data_sort
  # number of predictors for making b0 and Vb:
  ncovs <- length(which(!is.na(covs)))
  
  # inputs:
  data <- list(x = cbind(int = 1, data_sort[,covs[which(!is.na(covs))]]), # add column of 1's for intercept term
               y = data_sort$y, hot = data_sort$hs, 
               sites = nrow(data_sort), hs = length(unique(data_sort$hs)),
               b0 = as.vector(rep(0, ncovs + 1)), Vb = solve(diag(10000, ncovs + 1)),  # add +1 for the intercept term
               q0 = 1, qb = 1, 
               c = length(which(data_sort$c == "l")), 
               d = length(which(data_sort$c == "l")) + length(which(data_sort$c == "ld")), 
               e = nrow(data_sort))
  
  # run the test model:
  mort_jags <- jags.model(file = textConnection(model),
                          data = data,
                          n.chains = 3)
  jags_out <- coda.samples(model = mort_jags, 
                           variable.names = c("b", "q", "tau", "alpha", "y", "mu"),
                           n.iter = niter)
  
  # run DIC
  DIC <- dic.samples(mort_jags, n.iter = diter)
  sum <- sum(DIC$deviance, DIC$penalty)
  
  ### Make output list
  # track metadata
  metadata <- tibble::lst(model, data, run, ncovs, covars = model_list[run,])#, init)
  # model selection
  dic <- list(DIC, sum)
  # model output
  out <- as.matrix(jags_out)
  # combine output
  output <- tibble::lst(metadata, dic, jags_out, out)
  
  return(output)
}

# function for saving models:
### Function for saving model output:
##' @param out_path model outputs save file path
##' @param run_path model runs save file path
##' @param jags_model output from model run function
multi_model_save <- function(out_path, run_path, jags_model){
  # choose file path
  filepath_outputs <- out_path
  filepath_runs <- run_path
  # date
  date <- as.character(Sys.Date())
  # make file name
  filename_outputs <- paste0(filepath_outputs, 
                             date, 
                             "_multi_modelrun_", as.character(jags_model$metadata$run),
                             "_covs_", jags_model$metadata$covars[1],
                             "_", jags_model$metadata$covars[2],
                             "_", jags_model$metadata$covars[3],
                             "_output",".csv")
  filename_runs <- paste0(filepath_runs,
                          date,
                          "_modelrun_", as.character(jags_model$metadata$run),
                          "_covs_", jags_model$metadata$covars[1],
                          "_", jags_model$metadata$covars[2],
                          "_", jags_model$metadata$covars[3],
                          "_data",".RData")
  
  # save outputs to folder
  write.csv(jags_model$out, file = filename_outputs)
  # save model selection and metadata to folder
  model_info <- jags_model[c('jags_out', 'dic', 'metadata')]
  save(model_info, file = filename_runs)
}

### Loop for running multivariate models
# ran 4/24/2025
for (i in 1:nrow(mult_models)){
  # additional inputs to model run functions:
  model <- multmort_model
  niter = 150000
  diter = 35000
  run = i
  # model save function inputs:
  out_path = "Mortality_Model_Runs/model_outputs/"
  run_path = "Mortality_Model_Runs/model_runs/"
  
  # run the model:
  mort_model <- run_multi_mort_model(model_data, mult_models, model, niter, diter, run)
  multi_model_save(out_path, run_path, mort_model)
}

##### ARCHIVE #####
# testing the loop
# test <- run_multi_mort_model(model_data = model_data,
#                              model_list = mult_models,
#                              model = multmort_model,
#                              niter = 1000,
#                              diter = 1000,
#                              run = 13)

#vars <- varnames(jags_out)
#params <- jags_out[,grep("^alpha|b|q|tau", vars)]
