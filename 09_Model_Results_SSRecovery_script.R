# Checkin' out the recovery state space models

## Set up
# load libraries:
librarian::shelf(rjags, coda, dplyr)

# set correct working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

# model files:
models <- list.files(paste0(dir, "model_runs"))




# load model
load(paste0(dir, "model_runs/", models[3]))

# load jags output:
jags_out <- model_info$jags_out
vars <- varnames(jags_out)
params <- jags_out[,grep("r0|^tau", vars)]
R_samp <- sample(vars[grep("R", vars)], 6)
x_samp <- sample(vars[grep("x", vars)], 6)
R_params <- jags_out[,R_samp]
x_params <- jags_out[,x_samp]

# for the random effects model
atime_params <- jags_out[,grep("^at", vars)]
asite_samp <- sample(vars[grep("^as", vars)], 6)
asite_params <- jags_out[,asite_samp]


# collecting DICs
model_dics <- matrix(NA, nrow = length(models), ncol = 3)
for (i in 1:length(models)){
  # load model
  load(paste0(dir, "model_runs/", models[i]))
  # extract DIC
  
}