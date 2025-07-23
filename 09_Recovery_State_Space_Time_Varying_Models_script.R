# recovery rate state space model script for non-static (time-varying) recovery models
# this script contains the code for running a jags model for estimating recovery rates
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)

# setting task id for running array job on SCC:
task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

# load environment if needed:
load("Environments/2025_07_07_environment.RData")