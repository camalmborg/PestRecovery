### Identifying significant variables with confidence interval overlap
# Script for working with slope (beta) confidence intervals

## Load libraries
library(dplyr)
library(tidyverse)
library(rjags)

## set working directory
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Recovery_State_Space_Runs/"
setwd(dir)

## Load model output files
model_params <- read.csv("2025_07_31_all_base_uni_recov_models_param_means.csv")
load("2025_07_31_recov_models_outputs_list.RData")  # object is called model_outputs


## Calculating CIs for all models
# getting means, low CI (0.025 quantile), high CI (0.975 quantile):
model_means <- matrix(NA, nrow = length(model_outputs), ncol = 6)
model_lows <- matrix(NA, nrow = length(model_outputs), ncol = 6)
model_highs <- matrix(NA, nrow = length(model_outputs), ncol = 6)
for (i in 1:length(model_outputs)){
  # add model name:
  model_means[i,1] <- names(model_outputs)[i]
  model_lows[i,1] <- names(model_outputs)[i]
  model_highs[i,1] <- names(model_outputs)[i]
  # get results from list:
  result <- model_outputs[[i]]
  if (TRUE %in% grepl("beta", colnames(result))){
    # separate beta cols
    beta_means <- apply(result[grep("beta", colnames(result))], 2, mean)
    beta_lows <- apply(result[grep("beta", colnames(result))], 2, quantile, c(0.025))
    beta_highs <- apply(result[grep("beta", colnames(result))], 2, quantile, c(0.975))
    for (j in 1:length(beta_means)){
      model_means[i,j+1] <- round(beta_means[j], 3)
      model_lows[i,j+1] <- round(beta_lows[j], 3)
      model_highs[i,j+1] <- round(beta_highs[j], 3)
    }
  }
}
# remove unnecessary things:
rm(result)

# making one data frame with all the results compiled:
model_results_CIs <- as.data.frame(model_means) %>%
  rename(model = 1) %>%
  rename_with(~ paste0("mean_", seq_along(.)), .cols = -c(model)) %>%
  # add the lower CI
  bind_cols(model_lows[,2:ncol(model_lows)]) %>%
  rename_with(~ paste0("low_", seq_along(.)), .cols = -c(1:6)) %>%
  # add the higher CI
  bind_cols(model_highs[,2:ncol(model_highs)]) %>%
  rename_with(~ paste0("high_", seq_along(.)), .cols = -c(1:11))

# fix the categorcial:
model_results_CIs_cat <- model_results_CIs %>%
  filter(.[[1]] == "uni_nlcd_cat") %>%
  pivot_longer(
    cols = matches("^(mean|low|high)_\\d+"),
    names_to = c(".value", "set"),
    names_pattern = "(mean|low|high)_(\\d+)") %>%
  mutate(model = paste0("uni_nlcd_cat_", set), .before = 2) %>%
  select(-set)

# bind to original, fix column names:
model_results_CIs <- model_results_CIs %>%
  filter(!.[[1]] == "uni_nlcd_cat") %>%
  mutate(mean = mean_1) %>%
  mutate(low = low_1) %>%
  mutate(high = high_1) %>%
  select(model, mean, low, high) %>%
  bind_rows(model_results_CIs_cat)
#remove unnecessary things:
rm(model_results_CIs_cat)


## Ridge Plot for betas
# load libraries:
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)

# get the data together:
beta_list <- list()
for (i in seq_along(model_outputs)) {
  model_name <- names(model_outputs)[i]
  result <- model_outputs[[i]]
  
  # Find all columns with 'beta' in the name
  beta_cols <- grep("beta", colnames(result), value = TRUE)
  if (length(beta_cols) > 0) {
    beta <- result[, beta_cols, drop = FALSE]  # keep as data frame
    # Rename columns: beta_1_modelName, beta_2_modelName, etc.
    new_names <- paste0("beta_", seq_along(beta_cols), "_", model_name)
    colnames(beta) <- new_names
    beta_list[[i]] <- beta
  } else {
    # No beta columns → create NA columns with proper row count
    n_rows <- nrow(result)
    na_col <- data.frame(matrix(NA, nrow = n_rows, ncol = 1))
    colnames(na_col) <- model_name
    beta_list[[i]] <- na_col
  }
}
# unlist and combine:
beta_ridges <- do.call(cbind, beta_list)
rm(na_col, beta, n_rows, result)
# rename columns for models:
#colnames(beta_ridges)

# pivot longer to prepare for making ridges:
beta_ridges_long <- beta_ridges %>%
  select(contains("beta")) %>%
  pivot_longer(., 
               cols = everything(),
               names_to = "model",
               values_to = "beta_est") %>%
  arrange(model) %>%
  # to remove categorical...
  filter(!str_detect(model, "nlcd")) %>%
  # make the names better for plot:
  mutate(model = gsub("beta_1_uni_vp", "beta_1_uni_vpd", model)) %>%
  mutate(model = gsub("_", " ", model)) %>%
  mutate(model = gsub("beta 1 uni ", "", model)) %>%
  mutate(model = gsub("mean", "2014-15 mean", model)) %>%
  mutate(model = gsub("vpdd", "vpd", model)) %>%
  mutate(model = gsub("prcp", "precip", model)) %>%
  group_by(model) %>%
  # arrange by descending mean to make ridges look nice
  mutate(mean_beta = mean(beta_est)) %>% ungroup() %>%
  mutate(model = fct_reorder(model, mean_beta, .desc = TRUE)) %>%
  # make model a factor for plotting:
  mutate(model = as.factor(model)) %>%
  arrange(model) 

# Selecting different groups:
betas_dist_hist <- beta_ridges_long %>%
  filter(grepl("tcg", model))
betas_pre_dist <- beta_ridges_long %>%
  filter(grepl("20", model))
betas_post_dist <- beta_ridges_long %>%
  filter(grepl("tcg|20", model) == FALSE)

# Now let's make the ridge plot:
beta_ridge_plot <- ggplot(betas_post_dist, aes(x = beta_est, y = model, fill = after_stat(x))) +
  # title:
  labs(title = 'Beta Parameter (Slope) Estimates') +
  # axes titles:
  xlab("Estimated Value") +
  ylab("Model Covariates") +
  # add color scaling:
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.001, 
                               bandwidth = 0.001, lwd = 0.5) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0) +
  # add vertical line at 0:
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  # set limit x-axis:
  xlim(c(-0.1, 0.13)) + # for all + post-dist
  #xlim(c(-0.055, 0.05)) + # for pre-dist
  #xlim(c(-0.09, 0.001)) + # for dist hist
  theme_bw() +
  theme(legend.position = "none")

beta_ridge_plot

# save them:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/"
setwd(save_dir)
# all vars:
# Save the plot to a PNG file
ggsave("2025_08_06_ridges_post_dist_ESA.png", 
       plot = beta_ridge_plot,
       dpi = 600)
