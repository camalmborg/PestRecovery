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
#model_params <- read.csv("2025_10_06_all_recov_models_param_means.csv")
load("2026_02_09_uni_recov_models_outputs_list.RData")  # object is called model_outputs


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

# fix multivariate:
model_results_CIs_multi <- model_results_CIs %>%
  filter(., grepl("multi", model)) %>%
  pivot_longer(
    cols = matches("^(mean|low|high)_\\d+"),
    names_to = c(".value", "set"),
    names_pattern = "(mean|low|high)_(\\d+)") %>%
  filter(set %in% c(1, 2, 3, 4)) %>%
  mutate(model = paste0(model, "_", set), .before = 2) %>%
  select(-set)

# bind to original, fix column names:
model_results_CIs <- model_results_CIs %>%
  filter(!.[[1]] == "uni_nlcd_cat" & !grepl("multi", model)) %>%
  mutate(mean = mean_1) %>%
  mutate(low = low_1) %>%
  mutate(high = high_1) %>%
  select(model, mean, low, high) %>%
  bind_rows(model_results_CIs_cat) %>%
  bind_rows(model_results_CIs_multi) %>%
  # remove NA rows from 2-var multivariate models:
  filter(!if_all(c(mean, low, high), is.na))
#remove unnecessary things:
rm(model_results_CIs_cat)
rm(model_results_CIs_multi)
# save:
write.csv(model_results_CIs, "2026_02_09_all_model_results_CIs.csv")


## Ridge Plot for betas
# load libraries:
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)

# remove multis (fewer rows):
full_model_outputs <- model_outputs
model_outputs <- full_model_outputs[-c(grep("multi", names(full_model_outputs)))]
  
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
#rm(na_col, beta, n_rows, result)
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
beta_ridge_plot <- ggplot(beta_ridges_long, aes(x = beta_est, y = model, fill = after_stat(x))) +
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
  #xlim(c(-0.069, 0.04)) + # for all + post-dist
  #xlim(c(-0.055, 0.05)) + # for pre-dist
  #xlim(c(-0.09, 0.001)) + # for dist hist
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

beta_ridge_plot

# save them:
save_dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/"
setwd(save_dir)
# all vars:
# Save the plot to a PNG file
png(filename = paste0(save_dir, "2026_02_09_ridges_uni_vars.png"),
    height = 8,
    width = 6,
    units = "in",
    res = 600)
beta_ridge_plot
dev.off()
# ggsave("2026_02_09_ridges_uni_vars.png",
#        plot = beta_ridge_plot,
#        height = 10,
#        width = 8,
#        units = "in",
#        dpi = 600)

## Make a quick table of results
library(grid)
library(gridExtra)
library(gt)
library(gtExtras) 

# get model performance list:
model_beta_result <- data.frame(model = unique(beta_ridges_long$model), 
                                means = unique(beta_ridges_long$mean_beta)) %>%
  mutate(means = round(means, digits = 3))
# sort by absolute value of beta mean:
model_beta_result <- model_beta_result[order(abs(model_beta_result$means), decreasing = TRUE),]
model_beta_result$perform <- 1:nrow(model_beta_result) 

# make the table:
#colnames(model_beta_result) <- c("Model Covariate", "Parameter Mean", "Model Performance")
uni_param_table <- gt(model_beta_result) |>
  cols_label(
    model = html("Model<br>Covariate"),
    means = html("Parameter<br>Mean"),
    perform = html("Model<br>Performance")) |>
  tab_style(
    style = cell_text(align = "center", weight = "bold"),
    locations = cells_column_labels(columns = everything())) |>
  cols_align(
    align = "center",
    columns = everything())
  # tab_options(table.font.size = 16) #%>%
uni_param_table

# save table
gtsave(uni_param_table,
       file = "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/Figures/2026_02_09_uni_params_table.rtf")

# print to plots tab:
# png("2026_02_09_model_beta_results_uni_vars.png")
# grid.table(model_beta_result)
# dev.off()
