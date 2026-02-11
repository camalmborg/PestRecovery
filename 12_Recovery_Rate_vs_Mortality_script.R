### Script for Recovery Rate vs Mortality

## Load libraries:
librarian::shelf(dplyr, tidyverse, ggplot2, ggpubr, ggpmisc)

## Load data
# set working directory:
dir <- "/projectnb/dietzelab/malmborg/Ch2_PestRecovery/"
setwd(dir)
load("Environments/2025_04_07_environment.RData")

## Get recovery rates for Harvard Forest sites between 2017 and 2024
# select start and end years:
start <- grep(as.character(2017), colnames(hf_tcg))
end <- grep(as.character(2024), colnames(hf_tcg))
# number of sites:
nsite <- length(hf_tcg[,1])
# containing slope and recovery rates:
slope <- list()
recov_rate <- vector()

for (i in 1:nsite){
  recov <- as.vector(as.matrix(hf_tcg[i, start:end]))
  recov_length <- 1:length(recov)
  slope[[i]] <- lm(recov ~ recov_length)
  recov_rate[i] <- slope[[i]]$coefficients["recov_length"]
}
rm(slope, recov, recov_length)

# add recovery slopes to hf_tcg:
hf_tcg$recovery_rate <- recov_rate

## Recovery Rate vs. Mortality:
# select plots without timber harvest:
hf_tcg_pdba <- hf_tcg[hf_tcg$plot %in% resp$plot,]
# add column for dead bin:
hf_tcg_pdba$dead_bin <- resp$dead_bin
# add plot for percent mortality:
hf_tcg_pdba$pdead <- resp$pdead
# add plot for percent dead basal area:
hf_tcg_pdba$pdba <- resp$pdba*100

# Make a nice ggplot of it:
recov_mort_plot <- ggplot(hf_tcg_pdba, aes(x = pdba, y = recovery_rate)) +
  geom_point() +
  geom_smooth(method=lm , color="red", se=FALSE) +
  stat_fit_glance(method = "lm", method.args = list(formula = y ~ x), 
                  aes(label = paste("P-value =", signif(after_stat(p.value), 3))), 
                  label.x = 0.785, label.y = 0.90, size = 6) +
  stat_regline_equation(label.x = 49, label.y = 0.02, size = 6) +
  labs(x = "Percent Dead Basal Area in Plot",
       y = "Recovery Rate (TCG per Year)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
recov_mort_plot


## Box Plot version:
# get plot-level mortality:
tree_to_plot <- trees %>%
  select(hotspot, point, plot, tot_tree, tot_dead, pdead, tot_dba, pdba) %>%
  group_by(plot) %>% summarise(across(everything(), first))
# remove the one with no trees:
notree <- setdiff(plots$plot, tree_to_plot$plot)
plot_out <- which(plots$plot == notree)
# make data:
hf_tcg_box <- hf_tcg[-plot_out,]
hf_tcg_box$tot_dead <- tree_to_plot$tot_dead
hf_tcg_box$dead_bin <- ifelse(hf_tcg_box$tot_dead > 0, "D", "L")
# timber harvest
hf_tcg_box$timber_harvest <- plots$recent_timber_harvest[-plot_out]
hf_tcg_box$timber_harvest <- ifelse(hf_tcg_box$timber_harvest == 1, "T", "S")
# make three groups:
hf_tcg_box$box_groups <- paste0(hf_tcg_box$dead_bin, "-", hf_tcg_box$timber_harvest)
hf_tcg_box$box_groups <- ifelse(grepl("T", hf_tcg_box$box_groups), "Timber Harvest",
                                ifelse(grepl("D", hf_tcg_box$box_groups), "Dead", "Live"))

# ANOVA:
anova_box <- aov(hf_tcg_box$recovery_rate ~ hf_tcg_box$box_groups, data = hf_tcg_box)
p_value <- round(summary(anova_box)[[1]][["Pr(>F)"]][1], 9)

# Boxplot
hf_box <- ggplot(hf_tcg_box, aes(x = box_groups, y = recovery_rate, fill = box_groups)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, size = 0.75) +
  scale_fill_manual(values = c("Timber Harvest" = "#F0E442", "Dead" = "#0072B2", "Live" = "#D55E00")) +
  labs(title = "Mortality Outcomes vs Recovery Rates",
       x = "Group",
       y = "Recovery Rates") +
  #stat_compare_means(method = "anova", size = 6, label.y = -0.0075, label.x = 2.5) +
  annotate("text", x = 2.98, y = -0.009, label = paste0("p-value: ", p_value), size = 5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
hf_box

# Tukey's Honest Significant Difference test
tukey <- TukeyHSD(anova_box, conf.level = 0.95)
tukey <- as.data.frame(tukey$`hf_tcg_box$box_groups`)

# make a table for supplement:
tukey_table <- tukey |>
  # add group name to columns:
  mutate(Group = rownames(tukey), .before = 1) |>
  # rounding:
  mutate(across(where(is.numeric), function(x) round(x, 4))) |>
  # renaming columns:
  rename("Lower Hinge" = lwr) |>
  rename("Upper Hinge" = upr) |>
  rename("Diff" = diff) |>
  rename("p-value" = "p adj") |>
  # neater p-values:
  mutate("p-value" = ifelse(`p-value` < 0.001, "< 0.001", round(`p-value`, 3)))


