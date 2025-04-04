### This is the script for running mortality analyses
### It is for the first set of analyses, using harvard forest field plot data
### Assessing drivers of plot-level mortality, using remote sensing and field observations

#### ----- Libraries ----- ####
#install.packages("librarian")
#install.packages("rjags")
librarian::shelf(tidyverse, dplyr, rjags, ggplot2)  # removed: mgcv, AER, nlme, VGAM, lme4, remotes, plm, censReg

#### ----- Load Data (if not in environment) ----- ####
# disturbance magnitude data
# harvard forest plot data

#### ----- prepping data to run models ----- ####
### preparing plot-level information from tree data:
tree_to_plot <- trees %>%
  select(hotspot, point, plot, tot_tree, tot_dead, pdead, tot_dba, pdba) %>%
  group_by(plot) %>% summarise(across(everything(), first)) 

notree <- setdiff(plots$plot, tree_to_plot$plot)
plot_out <- which(plots$plot == notree)
tree_to_plot$harv = plots$recent_timber_harvest[-plot_out]
tree_to_plot$invas = plots$invasives[-plot_out]

### making the matrix for predictors variables:
## disturbance magnitude and tcg/scores during disturbance:
cols <- c(grep("^20", names(hf_dmag)), 
          grep("^dmag", names(hf_dmag)))
dmags <- cbind(hf_dmag[cols], hf_dmag_cs[cols])[-which(plots$plot == notree),]  # remote no-tree plot
colnames(dmags) <- c("tcg_y1", "tcg_y2", "dmag_y1", "dmag_y2",
                     "cs_y1", "cs_y2", "dmag_cs_y1", "dmag_cs_y2")

## oak dead basal area:
oak_to_plot <- oaktrees %>%
  select(hotspot, point, plot, tot_oak, tba_oak, pba) %>%
  group_by(plot) %>% summarise(across(everything(), first))
oak_to_plot <- oak_to_plot[-which(plots$plot == notree),]

oak_dat <- tree_to_plot %>%
  select(plot, hotspot, point) %>%
  left_join(., oak_to_plot, by = 'plot', na_matches = "na") %>%
  replace(is.na(.), 0)

## make predictor variable data set:
pred <- cbind(dmags,                           # disturbance magnitude and disturbance year tcg/scores data
              oak_dat$tba_oak, oak_dat$pba)#,    # oak tree density data
              #dayms)                           # daymet data

## make response variable data set:
resp <- cbind.data.frame(plot = tree_to_plot$plot, 
                         hotspot = tree_to_plot$hotspot,
                         dead = tree_to_plot$tot_dead > 0,
                         pdead = (tree_to_plot$pdead/100),      # make a percentage (btw 0-1 instead of 0-100)
                         pdba = (tree_to_plot$pdba/100)) %>%    # make a percentage 
  mutate(dead_bin = ifelse(dead == TRUE, 1, 0), .after = 2)

## remove salvage logged sites:
resp <- resp[tree_to_plot$harv == 0,]
pred <- pred[tree_to_plot$harv == 0,]


#### ----- running models----- ####
# sort by 0-1 mortality percentages before running

# model with logit link:
model_log <- "model{
### Loop over individual sites

	### Data Model:
	## left (0) censored:
	for (i in  1:c){
		y[i] ~ dbern(theta[i])
		theta[i] <- pnorm(0, mu[i], tau)  
	}
	
		## between 0-1:
	for (i in (c+1):d){
		y[i] ~ dnorm(mu[i], tau)
	}

	## right (1) censored:
	for (i in (d+1):e){
		y[i] ~ dbern(theta[i])
		theta[i] <- 1 - pnorm(1, mu[i], tau)
	}

	### Process Model:
	for (s in 1:sites) {
	logit(mu[s]) <- b[1] + b[2]*x[s] + alpha[hot[s]]
	}

  ### Random effect for hotspot:
  for (h in 1:hs) {
  	alpha[h] ~ dnorm(0, q)
	}

### Priors:
b ~ dmnorm(b0, Vb)
q ~ dgamma(q0, qb)
tau ~ dgamma(0.001, 1)
}"

# model without logit link:
model_nolog <- "model{
### Loop over individual sites

	### Data Model:
	## left (0) censored:
	for (i in  1:c){
		y[i] ~ dbern(theta[i])
		theta[i] <- pnorm(0, mu[i], tau)  
	}
	
		## between 0-1:
	for (i in (c+1):d){
		y[i] ~ dnorm(mu[i], tau)
	}

	## right (1) censored:
	for (i in (d+1):e){
		y[i] ~ dbern(theta[i])
		theta[i] <- 1 - pnorm(1, mu[i], tau)
	}

	### Process Model:
	for (s in 1:sites) {
	mu[s] <- b[1] + b[2]*x[s] + alpha[hot[s]]
	}

  ### Random effect for hotspot:
  for (h in 1:hs) {
  	alpha[h] ~ dnorm(0, q)
	}

### Priors:
b ~ dmnorm(b0, Vb)
q ~ dgamma(q0, qb)
tau ~ dgamma(0.001, 1)
}"

# make test data:
test_data <- cbind.data.frame(y = resp$pdba, x1 = pred$tcg_y1, x2 = pred$tcg_y2, hs = resp$hotspot)

# sorting criteria:
c <- vector()
for (i in 1:nrow(test_data)){
  if (test_data$y[i] == 0){
    c[i] <- "l"   
  } else if (test_data$y[i] > 0 & test_data$y[i] < 1) {
    c[i] <- "ld"
  } else if (test_data$y[i] == 1) {
    c[i] <- "d"
  }
}
# add to data frame:
test_data$c <- c
# sort by y:
test_data_sort <- test_data[order(test_data$y),]

# inputs:
data <- list(x = test_data_sort$x1, y = test_data_sort$y, hot = test_data_sort$hs, 
             sites = nrow(test_data_sort), hs = length(unique(test_data_sort$hs)),
             b0 = as.vector(c(0,0)), Vb = solve(diag(10000, 2)), 
             q0 = 1, qb = 1, 
             c = length(which(test_data_sort$c == "l")), 
             d = length(which(test_data_sort$c == "l")) + length(which(test_data_sort$c == "ld")), 
             e = nrow(test_data_sort))

# run the test model:
jags_test <- jags.model(file = textConnection(model_log),
                        data = data,
                        n.chains = 3)
jags_out <- coda.samples(model = jags_test, 
                         variable.names = c("b", "q", "tau", "alpha", "y"),
                         n.iter = 10000)
#plot(jags_out)

out <- as.matrix(jags_out)
pairs(out)
cor(out)
gelman.diag(jags_out)
gelman.plot(jags_out)
effectiveSize(jags_out)

# make y data for comparison plots
ymeans <- apply(out[,grep("y", colnames(out))], 2, mean)
plot(ymeans, resp$pdba)

# model DICs:
dic <- dic.samples(jags_test, n.iter = 20000)

#### Archive ####-----------------------------------------------------------------------####

#### Daymet Code: 
# ## daymet data:
# # making seasonal values:
# spring = c(3,4,5)
# summer = c(6,7,8)
# # spring seasonal aggregates:
# spr_daym <- hf_daym %>%
#   rename(spr_precip = prcp..mm.day., spr_tmax = tmax..deg.c., spr_tmin = tmin..deg.c., spr_vpd = vp..Pa.) %>%
#   group_by(site, year) %>% filter(month %in% spring) %>%
#   summarise(across(c(spr_precip, spr_tmax, spr_tmin, spr_vpd), mean))
# # summer seasonal aggregates:
# summ_daym <- hf_daym %>%
#   rename(summ_precip = prcp..mm.day., summ_tmax = tmax..deg.c., summ_tmin = tmin..deg.c., summ_vpd = vp..Pa.) %>%
#   group_by(site, year) %>% filter(month %in% summer) %>%
#   summarise(across(c(summ_precip, summ_tmax, summ_tmin, summ_vpd), mean))
# # spring/summer monthly means for 2016:
# daym16 <- hf_daym %>%
#   filter(., month %in% c(spring, summer) & year == 2016) %>%
#   rename(precip = prcp..mm.day., tmax = tmax..deg.c., tmin = tmin..deg.c., vpd = vp..Pa.) %>%
#   mutate(date = month.abb[month]) %>% 
#   pivot_wider(names_from = date, 
#               values_from = c(precip, tmax, tmin, vpd)) %>%
#   group_by(site) %>% summarise(across(everything(), sum, na.rm = T)) %>%
#   ungroup()%>%
#   rename_at(vars(-c(site, year, month)), ~paste0("2016_", .))
#   
# # combine data for seasonal 2014-2015:  
# daym <- cbind(spr_daym[which(spr_daym$year == 2014), grep("^spr", names(spr_daym))], 
#               summ_daym[which(summ_daym$year == 2014), grep("^summ", names(summ_daym))]) %>%
#                 rename_with(~paste0("2014_", .x)) %>%
#   bind_cols(spr_daym[which(spr_daym$year == 2015), grep("^spr", names(spr_daym))],
#             summ_daym[which(summ_daym$year == 2015), grep("^summ", names(summ_daym))]) %>%
#   rename_at(vars(starts_with('s')), ~paste0('2015_', .)) %>%
#   bind_cols(spr_daym[which(spr_daym$year == 2016), grep("^spr", names(spr_daym))],
#             summ_daym[which(summ_daym$year == 2016), grep("summ", names(summ_daym))]) %>%
#   rename_at(vars(starts_with('s')), ~paste0('2016_', .))
# 
# # combine daymet data:
# dayms <- cbind(daym, daym16[,grep("^2", names(daym16))])[-which(plots$plot == notree),]
# # clean up:
# rm(daym, daym16, spr_daym, summ_daym)




#### Model test code:
# run a test first
# test_data <- cbind.data.frame(y = resp$pdead, x1 = pred$tcg_y1, x2 = pred$tcg_y2, hs = resp$hotspot)
# test_data <- pdata.frame(test_data, index = "hs")

# # using jags:
# test_model <- "
# ## Jags code for mortality analyses
# model{
# 
# ### Loop over individual sites
# for (s in 1:sites) {
# 	
# 	### Data Model:
# 	y[s] ~ dnorm(mu[s], p)
# 
# 	### Process Model:
# 	logit(mu[s]) <- b[1] + b[2]*x[s] + alpha[hot[s]]
# }
# 
# 
# ### Random effect for hotspot:
# for (h in 1:hs) {
# 	alpha[h] ~ dnorm(0, tau)
# 	}
# 
# ### Priors:
# b ~ dmnorm(b0, Vb)
# p ~ dgamma(p0, pb)
# tau ~ dgamma(0.001, 1)
# }"
# 
# # data for the model:
# data <- list(x = pred$tcg_y1, y = resp$pdead, 
#              hot = resp$hotspot, sites = nrow(resp),
#              hs = length(unique(resp$hotspot)),
#              b0 = as.vector(c(0,0)), Vb = solve(diag(10000, 2)),
#              p0 = 0.1, pb = 0.1)
# 
# # run the test model:
# jags_test <- jags.model(file = textConnection(test_model),
#                         data = data,
#                         n.chains = 3)
# jags_out <- coda.samples(model = jags_test, 
#                          variable.names = c("b", "p", "alpha"),
#                          n.iter = 50000)
# plot(jags_out)
# 
# out <- as.matrix(jags_out)
# pairs(out)
# cor(out)
# gelman.diag(jags_out)
# gelman.plot(jags_out)
# effectiveSize(jags_out)



#install.packages("censReg")
#install.packages("AER")
#install.packages("remotes")
# using beta regression in mgcv package
# function example: gam(y ~ s(x), family = betar(link = "logit"), data = data)
#test <- gam(y ~ s(x), family = betar(link = "logit"), data = test_data)

# library(VGAM)
# model <- vgam(y_censored ~ x + vgam(1 | group), family = tobit(), data = data.frame(y_censored, x, group))
# 
# # Print the model summary
# summary(model)

# vglm(formula = apt ~ read + math + prog, family = tobit(Upper = 800), 
#      ##     data = dat)
# test_model <- vglm(y ~ x, family = tobit(Lower = 0, Upper = 1), data = test_data)

# # Install packages if not already installed
# # install.packages(c("AER", "lme4"))  # Or "nlme" depending on preference
# library(AER)
# library(lme4) # Or library(nlme)
# # Example with lme4
# # Assuming your data frame is called 'my_data'
# # 'y' is the dependent variable, 'x1', 'x2' are independent variables, and 'group' is the grouping variable
# # 'left' is the censoring threshold (if applicable)


# # Fit the Tobit model with a logit link and random effects
#model <- lmer(y ~ x1 + x2 + (1|group), data = my_data, family = binomial(link = "logit"))
# # Or with nlme
#model <- lme(y ~ x1 + x2, random = ~1|group, data = my_data, family = binomial(link = "logit"))

#test_model <- lme(y ~ x, random = ~1|hs, data = test_data, family = binomial(link = "logit"))
#test_model <- lmer(y ~ x + (1|hs), data = test_data, family = binomial(link = "logit"))
# 
# # Print the model summary
# #summary(model)

# try with VGAM:
# test_model <- vglm(y ~ x1, 
#                    family = tobit(Lower = 0, Upper = 1, lmu = "logitlink", type.fitted = c("censored")),
#                    data = test_data)   # works but does not include random effects
# fit <- fitted.values(test_model)

# test_model <- tobit(y ~ x1 + (1|hs),    # works but does not work with random effect added
#                     left = 0, right = 1,
#                     dist = "gaussian",
#                     data = test_data)
# 
# test_model_1 <- censReg(y ~ x1, left = 0, right = 1, method = "BHHH", data = test_data)   # works but does not easily have fitted values estimation
# summary(test_model_1)
# 
# test_model_2 <- censReg(y ~ x1 + x1* x2, left = 0, right = 1, method = "BHHH", data = test_data)
# summary(test_model)
# AIC(test_model)