# Model test script for censored regression with logit link and random effects

#### load libraries
library(dplyr)
library(rjags)

#### Simulate Data
# x and y variables:
y = c(rep(0,16), runif(24, min = 0, max = 1), rep(1,10))
x1 = rnorm(length(y), mean = 0, sd = 2)
x2 = rnorm(length(y), mean = 0, sd = 3)
# group simulation:
hs = c(rep(1,5), rep(2,12), rep(3,8), rep(4,9), rep(5,10), rep(6,6))
# make dataframe:
df <- data.frame(y = y, x = x1, hs = hs)

# sorting criteria:
c <- vector()
for (i in 1:nrow(df)){
  if (df$y[i] == 0){
    c[i] <- "l"   
  } else if (df$y[i] > 0 & df$y[i] < 1) {
    c[i] <- "ld"
  } else if (df$y[i] == 1) {
    c[i] <- "d"
  }
}
# add to data frame:
df$c <- c
# sort by y:
df <- df[order(df$y),]

#### Running JAGS model
# model object:
model <- "model{
### Loop over individual sites

	### Data Model:
	## left (0) censored:
	for (i in  1:c){
		y[i] ~ dbern(theta[i])
		theta[i] <- pnorm(0, mu[i], tau)  
	}
	
		## between 0-1:
	for (i in (c+1):d){
		y[i] ~ dnorm(mu[i], p)
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
p ~ dgamma(p0, pb)
q ~ dgamma(q0, qb)
tau ~ dgamma(0.001, 1)
}"

# inputs:
data <- list(x = df$x, y = df$y, hot = df$hs, sites = nrow(df), hs = length(unique(df$hs)),
             b0 = as.vector(c(0,0)), Vb = solve(diag(10000, 2)), 
             p0 = 0.01, pb = 0.001,  
             q0 = 0.01, qb = 0.001, 
             c = length(which(df$c == "l")), 
             d = length(which(df$c == "l")) + length(which(df$c == "ld")), 
             e = nrow(df))

# run the model:
jags_test <- jags.model(file = textConnection(model),
                        data = data,
                        n.chains = 3)
jags_out <- coda.samples(model = jags_test, 
                         variable.names = c("b", "p", "alpha", "y"),
                         n.iter = 10000)

# pushing changes