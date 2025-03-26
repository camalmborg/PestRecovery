# Model test script for censored regression with logit link and random effects

### load libraries
library(dplyr, rjags)

#### Simulate Data
# x and y variables:
y = c(rep(0,16), runif(24, min = 0, max = 1), rep(1,10))
x1 = rnorm(length(y), mean = 0, sd = 2)
x2 = rnorm(length(y), mean = 0, sd = 3)
# make dataframe:
df <- data.frame(y = y, x = x1)

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


