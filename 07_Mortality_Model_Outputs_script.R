# model output script

out <- as.matrix(jags_out)
pairs(out)
cor(out)
gelman.diag(jags_out)
gelman.plot(jags_out)
effectiveSize(jags_out)

# make y data for comparison plots
ymeans <- apply(out[,grep("y", colnames(out))], 2, mean)
plot(ymeans, resp$pdba)