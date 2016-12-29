### stat221 problem set 4
### Xiaolin Zhuo


library(coda)
library(networkTomography)
source("xiaolinzhuo_mcmc.R")

# read data
#dat <- read.csv("1router_allcount.dat")
data(bell.labs)
A <- bell.labs$A
reorder <- qr(A)$pivot
matX <- bell.labs$X # OD flows in matrix form
matY <- bell.labs$Y # link loads in matrix form
r <- nrow(A)
c <- ncol(A)

# use data from time point 5 only for this task
Y <- matY[5, ]


# draw ten chains
dataX <- NA
dataL <- NA

for (i in 1:10) {
  print(paste("chain", i))
  
  set.seed(i)
  res <- network_mcmc(Y, A, list(a=0, b=0))
  print("acceptance rates:")
  print(res[["accepts"]])
  
  # create autocorrelation plots
  pdf(paste("xiaolinzhuo_q2_acf_", i, ".pdf", sep=""))
  par(mfrow=c(3, 3))
  for (j in (r+1):c) {
    acf(res[["X"]][, j], lag.max=10000, main=paste("X", j, sep=""))
  }
  dev.off()
  
  # effective sizes
  for (j in (r+1):c) {
    print(paste("X", j, ":", as.numeric(effectiveSize(res[["X"]][, j]))))
  }
  
  # combine chains
  if (is.na(dataX)) {
    dataX <- res[["X"]]
    dataL <- res[["L"]]
  } else {
    dataX <- rbind(dataX, res[["X"]])
    dataL <- rbind(dataL, res[["L"]])
  }
  
  print("========")
}


# plot marginal posterior histograms of Xs
pdf("xiaolinzhuo_1router_t5_uniform_1.pdf")
par(mfrow=c(4, 4))
for (j in 1:c) {
  hist(dataX[, j], xlab=paste("X", j, sep=""), main="")
  abline(v=matX[5, reorder][j], col="red") # add a vertical line for true Xs
}
dev.off()


# plot marginal posterior densities of lambdas
pdf("xiaolinzhuo_1router_t5_uniform_2.pdf")
par(mfrow=c(4, 4))
for (j in 1:c) {
  densplot(mcmc(dataL[, j]), xlab=paste("L", j, sep=""))
}
dev.off()




# another way: use the R package networkTomography
# draw ten chains
dataX2 <- NA
dataL2 <- NA

for (i in 1:10) {
  print(paste("chain", i))
  
  set.seed(i)
  #res <- network_mcmc(Y, A, list(a=0, b=0))
  res <- twMCMC(Y, A, list(a=0, b=0))
  
  for (j in (r+1):c) {
    print(paste("X", j, ":", as.numeric(effectiveSize(res[["XDraws"]][, j]))))
  }
  
  # combine chains
  if (is.na(dataX)) {
    dataX2 <- res[["XDraws"]]
    dataL2 <- res[["lambdaDraws"]]
  } else {
    dataX2 <- rbind(dataX, res[["XDraws"]])
    dataL2 <- rbind(dataL, res[["lambdaDraws"]])
  }
  
  print("========")
}


# plot marginal posterior histograms of Xs
pdf("xiaolinzhuo_1router_t5_uniform_3.pdf")
par(mfrow=c(4, 4))
for (j in 1:c) {
  hist(dataX2[, j], xlab=paste("X", j, sep=""), main="")
  abline(v=matX[5, reorder][j], col="red") # add a vertical line for true Xs
}
dev.off()


# plot marginal posterior densities of lambdas
pdf("xiaolinzhuo_1router_t5_uniform_4.pdf")
par(mfrow=c(4, 4))
for (j in 1:c) {
  densplot(mcmc(dataL2[, j]), xlab=paste("L", j, sep=""))
}
dev.off()


# create horizontal boxplots
# data contains too many rows. sample 10% of data for plotting
n <- dim(dataX2)[1]
n
tmp <- dataX2[sample(n, size=n*0.1), ]
pdf("xiaolinzhuo_1router_t5_uniform_5.pdf")
boxplot(tmp, horizontal=T, las=1, names=1:c, main="Uniform Prior")
points(x=matX[5, reorder], y=1:c, col="red", pch=16)
dev.off()
