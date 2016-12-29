### stat221 problem set 4
### Xiaolin Zhuo

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
  #res <- network_mcmc(Y, A, list(a=0.02*matX[4, reorder], b=rep(0.02, c)))
  res <- twMCMC(Y, A, list(a=0.02*matX[4, reorder], b=rep(0.02, c)))
  
  #print("acceptance rates:")
  #print(res[["accepts"]])
  
  # create autocorrelation plots
  pdf(paste("xiaolinzhuo_q3_acf_", i, ".pdf", sep=""))
  par(mfrow=c(3, 3))
  for (j in (r+1):c) {
    #acf(res[["X"]][, j], lag.max=10000, main=paste("X", j, sep=""))
    acf(res[["XDraws"]][, j], lag.max=10000, main=paste("X", j, sep=""))
  }
  dev.off()
  
  # effective sizes
  for (j in (r+1):c) {
    #print(paste("X", j, ":", as.numeric(effectiveSize(res[["X"]][, j]))))
    print(paste("X", j, ":", as.numeric(effectiveSize(res[["XDraws"]][, j]))))
  }
  
  # combine chains
  if (is.na(dataX)) {
    dataX <- res[["XDraws"]]
    dataL <- res[["lambdaDraws"]]
  } else {
    dataX <- rbind(dataX, res[["XDraws"]])
    dataL <- rbind(dataL, res[["lambdaDraws"]])
  }
  
  print("========")
}


# plot marginal posterior histograms of Xs
pdf("xiaolinzhuo_1router_t5_informative_1.pdf")
par(mfrow=c(4, 4))
for (j in 1:c) {
  hist(dataX[, j], xlab=paste("X", j, sep=""), main="")
  abline(v=matX[5, reorder][j], col="red") # add a vertical line for true Xs
}
dev.off()


# plot marginal posterior densities of lambdas
pdf("xiaolinzhuo_1router_t5_informative_2.pdf")
par(mfrow=c(4, 4))
for (j in 1:c) {
  densplot(mcmc(dataL[, j]), xlab=paste("L", j, sep=""))
}
dev.off()


# create horizontal boxplots
# data contains too many rows. sample 10% of data for plotting
n <- dim(dataX)[1]
n
tmp <- dataX[sample(n, size=n*0.1), ]
pdf("xiaolinzhuo_1router_t5_informative_5.pdf")
boxplot(tmp, horizontal=T, las=1, names=1:c, main="Informative Prior")
points(x=matX[5, reorder], y=1:c, col="red", pch=16)
dev.off()


