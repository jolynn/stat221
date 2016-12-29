### stat221 problem set 4
### Xiaolin Zhuo


library(coda)
library(networkTomography)
source("xiaolinzhuo_mcmc.R")

# get time point (== slurm job task ID)
t <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
t


# read data
#dat <- read.csv("1router_allcount.dat")
data(bell.labs)
A <- bell.labs$A
reorder <- qr(A)$pivot
matX <- bell.labs$X # OD flows in matrix form
matY <- bell.labs$Y # link loads in matrix form
r <- nrow(A)
c <- ncol(A)

# use data from time point t
Y <- matY[t, ]

# draw ten chains
dataX <- NA
dataL <- NA

for (i in 1:10) {
  print(paste("chain", i))
  
  set.seed(i)
  #res <- network_mcmc(Y, A, list(a=0, b=0))
  res <- twMCMC(Y, A, list(a=0.02*matX[t-1, reorder], b=rep(0.02, c)))
  
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


output <- matrix(NA, nrow=c, ncol=6) 
for (j in 1:c) {
  print(j)
  sum(is.na(dataX[, j]))
  sum(is.na(dataL[, j]))
  output[j, ] <- c(median(dataX[, j], na.rm=T), 
                   quantile(dataX[, j], 0.025, na.rm=T), 
                   quantile(dataX[, j], 0.975, na.rm=T),
                   median(dataL[, j], na.rm=T), 
                   quantile(dataL[, j], 0.025, na.rm=T), 
                   quantile(dataL[, j], 0.975, na.rm=T))
  print("======")
}
write.csv(output, paste("informative_time_", t, "_output.csv", sep=""), row.names=F)



df <- read.csv(paste("informative_time_", 2, "_output.csv", sep=""))
df[, 7] <- 2
df
for (t in 3:10) {
  df0 <- read.csv(paste("informative_time_", t, "_output.csv", sep=""))
  df0[, 7] <- t
  df <- rbind(df, df0) 
}


# plot median and CI for X over time
pdf("xiaolinzhuo_1router_all_times_informative_1.pdf")
par(mfrow=c(4, 4))
for (i in 1:c) {
  ind <- seq(from=i, by=16, length.out=9)
  plot(df[ind, 1], type="l", ylim=c(min(df[ind, 2]), max(df[ind, 3])), 
       xlab=paste("X", i, sep=""), ylab="", las=1, xaxt = "n")
  lines(x=1:9, y=df[ind, 2], lty=2)
  lines(x=1:9, y=df[ind, 3], lty=2)
  axis(1, at=1:9, labels=paste("t", 2:10, sep=""))
}
dev.off()


# plot median and CI for L over time
pdf("xiaolinzhuo_1router_all_times_informative_2.pdf")
par(mfrow=c(4, 4))
for (i in 1:c) {
  ind <- seq(from=i, by=16, length.out=9)
  plot(df[ind, 4], type="l", ylim=c(min(df[ind, 2]), max(df[ind, 3])), 
       xlab=paste("lambda", i, sep=""), ylab="", las=1, xaxt = "n")
  lines(x=1:9, y=df[ind, 5], lty=2)
  lines(x=1:9, y=df[ind, 6], lty=2)
  axis(1, at=1:9, labels=paste("t", 2:10, sep=""))
}
dev.off()