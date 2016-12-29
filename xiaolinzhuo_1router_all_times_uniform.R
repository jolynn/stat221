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
  res <- twMCMC(Y, A, list(a=0, b=0))
  
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
write.csv(output, paste("uniform_time_", t, "_output.csv", sep=""), row.names=F)








