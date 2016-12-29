### stat221 pset 5
### Xiaolin Zhuo


source("xiaolinzhuo_functions.R")


###### 1.9 fit locally iid model to 2-router data ######
######     ONLY THIS PART RUN ON ODYSSEY          ######


df2 <- read.csv("2router_linkcount.dat")
unique(df2$nme)

nodes <- c("router5", "r4-local", "switch", "r4-others", "gw1", "gw2", "gw3", "gw-others")

# create matrix 'names' in the same order as the labels in figure 9
names <- matrix(paste("router5 >", nodes), nrow=1)
for (n in nodes[-1]) {
  names <- rbind(paste(n, ">", nodes), names)
}
names

# A is 15 x 64 matrix
if (file.exists("matrixA_2router.csv")) {
  A <- as.matrix(read.csv("matrixA_2router.csv"))
} else {
  A <- matrix(c(rep(0, 56), rep(1, 8), 
                rep(0, 48), rep(1, 8), rep(0, 8), 
                rep(0, 40), rep(1, 8), rep(0, 16), 
                rep(0, 32), rep(1, 8), rep(0, 24), 
                rep(0, 24), rep(1, 8), rep(0, 32), 
                rep(0, 16), rep(1, 8), rep(0, 40), 
                rep(0, 8), rep(1, 8), rep(0, 48),
                rep(1, 8), rep(0, 56),
                as.numeric(1:64 %in% seq(1, 64, 8)),
                as.numeric(1:64 %in% seq(2, 64, 8)), 
                as.numeric(1:64 %in% seq(3, 64, 8)), 
                as.numeric(1:64 %in% seq(4, 64, 8)), 
                as.numeric(1:64 %in% seq(5, 64, 8)), 
                as.numeric(1:64 %in% seq(6, 64, 8)), 
                as.numeric(1:64 %in% seq(7, 64, 8))), nrow=15, byrow=T)
  write.csv(A, "matrixA_2router.csv", row.names=F)
}

# data Y
if (file.exists("matrixY_2router.csv")) {
  Y <- as.matrix(read.csv("matrixY_2router.csv"))
} else {
  Y <- matrix(df2[1:15, "value"], nrow=1)
  for (i in seq(17, nrow(df2), 16)) {
    Y <- rbind(Y, df2[i:(i+14), "value"])
  }
  write.csv(Y, "matrixY_2router.csv", row.names=F)
}
dim(Y)
head(Y)


c <- 2
w <- 11
h <- 5 # w = 2*h + 1


# get slurm job task ID
id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#id <- 0
set.seed(id)

# fits model for 10 time points per job
# store parameters (lambdas, phi)
theta <- matrix(NA, nrow=10, ncol=ncol(A)+1)

for (i in 1:10) {
  cat(i, " ")
  t <- id + i
  
  # ensure the time point corresponds to a valid 11-unit time window
  if (t >= h+1 & t <= nrow(Y)-h) {
    data <- Y[(t-h):(t+h), ] # determine time window
    res <- locally_iid_EM(data, c, A)
    theta[i, ] <- c(res$lambda, res$phi)
  }
}


# remove empty rows in theta - the first job and the last job will have some empty rows
theta <- theta[rowSums(is.na(theta)) == 0, ]


# export theta to csv
write.csv(theta, paste("theta_2router_", id, ".csv", sep=""), row.names=F)




###### 1.9 fit refined model to 2-router data ######
######   ONLY THIS PART RUN ON ODYSSEY        ######

# A is 15 x 64 matrix
if (file.exists("matrixA_2router.csv")) {
  A <- as.matrix(read.csv("matrixA_2router.csv"))
} else {
  A <- matrix(c(rep(0, 56), rep(1, 8), 
                rep(0, 48), rep(1, 8), rep(0, 8), 
                rep(0, 40), rep(1, 8), rep(0, 16), 
                rep(0, 32), rep(1, 8), rep(0, 24), 
                rep(0, 24), rep(1, 8), rep(0, 32), 
                rep(0, 16), rep(1, 8), rep(0, 40), 
                rep(0, 8), rep(1, 8), rep(0, 48),
                rep(1, 8), rep(0, 56),
                as.numeric(1:64 %in% seq(1, 64, 8)),
                as.numeric(1:64 %in% seq(2, 64, 8)), 
                as.numeric(1:64 %in% seq(3, 64, 8)), 
                as.numeric(1:64 %in% seq(4, 64, 8)), 
                as.numeric(1:64 %in% seq(5, 64, 8)), 
                as.numeric(1:64 %in% seq(6, 64, 8)), 
                as.numeric(1:64 %in% seq(7, 64, 8))), nrow=15, byrow=T)
  write.csv(A, "matrixA_2router.csv", row.names=F)
}

# data Y
if (file.exists("matrixY_2router.csv")) {
  Y <- as.matrix(read.csv("matrixY_2router.csv"))
} else {
  Y <- matrix(df2[1:15, "value"], nrow=1)
  for (i in seq(17, nrow(df2), 16)) {
    Y <- rbind(Y, df2[i:(i+14), "value"])
  }
  write.csv(Y, "matrixY_2router.csv", row.names=F)
}
dim(Y)
head(Y)


# empirically derive prior var-cov matrix V
if (file.exists("matrixV_2router.csv")) {
  V <- as.matrix(read.csv("matrixV_2router.csv"))
} else {
  res <- NA
  for (i in seq(0, 287, 10)) {
    if (i == 0) {
      mat <- read.csv(paste("theta_2router_", i, ".csv", sep=""))
    } else {
      mat <- rbind(mat, read.csv(paste("theta_2router_", i, ".csv", sep="")))
    }
  }
  
  V <- cov(mat)
  write.csv(V, "matrixV_2router.csv", row.names=F)
}



c <- 2
w <- 11
h <- 5 # w = 2*h + 1


# get slurm job task ID
id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#id <- 0
set.seed(id)

# fit model to 10 time points per job
# store parameters (lambdas, phi)
theta <- matrix(NA, nrow=10, ncol=ncol(A)+1)

for (i in 1:10) {
  cat(i, " ")
  t <- id + i
  
  # ensure the time point corresponds to a valid 11-unit time window
  if (t >= h+1 & t <= nrow(Y)-h) {
    data <- Y[(t-h):(t+h), ] # determine time window
    res <- smoothed_EM(data, c, A, V)
    theta[i, ] <- c(res$lambda, res$phi)
  }
}

# remove empty rows in theta - the first job and the last job will have some empty rows
theta <- theta[rowSums(is.na(theta)) == 0, ]

# export theta to csv
write.csv(theta, paste("smoothed_theta_2router_", id, ".csv", sep=""), row.names=F)



###### ====================================================================

###### 1.9 replicate figure 6 to 2-router data ######
# estimates from refined model is unavailable
# plot only theta prior


thetaPrior <- NA
for (i in seq(0, 287, 10)) {
  if (i == 0) {
    thetaPrior <- read.csv(paste("theta_2router_", i, ".csv", sep=""))
  } else {
    thetaPrior <- rbind(thetaPrior, read.csv(paste("theta_2router_", i, ".csv", sep="")))
  }
}

# add 5 empty rows to top and bottom (because I only estimated time points 6-282)
nrow(thetaPrior)
thetaPrior <- rbind(matrix(NA, nrow=5, ncol=65), thetaPrior, matrix(NA, nrow=5, ncol=65))
nrow(thetaPrior)


nodes <- c("router5", "r4-local", "switch", "r4-others", "gw1", "gw2", "gw3", "gw-others")

# create matrix 'names' in the same order as the labels in figure 9
names <- matrix(paste("router5 >", nodes), nrow=1)
for (n in nodes[-1]) {
  names <- rbind(paste(n, ">", nodes), names)
}
names


pdf("xiaolinzhuo_fig6_2router.pdf")
par(mfrow=c(8, 8), mai=c(0, 0, 0, 0), oma=c(5, 5, 5, 4))

for (i in 1:64) {
  plot(thetaPrior[, i], type="l", ylim=c(0, 1e5), main="", xlab="", ylab="", xaxt="n", yaxt="n")
  
  # add x axis
  if (i > 56) {
    if (i %% 2 == 0) axis(1, at=seq(0, 288, 48), labels=seq(0, 24, 4))
    else axis(1, at=seq(0, 288, 48), labels=NA)
  }  
  
  # add y axis
  if (i %% 8 == 1) {
    if (i %in% c(9, 25, 41, 57)) axis(2, at=seq(0, 1e5, 2e4), las=2, labels=c("0", "20K", "40K", "60K", "80K", "100K"))
    else axis(2, at=seq(0, 1e5, 2e4), las=2, labels=NA)
  }
  
  # add title
  text(140, 0.95e5, as.vector(names)[i], font=2, cex=0.6)
}


# add main X axis
mtext("hour of day", side=1, outer=T, line=3, cex=1.1)

# add main Y axis
mtext("bytes/sex", side=2, outer=T, line=3, cex=1.1)
dev.off()


