### stat221 pset 5
### Xiaolin Zhuo


source("xiaolinzhuo_functions.R")


###### 1.1 replicate figure 2 of link loads ######
df <- read.csv("1router_allcount.dat")
names(df)
unique(df$nme)

nodes <- c("corp", "local", "switch", "fddi")

pdf("xiaolinzhuo_fig2.pdf")
par(mfrow=c(4, 1), mai=c(0.1, 0.5, 0.25, 0.5), oma=c(5, 1, 4, 1)) 
for (i in 1:length(nodes)) {
  plot(df[df$nme==paste("src", nodes[i]), "value"], type="l", col="red", 
       ylim=c(0, 1e6), main=nodes[i], xlab="", ylab="", xaxt="n", yaxt="n")
  lines(df[df$nme==paste("dst", nodes[i]), "value"], col="blue")
  
  # add y axis
  if (i %% 2 == 1) { # i is odd number
    axis(4, at=seq(0, 1e6, 2e5), las=2, 
         labels=c("0", "200K", "400K", "600K", "800K", "1M"))
  } else {
    axis(2, at=seq(0, 1e6, 2e5), las=2, 
         labels=c("0", "200K", "400K", "600K", "800K", "1M"))
  }
  
  # add x axis
  if (i==1) {
    axis(3, at=seq(0, 288, 48), labels=NA)
  } else if (i==4) {
    axis(1, at=seq(0, 288, 48), labels=seq(0, 24, 4))
  }
  
  
}

# add legend
# reset coordinates, allow plotting outside of plot region
op <- par(usr=c(0, 1, 0, 1), xpd=NA) 
legend(0.35, 5.5, c("origin", "destination"), lty=1, col=c("red", "blue"), 
       ncol=2, bty="n", lwd=2) 

# add title
mtext("hour of day", side=1, line=3, cex=1.1)

dev.off()


###### ====================================================================

###### 1.2 replicate figure 4 ######
unique(df$time)
# (02/22/99 11:32:43)
# (02/22/99 15:32:42)


pdf("xiaolinzhuo_fig4_1router.pdf")
par(mfrow=c(1, 2), mai=c(0.2, 0.2, 0.2, 0), oma=c(4, 3, 1.5, 1.5)) 
nodes <- c("fddi", "switch", "local", "corp")

# left plot for time 11:30
ind <- which(levels(df$time)=="(02/22/99 11:32:43)")
levels(df$time)[seq(ind-5, ind+5)]

# calculate log10(mean) and log10(var)
dat1 <- matrix(NA, nrow=8, ncol=2)
for (i in 1:4) {
  d <- df[(df$nme==paste("src", nodes[i]) & df$time %in% levels(df$time)[seq(ind-5, ind+5)]), ]
  dat1[i, ] <- c(log10(mean(d$value)), log10(var(d$value)))
  e <- df[(df$nme==paste("dst", nodes[i]) & df$time %in% levels(df$time)[seq(ind-5, ind+5)]), ]
  dat1[i+4, ] <- c(log10(mean(e$value)), log10(var(e$value)))
}

# plot
plot(dat1[, 1], dat1[, 2], type="n", xlab="", ylab="")
text(4.9,  10.8, "time 11:30", font=2)
abline(lm(dat1[, 2] ~ dat1[, 1]))
text(dat1[, 1], dat1[, 2], labels=1:8, cex= 0.7)
axis(3, at=seq(4.2, 5.4, 0.2), labels=NA)
abline(a=mean(dat1[, 2])-mean(dat1[, 1]), b=1, lty=2) 
abline(a=mean(dat1[, 2])-2*mean(dat1[, 1]), b=2, lty=2)


# right plot for time 15:30
ind <- which(levels(df$time)=="(02/22/99 15:32:42)")
levels(df$time)[seq(ind-5, ind+5)]

# calculate log10(mean) and log10(var)
dat2 <- matrix(NA, nrow=8, ncol=2)
for (i in 1:4) {
  d <- df[(df$nme==paste("src", nodes[i]) & df$time %in% levels(df$time)[seq(ind-5, ind+5)]), ]
  dat2[i, ] <- c(log10(mean(d$value)), log10(var(d$value)))
  e <- df[(df$nme==paste("dst", nodes[i]) & df$time %in% levels(df$time)[seq(ind-5, ind+5)]), ]
  dat2[i+4, ] <- c(log10(mean(e$value)), log10(var(e$value)))
}

# plot
plot(dat2[, 1], dat2[, 2], type="n", xlab="", ylab="", xaxt="n", yaxt="n", 
     ylim=range(dat1[, 2]), xlim=range(dat1[, 1]))
text(4.9,  10.8, "time 15:30", font=2)
abline(lm(dat2[, 2] ~ dat2[, 1]))
text(dat2[, 1], dat2[, 2], labels=1:8, cex= 0.7)
axis(1, at=seq(4.2, 5.4, 0.2), labels=NA)
axis(3, at=seq(4.2, 5.4, 0.2), labels=seq(4.2, 5.4, 0.2))
axis(4, at=c(7, 8, 9, 10), labels=NA)
abline(a=mean(dat2[, 2])-mean(dat2[, 1]), b=1, lty=2) 
abline(a=mean(dat2[, 2])-2*mean(dat2[, 1]), b=2, lty=2)
mtext("log10(var)", side=2, outer=T, line=1.5)
mtext("log10(mean)", side=1, outer=T, line=1.5)
dev.off()



###### ====================================================================

###### 1.3-1.4 implement EM algorithm for locally iid model ######
######           ONLY THIS PART RUN ON ODYSSEY              ######

if (file.exists("matrixA.csv")) {
  A <- as.matrix(read.csv("matrixA.csv"))
} else {
  library(networkTomography)
  data(bell.labs)
  A <- bell.labs$A
}

if (file.exists("matrixY.csv")) {
  Y <- as.matrix(read.csv("matrixY.csv"))
} else {
  library(networkTomography)
  data(bell.labs)
  Y <- bell.labs$Y
}



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
write.csv(theta, paste("theta_", id, ".csv", sep=""), row.names=F)



###### ====================================================================

###### 1.4 replicate figure 5 ######

# get estimated lambdas
theta <- NA
for (i in seq(0, 287, 10)) {
  if (i == 0) {
    theta <- read.csv(paste("theta_", i, ".csv", sep=""))
  } else {
    theta <- rbind(theta, read.csv(paste("theta_", i, ".csv", sep="")))
  }
}

head(theta)

# add 5 empty rows to top and bottom (because I only estimated time points 6-282)
nrow(theta)
theta <- rbind(matrix(NA, nrow=5, ncol=17), theta, matrix(NA, nrow=5, ncol=17))
nrow(theta)


names <- c("fddi->fddi", "fddi->switch", "fddi->local", "fddi->corp", 
           "switch->fddi", "switch->switch", "switch->local", "switch->corp", 
           "local->fddi", "local->switch", "local->local", "local->corp", 
           "corp->fddi", "corp->switch", "corp->local", "corp->corp")

orderInFig <- c("corp->fddi", "corp->switch", "corp->local", "corp->corp",
                "local->fddi", "local->switch", "local->local", "local->corp", 
                "switch->fddi", "switch->switch", "switch->local", "switch->corp", 
                "fddi->fddi", "fddi->switch", "fddi->local", "fddi->corp")

nodes <- c("corp", "local", "switch", "fddi")

pdf("xiaolinzhuo_fig5.pdf")
par(mfrow=c(5, 5), mai=c(0, 0, 0, 0), oma=c(5, 5, 5, 4))
###  row 1  ###
# plot destination
titles <- c("destination fddi", "destination switch", "destination local", "destination corp")
for (i in 1:4) {
  # plot estimates
  ind <- which(names %in% orderInFig[seq(i, 16, 4)])
  plot(apply(theta[, ind], 1, sum), type="l", ylim=c(0, 1e6), main="", xlab="", ylab="", xaxt="n", yaxt="n")
  
  # plot smoothed observations
  lines(smooth.spline(df[df$nme==paste("dst", nodes[5-i]), "value"], df=50), col="red")
  
  # add x axis
  if (i %% 2 == 1) axis(3, at=seq(0, 288, 48), labels=seq(0, 24, 4))
  if (i %% 2 == 0) axis(3, at=seq(0, 288, 48), labels=NA)
  
  # add y axis
  if (i == 1) axis(2, at=seq(0, 1e6, 2e5), las=2, labels=c("0", "200K", "400K", "600K", "800K", "1M"))
  
  # add title
  text(130, 9.7e5, titles[i], font=2, col="red")
}

# plot total
plot(apply(theta, 1, sum), type="l", ylim=c(0, 1e6), main="", xlab="", ylab="", xaxt="n", yaxt="n")

# create a dataframe with 4 columns of destination flows
d <- df[df$nme==paste("dst", nodes[1]), "value"]
for (i in 2:4) d <- cbind(d, df[df$nme==paste("dst", nodes[i]), "value"])

lines(smooth.spline(apply(d, 1, sum), df=50), col="red")
axis(3, at=seq(0, 288, 48), labels=seq(0, 24, 4)) # add x axis
axis(4, at=seq(0, 1e6, 2e5), las=2, labels=NA)   # add y axis
text(140, 9.7e5, "total", font=2, col="red")   # add title


###  row 2-4  ###
for (r in 1:3) {
  for (i in 1:4) {
    j <- which(names %in% orderInFig[i+(r-1)*4])
    plot(theta[, j], type="l", ylim=c(0, 1e6), main="", xlab="", ylab="", xaxt="n", yaxt="n")
    
    # add y axis
    if (i == 1) {
      if (r %% 2 == 1) axis(2, at=seq(0, 1e6, 2e5), las=2, labels=NA)
      else axis(2, at=seq(0, 1e6, 2e5), las=2, labels=c("0", "200K", "400K", "600K", "800K", "1M"))
    }
    
    # add title
    text(140, 9.7e5, orderInFig[i+(r-1)*4], font=2)
  }
  
  # plot origin corp
  ind <- which(names %in% orderInFig[1:4+(r-1)*4])
  plot(apply(theta[, ind], 1, sum), type="l", ylim=c(0, 1e6), main="", xlab="", ylab="", xaxt="n", yaxt="n")
  lines(smooth.spline(df[df$nme==paste("src", nodes[r]), "value"], df=50), col="red")
  if (r %% 2 == 1) axis(4, at=seq(0, 1e6, 2e5), las=2, labels=c("0", "200K", "400K", "600K", "800K", "1M")) # add y axis
  text(140, 9.7e5, paste("origin", nodes[r]), font=2, col="red") # add title
  
}

###  row 5  ###
for (i in 13:16) {
  j <- which(names == orderInFig[i])
  plot(theta[, j], type="l", ylim=c(0, 1e6), main="", xlab="", ylab="", xaxt="n", yaxt="n")
  
  # add x axis
  if (i %% 2 == 1) axis(1, at=seq(0, 288, 48), labels=seq(0, 24, 4))
  if (i %% 2 == 0) axis(1, at=seq(0, 288, 48), labels=NA)
  
  # add y axis
  if (i == 1) axis(2, at=seq(0, 1e6, 2e5), las=2, labels=c("0", "200K", "400K", "600K", "800K", "1M"))
  
  # add title
  text(140, 9.7e5, orderInFig[i], font=2)
}

# plot origin fddi
ind <- which(names %in% orderInFig[13:16])
plot(apply(theta[, ind], 1, sum), type="l", ylim=c(0, 1e6), main="", xlab="", ylab="", xaxt="n", yaxt="n")
lines(smooth.spline(df[df$nme==paste("src", nodes[4]), "value"], df=50), col="red")
axis(1, at=seq(0, 288, 48), labels=seq(0, 24, 4)) # add x axis
axis(4, at=seq(0, 1e6, 2e5), las=2, labels=NA) # add y axis
text(140, 9.7e5, paste("origin", nodes[4]), font=2, col="red") # add title


# add main X axis
mtext("hour of day", side=1, outer=T, line=3, cex=1.1)

# add main Y axis
mtext("bytes/sex", side=2, outer=T, line=3, cex=1.1)
dev.off()



###### ====================================================================

###### 1.5-1.6 implement SMOOTHETHED EM algorithm for 1-router datas ######
######           ONLY THIS PART RUN ON ODYSSEY                       ######

#data(bell.labs)
#A <- bell.labs$A
#Y <- bell.labs$Y
A <- as.matrix(read.csv("matrixA.csv"))
Y <- as.matrix(read.csv("matrixY.csv"))
c <- 2
w <- 11
h <- 5 # w = 2*h + 1

# empirically derive prior var-cov matrix V
if (file.exists("matrixV.csv")) {
  V <- as.matrix(read.csv("matrixV.csv"))
} else {
  res <- NA
  for (i in seq(0, 287, 10)) {
    if (i == 0) {
      mat <- read.csv(paste("theta_", i, ".csv", sep=""))
    } else {
      mat <- rbind(mat, read.csv(paste("theta_", i, ".csv", sep="")))
    }
  }
  
  V <- cov(mat)
  write.csv(V, "matrixV.csv", row.names=F)
}


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
write.csv(theta, paste("smoothed_theta_", id, ".csv", sep=""), row.names=F)



###### ====================================================================

###### 1.6 replicate figure 6 ######

# get estimated lambdas
theta <- NA
for (i in seq(0, 287, 10)) {
  if (i == 0) {
    theta <- read.csv(paste("smoothed_theta_", i, ".csv", sep=""))
  } else {
    theta <- rbind(theta, read.csv(paste("smoothed_theta_", i, ".csv", sep="")))
  }
}

head(theta)

# add 5 empty rows to top and bottom (because I only estimated time points 6-282)
nrow(theta)
theta <- rbind(matrix(NA, nrow=5, ncol=17), theta, matrix(NA, nrow=5, ncol=17))
nrow(theta)



thetaPrior <- NA
for (i in seq(0, 287, 10)) {
  if (i == 0) {
    thetaPrior <- read.csv(paste("theta_", i, ".csv", sep=""))
  } else {
    thetaPrior <- rbind(thetaPrior, read.csv(paste("theta_", i, ".csv", sep="")))
  }
}

# add 5 empty rows to top and bottom (because I only estimated time points 6-282)
nrow(thetaPrior)
thetaPrior <- rbind(matrix(NA, nrow=5, ncol=17), thetaPrior, matrix(NA, nrow=5, ncol=17))
nrow(thetaPrior)




names <- c("fddi->fddi", "fddi->switch", "fddi->local", "fddi->corp", 
           "switch->fddi", "switch->switch", "switch->local", "switch->corp", 
           "local->fddi", "local->switch", "local->local", "local->corp", 
           "corp->fddi", "corp->switch", "corp->local", "corp->corp")

orderInFig <- c("corp->fddi", "corp->switch", "corp->local", "corp->corp",
                "local->fddi", "local->switch", "local->local", "local->corp", 
                "switch->fddi", "switch->switch", "switch->local", "switch->corp", 
                "fddi->fddi", "fddi->switch", "fddi->local", "fddi->corp")

nodes <- c("corp", "local", "switch", "fddi")

pdf("xiaolinzhuo_fig6_1router.pdf")
par(mfrow=c(5, 5), mai=c(0, 0, 0, 0), oma=c(5, 5, 5, 4))
###  row 1  ###
# plot destination
titles <- c("destination fddi", "destination switch", "destination local", "destination corp")
for (i in 1:4) {
  # plot estimates
  ind <- which(names %in% orderInFig[seq(i, 16, 4)])
  plot(apply(theta[, ind], 1, sum), type="l", ylim=c(0, 5e4), main="", xlab="", ylab="", xaxt="n", yaxt="n")
  
  # plot smoothed observations
  lines(apply(thetaPrior[, ind], 1, sum), col="red")
  
  # add x axis
  if (i %% 2 == 1) axis(3, at=seq(0, 288, 48), labels=seq(0, 24, 4))
  if (i %% 2 == 0) axis(3, at=seq(0, 288, 48), labels=NA)
  
  # add y axis
  if (i == 1) axis(2, at=seq(0, 5e4, 1e4), las=2, labels=c("0", "10K", "20K", "30K", "40K", "50K"))
  
  # add title
  text(130, 4.8e4, titles[i], font=2)
}

# plot total
plot(apply(theta, 1, sum), type="l", ylim=c(0, 5e4), main="", xlab="", ylab="", xaxt="n", yaxt="n")
lines(apply(thetaPrior, 1, sum), col="red")
axis(3, at=seq(0, 288, 48), labels=seq(0, 24, 4)) # add x axis
axis(4, at=seq(0, 5e4, 1e4), las=2, labels=NA)   # add y axis
text(140, 4.8e4, "total", font=2)   # add title


###  row 2-4  ###
for (r in 1:3) {
  for (i in 1:4) {
    j <- which(names %in% orderInFig[i+(r-1)*4])
    plot(theta[, j], type="l", ylim=c(0, 5e4), main="", xlab="", ylab="", xaxt="n", yaxt="n")
    lines(thetaPrior[, j], col="red")
    
    # add y axis
    if (i == 1) {
      if (r %% 2 == 1) axis(2, at=seq(0, 5e4, 1e4), las=2, labels=NA)
      else axis(2, at=seq(0, 5e4, 1e4), las=2, labels=c("0", "10K", "20K", "30K", "40K", "50K"))
    }
    
    # add title
    text(140, 4.8e4, orderInFig[i+(r-1)*4], font=2)
  }
  
  # plot origin corp
  ind <- which(names %in% orderInFig[1:4+(r-1)*4])
  plot(apply(theta[, ind], 1, sum), type="l", ylim=c(0, 5e4), main="", xlab="", ylab="", xaxt="n", yaxt="n")
  lines(apply(thetaPrior[, ind], 1, sum), col="red")
  if (r %% 2 == 1) axis(4, at=seq(0, 5e4, 1e4), las=2, labels=c("0", "10K", "20K", "30K", "40K", "50K")) # add y axis
  text(140, 4.8e4, paste("origin", nodes[r]), font=2) # add title
  
}

###  row 5  ###
for (i in 13:16) {
  j <- which(names == orderInFig[i])
  plot(theta[, j], type="l", ylim=c(0, 5e4), main="", xlab="", ylab="", xaxt="n", yaxt="n")
  lines(thetaPrior[, j], col="red")
  
  # add x axis
  if (i %% 2 == 1) axis(1, at=seq(0, 288, 48), labels=seq(0, 24, 4))
  if (i %% 2 == 0) axis(1, at=seq(0, 288, 48), labels=NA)
  
  # add y axis
  if (i == 1) axis(2, at=seq(0, 5e4, 1e4), las=2, labels=c("0", "10K", "20K", "30K", "40K", "50K"))
  
  # add title
  text(140, 4.8e4, orderInFig[i], font=2)
}

# plot origin fddi
ind <- which(names %in% orderInFig[13:16])
plot(apply(theta[, ind], 1, sum), type="l", ylim=c(0, 5e4), main="", xlab="", ylab="", xaxt="n", yaxt="n")
lines(apply(thetaPrior[, ind], 1, sum), col="red")
axis(1, at=seq(0, 288, 48), labels=seq(0, 24, 4)) # add x axis
axis(4, at=seq(0, 5e4, 1e4), las=2, labels=NA) # add y axis
text(140, 4.8e4, paste("origin", nodes[4]), font=2) # add title


# add main X axis
mtext("hour of day", side=1, outer=T, line=3, cex=1.1)

# add main Y axis
mtext("bytes/sex", side=2, outer=T, line=3, cex=1.1)
dev.off()
