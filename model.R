###################################
###### Stat221 final project ######
######      Xiaolin Zhuo     ######
######    Due: 12/19/2016    ######
###################################


library(RSiena)
library(latentnet)
library(sna)
source("simulate.R")


# set true parameters
n <- 50
w <- 4
rate <- 2
dens <- -1.9
rec <- 2
sim.z <- 0.8
rate.z <- 1
lin.z <- 0.1
avalt.z <- 1.2

dis.pair.coef <- -1
dis0.coef <- 0.9


###### ====================================================== ######

# model 1 - basic model

# create empty matrices to store estimates
mat.theta <- matrix(NA, nrow=200, ncol=11)
mat.se <- matrix(NA, nrow=200, ncol=11)
mat.tconv <- matrix(NA, nrow=200, ncol=11)

for (i in 1:200) {
  # simulate network and behavior
  res <- simulate(n, w, rate, dens, rec, sim.z, rate.z, lin.z, avalt.z)
  
  # set up Siena model
  X <- sienaDependent(res$networks, allowOnly=F) 
  Z <- sienaDependent(res$behaviors, type="behavior", allowOnly=F)
  mydata <- sienaDataCreate(X, Z)
  myeff <- getEffects(mydata)
  
  # sink to prevent from printing to screen
  sink("eff1.txt")
  myeff <- includeEffects(myeff, sameX, interaction1="Z")
  myeff <- includeEffects(myeff, name = "Z", quad, include=F)
  myeff <- includeEffects(myeff, name="Z", avAlt, interaction1="X")
  sink()
  
  # set up algorithm and estimate
  model1 <- sienaAlgorithmCreate(projname="twitter_troll", diagonalize=0.2, doubleAveraging=0)
  ans <- siena07(model1, data=mydata, effects=myeff, batch=T, silent=T)
  
  # store estimates
  mat.theta[i, ] <- ans$theta
  mat.se[i, ] <- ans$se
  mat.tconv[i, ] <- ans$tconv
  
  # print total convergence rate
  cat(i, " ", ans$tconv.max[1, 1], "\n")
}


### post-processing of outputs
# create histogram of estimate distribution 
pars <- c(rate, rate, rate, dens, rec, sim.z, rate.z, rate.z, rate.z, lin.z, avalt.z)

pdf("fig1.pdf")
par(mfrow=c(4, 3))
for (i in 1:11) {
  hist(mat.theta[, i], breaks=20)
  abline(v=pars[i], col="red")
  abline(v=mean(mat.theta[, i]), lty=2)
}
dev.off()


# print mean and 95% confidence interval
for (i in 1:11) {
  print(i)
  print(mean(mat.theta[, i]))
  print(as.numeric(quantile(mat.theta[, i], 0.025)))
  print(as.numeric(quantile(mat.theta[, i], 0.975)))
}



###### ====================================================== ######

# simulate latent variables
set.seed(1219)
lat <- matrix(rnorm(n*2), ncol=2)

pdf("latent_space.pdf")
plot(lat, xlim=c(-3, 3), ylim=c(-3, 3), pch=19, xlab="Z1", ylab="Z2")
abline(h=0, lty=2)
abline(v=0, lty=2)
dev.off()

# Euclidean pairwise distance between actors
dis <- as.matrix(dist(lat, method="euclidean"))

# Euclidean distance to origin
dis0 <- apply(lat, 1, function(row) sqrt(row[1]**2 + row[2]**2))

# for exmaple
res <- simulate(n, w, rate, dens, rec, sim.z, rate.z, lin.z, avalt.z, w.lat=T, dis, dis0, dis.pair.coef, dis0.coef)



###### ====================================================== ######

# model 2 - estimate without latent variables

mat2.theta <- matrix(NA, nrow=200, ncol=11)
mat2.se <- matrix(NA, nrow=200, ncol=11)
mat2.tconv <- matrix(NA, nrow=200, ncol=11)

for (i in 1:200) {
  res <- simulate(n, w, rate, dens, rec, sim.z, rate.z, lin.z, avalt.z, w.lat=T, dis, dis0, dis.pair.coef, dis0.coef)
  
  X <- sienaDependent(res$networks, allowOnly=F) 
  Z <- sienaDependent(res$behaviors, type="behavior", allowOnly=F)
  mydata <- sienaDataCreate(X, Z)
  myeff <- getEffects(mydata)
  
  # sink to prevent from printing to screen
  sink("eff2.txt")
  myeff <- includeEffects(myeff, sameX, interaction1="Z")
  myeff <- includeEffects(myeff, name="Z", quad, include=F)
  myeff <- includeEffects(myeff, name="Z", avAlt, interaction1="X")
  sink()
  
  model1 <- sienaAlgorithmCreate(projname="twitter_troll", diagonalize=0.2, doubleAveraging=0)
  ans <- siena07(model1, data=mydata, effects=myeff, batch=T, silent=T)
  
  mat2.theta[i, ] <- ans$theta
  mat2.se[i, ] <- ans$se
  mat2.tconv[i, ] <- ans$tconv
  
  cat(i, " ", ans$tconv.max[1, 1], "\n")
}


### post-processing of outputs
# create histogram of estimate distribution 
pars <- c(rate, rate, rate, dens, rec, sim.z, rate.z, rate.z, rate.z, lin.z, avalt.z)

pdf("fig2.pdf")
par(mfrow=c(4,3))
for (i in 1:11) {
  hist(mat2.theta[, i], breaks=20)
  abline(v=pars[i], col="red")
  abline(v=mean(mat2.theta[, i]), lty=2)
}
dev.off()


# print mean and 95% confidence interval
for (i in 1:11) {
  print(i)
  print(mean(mat.theta[, i]))
  print(as.numeric(quantile(mat.theta[, i], 0.025)))
  print(as.numeric(quantile(mat.theta[, i], 0.975)))
}


###### ====================================================== ######

# model 3 - estimate with latent variables
mat3.theta <- matrix(NA, nrow=200, ncol=13)
mat3.se <- matrix(NA, nrow=200, ncol=13)
mat3.tconv <- matrix(NA, nrow=200, ncol=13)

for (i in 1:200) {
  # simulate network and behavior
  res <- simulate(n, w, rate, dens, rec, sim.z, rate.z, lin.z, avalt.z, w.lat=T, dis, dis0, dis.pair.coef, dis0.coef)
  
  # estimate latent variables from last observation 
  fit <- ergmm(as.network(res$networks[,,4]) ~ euclidean(d=2))
  lat.hat <- cbind(colMeans(fit$sample$Z[,,1]), colMeans(fit$sample$Z[,,2]))
  dis.hat <- as.matrix(dist(lat.hat, method="euclidean"))
  dis0.hat <- apply(lat.hat, 1, function(row) sqrt(row[1]**2 + row[2]**2))
  
  # set up Siena model
  X <- sienaDependent(res$networks, allowOnly=F) 
  Z <- sienaDependent(res$behaviors, type="behavior", allowOnly=F)
  Vdis <- coDyadCovar(dis.hat)
  Vdis0 <- coCovar(dis0.hat)
  
  # create Siena data
  mydata <- sienaDataCreate(X, Z, Vdis, Vdis0)
  myeff <- getEffects(mydata)
  
  # sink to prevent from printing to screen
  sink("eff3.txt")
  myeff <- includeEffects(myeff, sameX, interaction1="Z")
  myeff <- includeEffects(myeff, name = "Z", quad, include=F)
  myeff <- includeEffects(myeff, name="Z", avAlt, interaction1="X")
  myeff <- includeEffects(myeff, X, interaction1="Vdis")
  myeff <- includeEffects(myeff, effFrom, name="Z", interaction1="Vdis0")
  sink()
  
  model1 <- sienaAlgorithmCreate(projname="twitter_troll", diagonalize=0.2, doubleAveraging=0)
  ans <- siena07(model1, data=mydata, effects=myeff, batch=T, silent=T)
  
  mat3.theta[i, ] <- ans$theta
  mat3.se[i, ] <- ans$se
  mat3.tconv[i, ] <- ans$tconv
  
  cat(i, " ", ans$tconv.max[1, 1], "\n")
}


### post-processing of outputs
# create histogram of estimate distribution 
pars <- c(rate, rate, rate, dens, rec, dis.pair.coef, sim.z, rate.z, rate.z, rate.z, lin.z, avalt.z, dis0.coef)

pdf("fig3.pdf")
par(mfrow=c(4,4))
for (i in 1:13) {
  hist(mat3.theta[, i], breaks=20)
  abline(v=pars[i], col="red")
  abline(v=mean(mat3.theta[, i]), lty=2)
}
dev.off()


# print mean and 95% confidence interval
for (i in 1:13) {
  print(i)
  print(mean(mat3.theta[, i]))
  print(as.numeric(quantile(mat3.theta[, i], 0.025)))
  print(as.numeric(quantile(mat3.theta[, i], 0.975)))
}
