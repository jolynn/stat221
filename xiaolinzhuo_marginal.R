#!/usr/bin/env Rscript

### Stat 221 Pset 3 
### Question 1.5 & 1.6
### Xiaolin Zhuo


library(Rmpfr)


post.marginal <- function(N, Y) {
  # return posterior marginal probability of N (without normalizing constant)
  
  n <- length(Y)
  S <- sum(Y)
  #return(prod(choose(N, Y))/prod(seq(n*N-S+1, n*N+1))/N)
  return(exp(mpfr(sum(log(choose(N, Y))) - sum(log(seq(n*N-S+1, n*N+1))) - log(N), 128)))
}


findConst <- function(Y) {
  # given data Y, return normalizing constant 
  
  s <- 0
  N <- max(Y)
  for (N in max(Y):1000) {
    #print(post.marginal(N, Y))
    s <- s + post.marginal(N, Y)  
  }
  return(1/s)
}



###### impala ######
Y <- read.table("impala.txt", stringsAsFactors=F)
Y
Y <- as.numeric(Y[-1,])
Y

# find posterior prob that N > 100 analytically
C <- findConst(Y) # find normalizing constant
C

prob <- 0
for (N in 101:1000) {
  prob <- prob + C*post.marginal(N, Y)
}
prob


# find posterior prob that N > 100 using MCMC
for (i in 1:10) {
  df <- read.csv(paste("impala_", i, ".csv", sep=""))
  print(mean(df$N>100))
}

# compute mean and sd for 10 posterior simulations
for (i in 1:10) {
  df <- read.csv(paste("impala_", i, ".csv", sep=""))
  print(c(mean(df$N), sd(df$N), mean(df$theta), sd(df$theta)))
}


###### waterbuck ######
Y <- read.table("waterbuck.txt", stringsAsFactors=F)
Y
Y <- as.numeric(Y[-1,])
Y


# find posterior prob that N > 100 analytically
C <- findConst(Y) # find normalizing constant
C
prob <- 0
for (N in 101:1000) {
  prob <- prob + C*post.marginal(N, Y)
}
prob


# find posterior prob that N > 100 using MCMC
for (i in 1:10) {
  df <- read.csv(paste("waterbuck_", i, ".csv", sep=""))
  print(mean(df$N>100))
}


# compute mean and sd for 10 posterior simulations
for (i in 1:10) {
  df <- read.csv(paste("waterbuck_", i, ".csv", sep=""))
  print(c(mean(df$N), sd(df$N), mean(df$theta), sd(df$theta)))
}

