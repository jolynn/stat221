#!/usr/bin/env Rscript

### Stat 221 Pset 3 
### Question 1.4
### Xiaolin Zhuo

library(MASS)
library(coda)

proposeN <- function(Y, N, s) {
  # given data Y, N at current iteration, and std of proposal ditribution s, propose a new value for N
  
  return(max(max(Y), round(rnorm(1, mean=N, sd=s), 0))) # N must be integer and >= max[Y]
}


prob <- function(Y, x, mu, s) {
  # return probability of x in Normal with mu and std s
  
  if (x==max(Y)) {
    return(pnorm(max(Y)+0.5, mean=mu, sd=s))
  } else {
    return(pnorm(x+0.5, mean=mu, sd=s) - pnorm(x-0.5, mean=mu, sd=s))
  }
}


getAcceptRatio <- function(Y, N, newN, theta, sigma) {
  # given data Y, N and theta at current iteration, new proposed value newN, and std of proposal ditribution s, 
  # return acceptance ratio
  
  # calculate the factor proportional to log-likelihood
  n <- length(Y)
  fN <- sum(log(choose(N, Y))) + n*N*log(1-theta) - log(N)
  fnewN <- sum(log(choose(newN, Y))) + n*newN*log(1-theta) - log(newN)
  

  return(exp(fnewN-fN)*prob(Y, x=N, mu=newN, s=sigma)/prob(Y, x=newN, mu=N, s=sigma))
}


drawSimul <- function(Y, ndraws=1e5, burnin=1e4, sigma=10) {
  # adj is adjustment factor in proposal distribution in MH step
  
  # set up data structure
  M <- ndraws + burnin
  chain <- matrix(NA, ncol=2, nrow=M)  # first col: N; second col: theta
  
  n <- length(Y)
  N0 <- max(Y) + rpois(1, 10) # start each chain from a different init value
  theta0 <- rbeta(1, sum(Y)+1, n*N0-sum(Y)+1)
  chain[1,] <- c(N0, theta0)
  print(paste("chain starts from N0=", N0, sep=""))
  
  # keep track of acceptance rate and correction rate 
  # (cases in which if a new draw of N < max(Y), set N = max(Y))
  naccept <- 0
  ncorrect <- 0
  
  # loop for mcmc iterations
  for (m in 2:M) {
    
    ###### update N ######
    
    # draw new N
    #newN <- max(max(Y), rpois(1, round(adj*chain[m-1,1], 0))) # N must be integer and >= max[Y]  
    newN <- proposeN(Y, chain[m-1,1], sigma)

    # calculate acceptance ratio
    prob.accept <- getAcceptRatio(Y, chain[m-1,1], newN, chain[m-1, 2], sigma)
    #print(c(m, round(chain[m-1,1], 0), round(chain[m-1,2], 3), prob.accept))
  
    if (runif(1) < prob.accept) {
      chain[m,1] <- newN

      # only count numbers of acceptance and corrections after burnin period
      if (m > burnin) {
        naccept <- naccept + 1  
        ncorrect <- ncorrect + (newN == max(Y))
      }
      
    } else {
      chain[m,1] <- chain[m-1,1]
    }
    
    
    ###### update theta ######
    chain[m,2] <- rbeta(1, sum(Y)+1, n*chain[m,1]-sum(Y)+1)
    
  }
  
  # set up return values
  return(list(N=tail(chain[,1], ndraws), theta=tail(chain[,2], ndraws), 
              burninN=head(chain[,1], burnin), burninTheta=head(chain[,2], burnin),
              accept=naccept/ndraws, correct=ncorrect/ndraws))
  
}


# impala
Y <- read.table("impala.txt", stringsAsFactors=F)
Y
Y <- as.numeric(Y[-1,])
Y


# adjust value of sigma of normal proposal distribution to ensure acceptance rate between 40%-60%

adjustSigma <- function(Y, target=0.5) {
  # run simulations using a range of sigma values and regress acceptance rate on sigma (i.e. step size)
  # find the sigma value that corresponds to acceptance rate of 50%
  
  S <- c(1, seq(10, 100, 10)) # a range of sigma values
  m <- matrix(NA, ncol=3, nrow=length(S)) # column names: sigma, acceptance rate, correction rate
  
  for (i in 1:length(S)) {
    cat(i, " ")
    set.seed(i)
    res <- drawSimul(Y, sigma=S[i])
    m[i,] <- c(S[i], res[["accept"]], res[["correct"]])
  }
  
  print(m)
  model <- lm(m[,2] ~ m[,1])
  print(summary(model))
  
  # solve for sigmal that corresponds to 50% acceptance rate
  print(as.numeric((target - coef(model)[1])/coef(model)[2]))
}

adjustSigma(Y)



# draw 10 chains. draw scatterplots of posterior simulations
for (i in 1:10) {
  cat(i, " ")
  
  set.seed(i)
  res <- drawSimul(Y, sigma=15)
  print(res[["accept"]])
  
  ### plot scatterplot of posterior distribution
  pdf(paste("xiaolinzhuo_ps3_fig", i, ".pdf", sep=""))
  plot(res[["N"]], res[["theta"]], xlab="N", ylab=expression(theta), ylim=c(0,1),
       main=bquote(paste("Posterior simulations of (N, ", theta, "), chain #", .(i))))
  contour(kde2d(res[["N"]], res[["theta"]]), nlevels=10, col="red", lwd=2, add=T) # add contour lines
  dev.off()
  
  ### diagnostic plots: traceplots and autocorrelation plots
  
  # traceplots with burnin
  plot(c(res[["burninN"]], res[["N"]]), type="l", xlab="Iteration", ylab="N", main="Traceplot of N")
  plot(c(res[["burninTheta"]], res[["theta"]]), type="l", xlab="Iteration", ylab=expression(theta), 
       main=expression(paste("Traceplot of ", theta))) 
  
  # autocorrelation plots
  plot(acf(res[["N"]], lag.max=5000), main="Autocorrletion plot of N")
  plot(acf(res[["theta"]], lag.max=5000), main=expression(paste("Autocorrletion plot of ", theta)))
  
  # save posterior distribution
  write.csv(data.frame(N=res[["N"]], theta=res[["theta"]]), file=paste("impala_", i, ".csv", sep=""), row.names=F)
}




# waterbuck
Y <- read.table("waterbuck.txt", stringsAsFactors=F)
Y
Y <- as.numeric(Y[-1,])
Y

# adjust sigma
adjustSigma(Y)

# draw 10 chains. draw scatterplots of posterior simulations
for (i in 1:10) {
  cat(i, " ")
  
  set.seed(i)
  res <- drawSimul(Y, sigma=40)
  print(res[["accept"]])
  
  ### plot scatterplot of posterior distribution
  pdf(paste("xiaolinzhuo_ps3_fig", i+10, ".pdf", sep=""))
  plot(res[["N"]], res[["theta"]], xlab="N", ylab=expression(theta), ylim=c(0,1),
       main=bquote(paste("Posterior simulations of (N, ", theta, "), chain #", .(i))))
  contour(kde2d(res[["N"]], res[["theta"]]), nlevels=10, col="red", lwd=2, add=T) # add contour lines
  dev.off()
    
  ### diagnostic plots: traceplots and autocorrelation plots
  
  # traceplots with burnin
  plot(c(res[["burninN"]], res[["N"]]), type="l", xlab="Iteration", ylab="N", main="Traceplot of N")
  plot(c(res[["burninTheta"]], res[["theta"]]), type="l", xlab="Iteration", ylab=expression(theta), 
       main=expression(paste("Traceplot of ", theta))) 
  
  # autocorrelation plots
  plot(acf(res[["N"]], lag.max=5000), main="Autocorrletion plot of N")
  plot(acf(res[["theta"]], lag.max=5000), main=expression(paste("Autocorrletion plot of ", theta)))
  
  # save posterior distribution
  write.csv(data.frame(N=res[["N"]], theta=res[["theta"]]), file=paste("waterbuck_", i, ".csv", sep=""), row.names=F)
}





