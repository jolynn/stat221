### stat221 problem set 4
### Xiaolin Zhuo


logP2 <- function(X1, X2, L, a, r) { 
  # given vectors X1, X2, L, index a (scalar), and dimension r (scalar)
  # return the number proportional to conditional log likelihood of Xa 
  
  return(X2[a] * log(L[a+r]) - lfactorial(X2[a]) + sum(X1 * log(L[1:r])) - sum(lfactorial(X1)))
}

logPois <- function(x, lambda) {
  # return the number proportional to log likelihood of Poisson(x, lambda)
  return(x * log(lambda) - lfactorial(x))
}

network_mcmc <- function(Y, A, prior, iter=1.2e5, burnin=2e4, verbose=FALSE) {
  
  r <- nrow(A)
  c <- ncol(A)
  
  # decompose matrix A
  reorder <- qr(A)$pivot
  A <- A[, reorder]
  A1 <- A[, 1:r]
  invA1 <- solve(A1)
  A2 <- A[, (r+1):c]
  
  
  # initialize vals of chain
  X0 <- ipfp(Y, A, rgamma(c, 10, 10))
  # all(A %*% X0 == Y)
  L0 <- rgamma(c, shape=prior[["a"]]+X0+1, scale=prior[["b"]]+1)
  
  
  # create empty matrices to store chains
  drawsX <- matrix(NA, nrow=iter+burnin, ncol=c)
  drawsL <- matrix(NA, nrow=iter+burnin, ncol=c) 
  
  # enter initial vals
  drawsX[1, ] <- X0
  drawsL[1, ] <- L0
  
  
  if (verbose) {
    print("X0")
    print(drawsX[1, ])
    print("L0")
    print(drawsL[1, ])
  }
  
  
  # store acceptance rates for each element of X2
  accepts <- rep(0, c-r)
  
  
  # loop for mcmc iterations
  for (t in 2:(iter+burnin)) {
    #if (t %% 1000 == 0) cat(t, " ")
    
    ###### draw X2 from uniform distribution ######
    
    # retrieve lambdas from prev iteration
    L <- drawsL[t-1, ]
    
    # place holders for current Xs
    X1 <- rep(NA, r) 
    X2 <- rep(NA, c-r)
    
    # place holders for newly drawn Xs
    newX1 <- X1
    newX2 <- X2
    
    # iterate through X2 elements
    for (a in 1:(c-r)) { 
      
      # if we are drawing the first element in X2, copy Xs from prev iteration
      # otherwise, copy the most recently updated Xs
      if (is.na(X1[1])) { 
        X1 <- drawsX[t-1, 1:r]
        X2 <- drawsX[t-1, (r+1):c]
        newX1 <- X1
        newX2 <- X2
      } else {
        X1 <- newX1
        X2 <- newX2
      }
      
      # derive bounds for Xa (ath element in X2)
      A0 <- A2
      A0[, a] <- 0
      i <- which(A2[, a] == 1) # index i runs over the set of links whose counts include Xa
      ub <- min(Y[i] - A0[i, ] %*% X2)
      
      # draw a new val of Xa
      newX2[a] <- rpois(1, L[a+r]) 
      #newX2[a]
      
      # compute X1 
      newX1 <- invA1 %*% (Y - A2 %*% newX2)
      #newX1
      
      # redraw Xa if its new value falls out of bounds or corresponding X1 are not all nonnegative
      while (!(newX2[a] >= 0 & newX2[a] <= ub) | !(all(newX1 >= 0))) {
        if (verbose & t %% 1000 == 0) {
          print(paste(t, "redraw X2", a, round(X2[a], 2), "ub:", round(ub, 2), "new val:", round(newX2[a], 2)))
        } 

        newX2[a] <- rpois(1, L[a+r])
        newX1 <- invA1 %*% (Y - A2 %*% newX2)
      }
      

            
      ###### accept or reject new val of Xa ######
      
      # calculate acceptance ratio
      logRatio <- logP2(newX1, newX2, L, a, r) - logP2(X1, X2, L, a, r) + logPois(X2[a], L[a+r]) - logPois(newX2[a], L[a+r])
                  
      if (log(runif(1)) < logRatio) {
        # accept: new val of X2[a] already in newX2
        
        if (t > burnin) accepts[a] <- accepts[a] + 1 # only record acceptance rates after burnin
        
      } else {
        # reject: replace Xa with existing val
        newX2[a] <- X2[a]
        newX1 <- invA1 %*% (Y - A2 %*% newX2)
        
      }
      
    }
    
    # recalcualte new X1 after updating all vals in X2
    newX1 <- invA1 %*% (Y - A2 %*% newX2)
    
    # enter new X1 and new X2 into matrix
    drawsX[t, ] <- c(newX1, newX2)
    
    
    
    ###### draw lambdas from Gamma distribution ######
    drawsL[t, ] <- rgamma(c, shape=prior[["a"]]+drawsX[t, ]+1, scale=prior[["b"]]+1)
    
    if (verbose & t %% 1000 == 0) {
      print(paste("stored at iter", t))
      print(drawsX[t, ])
      print(drawsL[t, ])
    }
    
        
  }
  
  return(list(X=tail(drawsX, iter), L=tail(drawsL, iter), accepts=accepts/iter))
  
}




