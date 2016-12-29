### stat221 pset 5
### Xiaolin Zhuo


loglik <- function(Y, l, phi, A) {
  sigma <- phi*diag(l**c)
  mat <- A %*% sigma %*% t(A)
  #print(mat)
  return(-nrow(Y)/2*log(det(mat)) - 1/2*sum(apply(Y, 1, function(y) t(y - A %*% l) %*% solve(mat) %*% (y - A %*% l))))
}


tr <- function(m) {
  # trace function 
  return(sum(diag(m)))
}


# E step
Q <- function(logTheta, m_k, v_k, w) {
  l <- exp(logTheta[1:I])
  phi <- exp(logTheta[I+1])
  sigma <- phi*diag(l**c)
  invSigma <- diag(1/(phi*l**c))
  q <- -w/2*(as.numeric(determinant(sigma)$modulus) + tr(invSigma %*% v_k)) - 1/2*sum(apply(m_k, 1, function(m) t(m - l) %*% invSigma %*% (m - l)))
  return(-q) # because optim minimizes
}


locally_iid_EM <- function(data, c, A, w=11, threshold=1e-4, n=1e5, verbose=F) {
  # the algorithm will stop if either increase in loglikelihood < threshold or number of iterations exceeds n
  
  I <- ncol(A)
  J <- nrow(A)
  
  # initialize lambdas
  a0 <- sum(apply(data, 2, mean))/sum(A)
  l0 <- rep(a0, I)
  #l0 <- ipfp(apply(data, 2, mean), A, rgamma(I, 10, 10))
  #all.equal(rep(1, J) %*% (A %*% l0), rep(1, J) %*% apply(data, 2, mean))
  
  # initialize phi
  j <- 1
  phi0 <- var(data[, j])/(sum(data[, j]^2)/w)
  
  if (verbose) {
    print(l0)
    print(phi0)
  }
  
  # start loop from k=2
  k <- 2
  l_k <- l0
  phi_k <- phi0
  
  # store loglikelihood
  ll <- rep(NA, n) 
  ll[1] <- loglik(data, l0, phi0, A)
  
  # store parameters
  #theta <- matrix(NA, nrow=n, ncol=I+1) 
  
  repeat {
  
    # E step - derive Q(theta, theta(k-1))
    sigma_k <- phi_k*diag(l_k**c)
    mat <- solve(A %*% sigma_k %*% t(A))
    m_k <- t(apply(data, 1, function(y) l_k + sigma_k %*% t(A) %*% mat %*% (y - A %*% l_k))) # each row representing E(x(t) | y(t), theta(k-1))
    v_k <- sigma_k - sigma_k %*% t(A) %*% mat %*% A %*% sigma_k  # var(x(t) | y(t), theta(k-1))
    
    # M step - maximize Q
    newLogTheta <- optim(par=c(log(l_k), log(phi_k)), fn=Q, m_k=m_k, v_k=v_k, w=w)$par
    newTheta <- exp(newLogTheta)
    ll[k] <- loglik(data, newTheta[1:I], newTheta[I+1], A)
    
    if (verbose) print(ll[k])
    
    if (k == n) {
      print("reached maximal number of iterations")
      print(paste("last recorded diff in loglik:", ll[k]-ll[k-1]))
      break
    } else if (ll[k]-ll[k-1] < threshold) {
      print("converged")
      break
    } else {
      l_k <- newTheta[1:I]
      phi_k <- newTheta[I+1]
      #theta[k, ] <- c(l_k, phi_k)
      k <- k + 1
    }
  }
  
  return(list(lambda=newTheta[1:I], phi=newTheta[I+1], loglik=na.omit(ll))) # remove empty elements in ll
}

  
# E step with prior in smoothed model
QWithPrior <- function(logTheta, m_k, v_k, prior_m_k, prior_v_k, w) {
  l <- exp(logTheta[1:I])
  phi <- exp(logTheta[I+1])
  
  sigma <- phi*diag(l**c)
  invSigma <- diag(1/(phi*l**c))
  
  q <- -w/2*(as.numeric(determinant(sigma)$modulus) + tr(invSigma %*% v_k)) - 1/2*sum(apply(m_k, 1, function(m) t(m - l) %*% invSigma %*% (m - l)))
  
  # I am not sure of the data type/dimension of (logTheta-prior_m_k)
  if (dim(as.matrix(logTheta-prior_m_k))[1] == 1) {
    logPrior <- -1/2*(as.numeric(determinant(prior_v_k)$modulus)) - 1/2*as.matrix(logTheta-prior_m_k) %*% solve(prior_v_k) %*% t(as.matrix(logTheta-prior_m_k))
  } else { # if dim(as.matrix(logTheta-prior_m_k))[1] == 17
    logPrior <- -1/2*(as.numeric(determinant(prior_v_k)$modulus)) - 1/2*t(as.matrix(logTheta-prior_m_k)) %*% solve(prior_v_k) %*% (as.matrix(logTheta-prior_m_k))
  }
  
  return(-q-logPrior[1, 1]) # because optim minimizes
}


smoothed_EM <- function(data, c, A, V, w=11, threshold=1e-4, n=1e5, verbose=F) {
  # the algorithm will stop if either increase in loglikelihood < threshold or number of iterations exceeds n
  
  I <- ncol(A)
  J <- nrow(A)
  
  # initialize lambdas
  a0 <- sum(apply(data, 2, mean))/sum(A)
  l0 <- rep(a0, I)
  #all.equal(rep(1, J) %*% (A %*% l0), rep(1, J) %*% apply(data, 2, mean))
  
  # initialize phi
  j <- 1
  phi0 <- var(data[, j])/(sum(data[, j]^2)/w)

  # initialize prior of log theta
  est <- read.csv("theta_0.csv")
  logTheta0 <- log(est[1, ])
  sigma0 <- matrix(1e5, nrow=length(logTheta0), ncol=length(logTheta0)) # choose a large sigma0
  
  if (verbose) {
    print(c(log(l0), log(phi0)))
    print(logTheta0)
    print(sigma0)
  }
  
  # start loop from k=2
  k <- 2
  l_k <- l0
  phi_k <- phi0
  prior_m_k <- logTheta0
  prior_v_k <- sigma0 + V
  
  # store loglikelihood
  ll <- rep(NA, n) 
  ll[1] <- loglik(data, l0, phi0, A)
  
  # store parameters
  #theta <- matrix(NA, nrow=n, ncol=I+1) 
  
  repeat {
  
    # E step - derive Q(theta, theta(k-1))
    sigma_k <- phi_k*diag(l_k**c)
    mat <- solve(A %*% sigma_k %*% t(A))
    m_k <- t(apply(data, 1, function(y) l_k + sigma_k %*% t(A) %*% mat %*% (y - A %*% l_k))) # each row representing E(x(t) | y(t), theta(k-1))
    v_k <- sigma_k - sigma_k %*% t(A) %*% mat %*% A %*% sigma_k  # var(x(t) | y(t), theta(k-1))

    # M step - maximize Q
    # I think current vals of theta are good starting points for optim, but
    # optim does not work when initial vals == prior_m_k.
    # So I add 0.01 to initial vals
    res <- optim(par=c(log(l_k+0.01), log(phi_k+0.01)), fn=QWithPrior, m_k=m_k, v_k=v_k, 
                 prior_m_k=prior_m_k, prior_v_k=prior_v_k, w=w, hessian=T)
    newLogTheta <- res$par
    newTheta <- exp(newLogTheta)
    ll[k] <- loglik(data, newTheta[1:I], newTheta[I+1], A)
    hes <- res$hessian # store Hessian matrix for calculating new var-cov matrix of prior
    
    if (verbose) print(ll[k])
    
    if (k == n) {
      print("reached maximal number of iterations")
      print(paste("last recorded diff in loglik:", ll[k]-ll[k-1]))
      break
    } else if (ll[k]-ll[k-1] < threshold) {
      print("converged")
      break
    } else {
      l_k <- newTheta[1:I]
      phi_k <- newTheta[I+1]
      # derive Hessian matrix using function in package numDeriv
      #prior_v_k <- V + solve(-solve(prior_v_k) + hessian(func=function(logTheta) QWithPrior(logTheta, m_k, v_k, prior_m_k, prior_v_k, w), x=newLogTheta))
      prior_m_k <- newLogTheta
      prior_v_k <- V + solve(-solve(prior_v_k) + hes)
      #theta[k, ] <- c(l_k, phi_k)
      k <- k + 1
    }
  }
  
  return(list(lambda=newTheta[1:I], phi=newTheta[I+1], loglik=na.omit(ll))) # remove empty elements in ll
}

