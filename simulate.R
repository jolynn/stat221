###################################
###### Stat221 final project ######
######      Xiaolin Zhuo     ######
######    Due: 12/19/2016    ######
###################################


library(RSiena)
library(latentnet)
library(sna)



simulate <- function(n, w, rate, dens, rec, sim.z, rate.z, lin.z, avalt.z, 
                     w.lat=F, dis.pair=NA, dis0=NA, dis.pair.coef=NA, dis0.coef=NA) {
  
  #### modified based on Tom resijders's sample codes at 
  #### https://www.stats.ox.ac.uk/~resijders/siena/NetworkSimulation.R

  # simulates w consecutive network waves according to a stochastic actor-oriented model, with
  # n: number of actors
  # w: number of waves
  # rate: speed at which actors get opportunities to change outgoing ties
  # indeg: in-degree related popularity effect; nodes with higher in-degree will get more incoming ties
  # rec: reciprocity effect
  # sim.z: homophily effect of behavior on network formation
  # rate.z: rate for behavioral changes
  # lin.z: linear tendency of behavioral changes
  # avalt.z: social influence effect based on average behavior in one's neighborhood
  
  #### optional parameters
  # w.lat: whether including latent variables/latent homophily effect in simulations
  # dis.pair: pairwise distance between actors based on estimated positions in latent space
  # dis0: estimated distance to origin on latent space
  # dis.pair.coef: effect of latent space distance on tie formation
  # dis0.coef: effect of latent space position on behavioral changes
  ####
  
  # first, simulate an initial wave of network and behavior
  # next, implement a sequence of two-wave (one-period) simulations
  # in which wave m+1 start at the simulated network from wave m

  # create initial 2-wave data to get suitable data structure
  # this initial network is arbitrary set to have an average degree of 3
  X0 <- matrix(rbinom(n*n, 1, 3/(n-1)), n, n) 
  diag(X0) <- 0
  X1 <- X0
  
  # X0 and X1 should not be identical for use in sienaDependent
  X0[1,2] <- 0
  X0[2,1] <- 1
  X1[1,2] <- 1
  X1[2,1] <- 0
  
  # create network DV
  XX <- array(c(X0, X1), dim=c(n, n, 2))
  X <- sienaDependent(XX, allowOnly=T)
  
  # create behavior DV
  # initially generated from Bernoulli(0.5)
  ZZ <- matrix(as.numeric(runif(n*2) < 0.5), nrow=n)
  Z <- sienaDependent(ZZ, type="behavior", allowOnly=F)
  
  
  if (w.lat) {
    # create covariates
    Vdis <- coDyadCovar(dis.pair)
    Vdis0 <- coCovar(dis0)
    
    # create Siena data
    InitData <- sienaDataCreate(X, Z, Vdis, Vdis0)
  } else {
    # create Siena data
    InitData <- sienaDataCreate(X, Z)
  }
  
  
  InitEff0 <- getEffects(InitData)
  
  # sink to avoid printing to the screen
  sink("eff.txt")
  
  # edit effects and specify effect parameters
  
  # rate parameters are first multiplied by 10 in order to 
  # simulate the first wave from the totally random initial network
  # rate parameters will be set back to input values for subsequent waves
  InitEff0 <- setEffect(InitEff0, Rate, type="rate", initialValue=10*rate)
  InitEff0 <- setEffect(InitEff0, name="Z", Rate, type="rate", initialValue=10*rate.z)
  
  InitEff0 <- setEffect(InitEff0, density, initialValue=dens)
  InitEff0 <- setEffect(InitEff0, recip, initialValue=rec)
  InitEff0 <- setEffect(InitEff0, sameX, interaction1="Z", initialValue=sim.z)
  InitEff0 <- setEffect(InitEff0, name="Z", linear, initialValue=lin.z)
  InitEff0 <- setEffect(InitEff0, name="Z", avAlt, interaction1="X", initialValue=avalt.z)
  
  if (w.lat) {
    InitEff0 <- setEffect(InitEff0, X, interaction1="Vdis", initialValue=dis.pair.coef)
    InitEff0 <- setEffect(InitEff0, effFrom, name="Z", interaction1="Vdis0", initialValue=dis0.coef)
  }
  
  
  # set up algorithm for simulation
  nthree <- sum(InitEff0$include)	+ 5   # n3 should be larger than sum(InitEff0$include)
  InitAlg <- sienaAlgorithmCreate(projname="Init", useStdInits=FALSE, 
                                  cond=FALSE, nsub=0, n3=nthree, simOnly=TRUE)
  
  # simulate first wave
  InitSim <- siena07(InitAlg, data=InitData, eff=InitEff0, returnDeps=TRUE, batch=TRUE, silent=TRUE)
  
  
  # simulate waves from 2 to w
  # create empty network and behavior matrices to store simulated outputs
  Xs <- array(NA, dim=c(n, n, w))
  Zs <- array(NA, dim=c(n, w))
  
  # reset rate parameter values in InitEff
  InitEff <- InitEff0
  InitEff <- setEffect(InitEff, Rate, type="rate", initialValue=rate)
  InitEff <- setEffect(InitEff, name="Z", Rate, type="rate", initialValue=rate.z)
  sink()
  
  for (i in 1:w) {
    # start the loop with the previously simulated network
    
    # transform the previously simulated network from edge list into adjacency matrix
    XXsim <- matrix(0, n, n)
    nsim  <- InitAlg$n3
    XXsim[InitSim$sims[[nsim]][[1]]$X[[1]][, 1:2]] <- InitSim$sims[[nsim]][[1]]$X[[1]][, 3]
    Zsim <- InitSim$sims[[nsim]][[1]]$Z[[1]]
    
    # add simulated network and behavior into output matrices
    Xs[,,i] <- XXsim
    Zs[,i] <- Zsim
    
    if (i < w) {
      # start next simluation from the simulated network and behavior 
      # i.e. put the simulated network and behavior in 1st layer in XX and ZZ
      XX[,,2] <- XX[,,1]
      XX[,,1] <- XXsim
      ZZ[, 2] <- ZZ[, 1]
      ZZ[, 1] <- Zsim
      
      # again, ensure that network and behavior differ between the two waves
      if (identical(XX[,,1], XX[,,2])) {
        XX[1, 2, 1] <- 1 - XX[1, 2, 2]
      }
      
      if (identical(ZZ[, 1], ZZ[, 2])) {
        ZZ[1, 1] <- ifelse((ZZ[1, 1] == 1), 0, 1)
      }
      

      # re-define inputs with newly simulated network and behavior as starting point
      X <- sienaDependent(XX, allowOnly=F)
      Z <- sienaDependent(ZZ, type="behavior", allowOnly=F)
      
      if (w.lat) {
        InitData <- sienaDataCreate(X, Z, Vdis, Vdis0)
      } else {
        InitData <- sienaDataCreate(X, Z)
      }
    
      
      # simulate wave i+1 starting at wave i
      InitSim <- siena07(InitAlg, data=InitData, eff=InitEff, returnDeps=TRUE, batch=TRUE, silent=TRUE)
    }
  }
  
  # print properties of simlated network and behavior
  cat("average degrees ", round(colSums(Xs, dims=2)/n, digits=2), "\n")
  cat("average behavior ", round(colSums(Zs)/n, digits=2), "\n")
  
  # return output as a list:
  # networks is an array of dimension n x n x w
  # behaviors is a matrix of dimension n x w
  return(list(networks=Xs, behaviors=Zs))
  
}



# for exmaple

# parameters
n <- 50
w <- 4
rate <- 2
dens <- -1.9
rec <- 2
sim.z <- 0.8
rate.z <- 1
lin.z <- 0.1
avalt.z <- 1.2

# model 1 - baseline model
res <- simulate(n, w, rate, dens, rec, sim.z, rate.z, lin.z, avalt.z)

for (i in 1:4) {
  print(network.density(as.network(res$networks[,,i])))
}



