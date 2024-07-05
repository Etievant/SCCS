## -----------------------------------------------------------------------------
## Description: This script illustrates Scenario 1 with violation of Condition 
##              (i) presented in Section 5 in Etievant, Gail and Follmann (2024)
##
##              The simulation is in Section 7.1 in the Main Document
## -----------------------------------------------------------------------------

### load packages --------------------------------------------------------------
library(tidyr)
library(parallel)
library(ggplot2)

### Set fixed parameters -------------------------------------------------------
n         <- 10^6             # (overall) population size
pW        <- 0.6              # P(W = 1)
alpha.T1  <- log(0.00006)     # log baseline risk in period 1 
alpha.T2  <- log(0.00006)     # log baseline risk in period 2 
gamma.T1  <- 0.3              # we use g_1: w -> gamma.T1 * w
gamma.T2  <- 0.3              # we use g_2: w -> gamma.T2 * w
eta.T1    <- 0.2              # we use h_1: w -> eta.T1 * w
#eta.T2                       # we use h_2: w -> eta.T2 * w
beta.T1   <- 0.6              
#beta.T2    
p.Z.W1    <- 0.2              # P(Z = 1 | W = 1)
p.Z.W0    <- 0.35             # P(Z = 1 | W = 0)

### Compute counterfactual probabilities analytically --------------------------
P.T1.doZ0.W1    <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * 1 + eta.T1 * 0 * 1) # P(T^{Z=0} = 1 | W = 1)
P.T1.doZ0.W0    <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * 0 + eta.T1 * 0 * 0) # P(T^{Z=0} = 1 | W = 0)
P.T1.doZ1.W1    <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 1 + eta.T1 * 1 * 1) # P(T^{Z=1} = 1 | W = 1)
P.T1.doZ1.W0    <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 0 + eta.T1 * 1 * 0) # P(T^{Z=1} = 1 | W = 0)

P.T1.doZ1       <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 1 + eta.T1 * 1 * 1) * pW + exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 0 + eta.T1 * 1 * 0) * (1 - pW) # P(T^{Z=1} = 1)
P.T1.doZ0       <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * 1 + eta.T1 * 0 * 1) * pW + exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * 0 + eta.T1 * 0 * 0) * (1 - pW) # P(T^{Z=0} = 1)

exp.beta        <- P.T1.doZ1 / P.T1.doZ0        # marginal causal effect
exp.beta.W1     <- P.T1.doZ1.W1 / P.T1.doZ0.W1  # stratum-specific causal effect in W = 1
exp.beta.W0     <- P.T1.doZ1.W0 / P.T1.doZ0.W0  # stratum-specific causal effect in W = 0

### Varying parameters ---------------------------------------------------------
BETA.T2         <- seq(from = 0, to = beta.T1, length.out = 4)
ETA.T2          <- c(0, 0.05, 0.1)

### Set replication parameters -------------------------------------------------
PART      <- 1:20   # splitting the replications
nreplic   <- 250    # total number of replications: Nreplic = nreplic * max(PART)
set.seed(1234)
seed      <- round(abs(rnorm(max(PART)) * 10^5))
param     <- as.data.frame(expand_grid(beta.T2 = BETA.T2, eta.T2 = ETA.T2, 
                                       part = PART))

### Function to replicate the example over the replications --------------------
onerun <- function(p) {
  beta.T2   <- param[p, 1]
  eta.T2    <- param[p, 2]
  part      <- param[p, 3]
  
  set.seed(seed[part])
  res <- NULL
  
  # Compute other counterfactual probabilities analytically --------------------
  P.T2.doZ1.W1    <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 1 + eta.T2 * 1 * 1) # P(T^{Z=1} = 2 | W = 1)
  P.T2.doZ0.W1    <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 1 + eta.T2 * 0 * 1) # P(T^{Z=0} = 2 | W = 1) 
  P.T2.doZ1.W0    <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 0 + eta.T2 * 1 * 0) # P(T^{Z=1} = 2 | W = 0)
  P.T2.doZ0.W0    <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 0 + eta.T2 * 0 * 0) # P(T^{Z=0} = 2 | W = 0)
  
  P.T2.doZ1       <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 1 + eta.T2 * 1 * 1) * pW + exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 0 + eta.T2 * 1 * 0) * (1 - pW) # P(T^{Z=1} = 2)
  P.T2.doZ0       <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 1 + eta.T2 * 0 * 1) * pW + exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 0 + eta.T2 * 0 * 0) * (1 - pW) # P(T^{Z=0} = 2)
  
  # Focus is on the stratum-specific contrasts
  W1.contrast         <- P.T1.doZ1.W1 / P.T2.doZ1.W1 # stratum-specific quantity estimated in practice from the SCCS data (analytic value)
  W0.contrast         <- P.T1.doZ1.W0 / P.T2.doZ1.W0 
  marginal.contrast   <- (P.T1.doZ1.W1 * p.Z.W1 * pW / 
                            (p.Z.W1 * pW + p.Z.W0 * (1 - pW)) + 
                            P.T1.doZ1.W0 * p.Z.W0 * (1 - pW) / 
                            (p.Z.W1 * pW + p.Z.W0 * (1 - pW))) / 
    (P.T2.doZ1.W1 * p.Z.W1 * pW / (p.Z.W1 * pW + p.Z.W0 * (1 - pW)) + 
       P.T2.doZ1.W0 * p.Z.W0 * (1 - pW) / (p.Z.W1 * pW + p.Z.W0 * (1 - pW))) # marginal quantity estimated in practice from the SCCS data (analytic value)
  
  for (nrep in 1:nreplic){
    
    # Simulate confounder W ----------------------------------------------------
    W         <- rbinom(n, 1, pW) # binary confounder
    n.W1      <- sum(W) 
    n.W0      <- sum(1 - W)
    
    # Simulate exposure Z ------------------------------------------------------
    Z         <- rep(NA, n) # binary exposure
    Z[W == 1] <- rbinom(n.W1, 1, p.Z.W1)
    Z[W == 0] <- rbinom(n.W0, 1, p.Z.W0)
    
    # Simulate the time to event outcome ---------------------------------------
    p1 <- exp(alpha.T1 + beta.T1 * Z + gamma.T1 * W + eta.T1 * Z * W)
    p2 <- exp(alpha.T2 + beta.T2 * Z + gamma.T2 * W + eta.T2 * Z * W)
    
    p1.Z0 <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * W + eta.T1 * 0 * W)
    p1.Z1 <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * W + eta.T1 * 1 * W)
    
    p2.Z0 <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * W + eta.T2 * 0 * W)
    p2.Z1 <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * W + eta.T2 * 1 * W)
    
    eT <- runif(n)
    T <- 2 * (eT <= p2) + 1 * ((eT > p2)&(eT <= (p1 + p2)))
    T.Z1 <- 2 * (eT <= p2.Z1) + 1 * ((eT > p2.Z1)&(eT <= (p1.Z1 + p2.Z1)))
    T.Z0 <- 2 * (eT <= p2.Z0) + 1 * ((eT > p2.Z0)&(eT <= (p1.Z0 + p2.Z0)))
    
    # SCCS data ----------------------------------------------------------------
    n.sccs          <- sum(T != 0 & Z == 1)
    n.sccs.W1       <- sum(T != 0 & Z == 1 & W == 1)
    n.sccs.W0       <- sum(T != 0 & Z == 1 & W == 0)
    n.sccs.T1.W1    <- sum(T == 1 & Z == 1 & W == 1)
    n.sccs.T1.W0    <- sum(T == 1 & Z == 1 & W == 0)
    
    while ((n.sccs.W0 == n.sccs.T1.W0)|
           (n.sccs.W1 == n.sccs.T1.W1)|
           (n.sccs.T1.W0 == 0)|
           (n.sccs.T1.W1 == 0)|
           (sum(T.Z0[W == 1] == 1) == 0)|
           (sum(T.Z1[W == 1] == 1) == 0)|
           (sum(T.Z0[W == 1] == 2) == 0)|
           (sum(T.Z1[W == 1] == 2) == 0)) {

      W         <- rbinom(n, 1, pW)
      n.W1      <- sum(W) 
      n.W0      <- sum(1 - W)
      
      Z         <- rep(NA, n) 
      Z[W == 1] <- rbinom(n.W1, 1, p.Z.W1)
      Z[W == 0] <- rbinom(n.W0, 1, p.Z.W0)
      
      p1 <- exp(alpha.T1 + beta.T1 * Z + gamma.T1 * W + eta.T1 * Z * W)
      p2 <- exp(alpha.T2 + beta.T2 * Z + gamma.T2 * W + eta.T2 * Z * W)
      
      p1.Z0 <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * W + eta.T1 * 0 * W)
      p1.Z1 <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * W + eta.T1 * 1 * W)
      
      p2.Z0 <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * W + eta.T2 * 0 * W)
      p2.Z1 <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * W + eta.T2 * 1 * W)
      
      eT <- runif(n)
      T <- 2 * (eT <= p2) + 1 * ((eT > p2)&(eT <= (p1 + p2)))
      T.Z1 <- 2 * (eT <= p2.Z1) + 1 * ((eT > p2.Z1)&(eT <= (p1.Z1 + p2.Z1)))
      T.Z0 <- 2 * (eT <= p2.Z0) + 1 * ((eT > p2.Z0)&(eT <= (p1.Z0 + p2.Z0)))
      
      n.sccs          <- sum(T != 0 & Z == 1)
      n.sccs.W1       <- sum(T != 0 & Z == 1 & W == 1)
      n.sccs.W0       <- sum(T != 0 & Z == 1 & W == 0)
      n.sccs.T1.W1    <- sum(T == 1 & Z == 1 & W == 1)
      n.sccs.T1.W0    <- sum(T == 1 & Z == 1 & W == 0)
    }
    
    pop             <- cbind(W, Z, p1, p2, p1.Z0, p1.Z1, T, T.Z1, T.Z0)       # general population of interest
    sccs            <- as.data.frame(pop[which((Z == 1)&(T != 0)),]) # SCCS population
    
    # Quantities estimated in practice, using the SCCS data --------------------
    W1.contrast.est       <- mean(T[W == 1 & Z == 1 & T != 0] == 1) / mean(T[W == 1 & Z == 1 & T != 0] == 2)
    # mean(T[W == 1 & Z == 1] == 1) / mean(T[W == 1 & Z == 1] == 2)
    
    W0.contrast.est       <- mean(T[W == 0 & Z == 1 & T != 0] == 1) / mean(T[W == 0 & Z == 1 & T != 0] == 2)
    # mean(T[W == 0 & Z == 1] == 1) / mean(T[W == 0 & Z == 1] == 2)
    
    marginal.contrast.est <- mean(sccs$T == 1) / mean(sccs$T == 2)
    
    res <- rbind(res, cbind(marginal.contrast.est = marginal.contrast.est,
                            W1.contrast.est = W1.contrast.est,
                            W0.contrast.est = W0.contrast.est,
                            n = n, n.W1 = n.W1, n.W0 = n.W0,
                            n.sccs = n.sccs, 
                            n.sccs.T1.W1 = n.sccs.T1.W1,
                            n.sccs.T1.W0 = n.sccs.T1.W0, 
                            n.sccs.W1 = n.sccs.W1,
                            n.sccs.W0 = n.sccs.W0,
                            pW = pW,
                            p.Z.W1 = p.Z.W1,
                            p.Z.W0 = p.Z.W0,
                            alpha.T1 = alpha.T1, 
                            alpha.T2 = alpha.T2,
                            gamma.T1 = gamma.T1,
                            gamma.T2 = gamma.T2,
                            eta.T1 = eta.T1,
                            eta.T2 = eta.T2,
                            beta.T1 = beta.T1,         
                            beta.T2 = beta.T2,
                            exp.beta = exp.beta,
                            exp.beta.W1 = exp.beta.W1,
                            exp.beta.W0 = exp.beta.W0,
                            marginal.contrast = marginal.contrast,
                            W0.contrast = W0.contrast,
                            W1.contrast = W1.contrast,
                            part = part))
  }
  myfile  <- paste0("simulres-example1-n", n, "-beta.T2", beta.T2, 
                    "-eta.T2", eta.T2, "-part", part, ".Rdata")
  save(res, file = myfile)
}

result <- mclapply(1:nrow(param), onerun, mc.cores = 64)

### Saving the simulation results ----------------------------------------------
RES <- NULL
for (p in 1:nrow(param)) {
  beta.T2   <- param[p, 1]
  eta.T2    <- param[p, 2]
  part      <- param[p, 3]
  
  load(paste0("simulres-example1-n", n, "-beta.T2", beta.T2, 
              "-eta.T2", eta.T2, "-part", part, ".Rdata"))
  RES <- rbind(RES, res)
}
RECAP           <- as.data.frame(RES)
ColNames        <- colnames(RECAP[,c(1:29)])
RECAP[ColNames] <- sapply(RECAP[ColNames], as.numeric)
RECAP$beta.T2   <- as.factor(RECAP$beta.T2)
RECAP$eta.T2    <- as.factor(RECAP$eta.T2)
myfile          <- paste0("SimulationResults-example1-n", n, ".Rdata")
save(RECAP, file = myfile)

### Details of the results  ----------------------------------------------------
load(paste0("SimulationResults-example1-n", n, ".Rdata"))

Nreplic   <- nreplic * max(PART)
param     <- param[,-3]
param     <- param[which(duplicated.matrix(param)==FALSE),]

details.contrast.W1    <- NULL
details.contrast.W0    <- NULL
details.contrast.marginal   <- NULL
for (i in 1:nrow(param)) {
  
  RECAP1 <- RECAP[((i-1) * (Nreplic) + 1):(i * (Nreplic)), ]
  
  # mean of estimated stratum-specific contrast for W = 1
  mean.est.log.W1.contrast  <- mean(log(RECAP1$W1.contrast.est))
  mean.est.W1.contrast      <- mean(RECAP1$W1.contrast.est)
  
  beta.W1     <- mean(log(RECAP$exp.beta.W1))
  relbias.W1  <- (mean.est.log.W1.contrast - beta.W1) / beta.W1
  
  # mean of estimated stratum-specific contrast for W = 0
  mean.est.log.W0.contrast  <- mean(log(RECAP1$W0.contrast.est))
  mean.est.W0.contrast      <- mean(RECAP1$W0.contrast.est)
  
  beta.W0     <- mean(log(RECAP$exp.beta.W0))
  relbias.W0  <- (mean.est.log.W0.contrast - beta.W0) / beta.W0
  
  # mean of estimated marginal contrast
  mean.est.log.marginal.contrast  <- mean(log(RECAP1$marginal.contrast.est))
  mean.est.marginal.contrast      <- mean(RECAP1$marginal.contrast.est)
  
  beta              <- mean(log(RECAP1$exp.beta))
  relbias.marginal  <- (mean.est.log.marginal.contrast - beta) / beta
  
  details.contrast.W1 <- rbind(details.contrast.W1, 
                                    c(mean.est.log.W1.contrast = mean.est.log.W1.contrast,
                                      mean.est.W1.contrast = mean.est.W1.contrast,
                                      relbias.W1 = relbias.W1,
                                      n = as.numeric(as.character(RECAP1$n[1])), 
                                      n.W1 = mean(RECAP1$n.W1),
                                      n.W0 = mean(RECAP1$n.W0), 
                                      n.sccs = mean(RECAP1$n.sccs), 
                                      n.sccs.W1 = mean(RECAP1$n.sccs.W1),
                                      n.sccs.W0 = mean(RECAP1$n.sccs.W0), 
                                      n.sccs.T1.W1 = mean(RECAP1$n.sccs.T1.W1),
                                      n.sccs.T1.W0 = mean(RECAP1$n.sccs.T1.W0), 
                                      pW = mean(RECAP1$pW), 
                                      p.Z.W1 = mean(RECAP1$p.Z.W1),
                                      p.Z.W0 = mean(RECAP1$p.Z.W0),
                                      beta.T1 = mean(RECAP1$beta.T1), 
                                      beta.T2 = as.numeric(as.character(RECAP1$beta.T2[1])),
                                      alpha.T1 = mean(RECAP1$alpha.T1), 
                                      alpha.T2 = mean(RECAP1$alpha.T2), 
                                      gamma.T1 = mean(RECAP1$gamma.T1), 
                                      gamma.T2 = mean(RECAP1$gamma.T2),
                                      eta.T1 = mean(RECAP1$eta.T1), 
                                      eta.T2 = as.numeric(as.character(RECAP1$eta.T2[1])),
                                      exp.beta = mean(RECAP1$exp.beta), 
                                      exp.beta.W1 = mean(RECAP1$exp.beta.W1), 
                                      exp.beta.W0 = mean(RECAP1$exp.beta.W0), 
                                      beta = beta,
                                      beta.W1 = beta.W1, 
                                      beta.W0 = beta.W0, 
                                      marginal.contrast = mean(RECAP1$marginal.contrast),
                                      W0.contrast = mean(RECAP1$W0.contrast),
                                      W1.contrast = mean(RECAP1$W1.contrast),
                                      log.marginal.contrast = mean(log(RECAP1$marginal.contrast)),
                                      log.W1.contrast = mean(log(RECAP1$W1.contrast)), 
                                      log.W0.contrast = mean(log(RECAP1$W0.contrast))))
  
  details.contrast.W0 <- rbind(details.contrast.W0, 
                               c(mean.est.log.W0.contrast = mean.est.log.W0.contrast,
                                 mean.est.W0.contrast = mean.est.W0.contrast,
                                 relbias.W0 = relbias.W0,
                                 n = as.numeric(as.character(RECAP1$n[1])), 
                                 n.W1 = mean(RECAP1$n.W1),
                                 n.W0 = mean(RECAP1$n.W0), 
                                 n.sccs = mean(RECAP1$n.sccs), 
                                 n.sccs.W1 = mean(RECAP1$n.sccs.W1),
                                 n.sccs.W0 = mean(RECAP1$n.sccs.W0), 
                                 n.sccs.T1.W1 = mean(RECAP1$n.sccs.T1.W1),
                                 n.sccs.T1.W0 = mean(RECAP1$n.sccs.T1.W0), 
                                 pW = mean(RECAP1$pW), 
                                 p.Z.W1 = mean(RECAP1$p.Z.W1),
                                 p.Z.W0 = mean(RECAP1$p.Z.W0),
                                 beta.T1 = mean(RECAP1$beta.T1), 
                                 beta.T2 = as.numeric(as.character(RECAP1$beta.T2[1])),
                                 alpha.T1 = mean(RECAP1$alpha.T1), 
                                 alpha.T2 = mean(RECAP1$alpha.T2), 
                                 gamma.T1 = mean(RECAP1$gamma.T1), 
                                 gamma.T2 = mean(RECAP1$gamma.T2),
                                 eta.T1 = mean(RECAP1$eta.T1), 
                                 eta.T2 = as.numeric(as.character(RECAP1$eta.T2[1])),
                                 exp.beta = mean(RECAP1$exp.beta), 
                                 exp.beta.W1 = mean(RECAP1$exp.beta.W1), 
                                 exp.beta.W0 = mean(RECAP1$exp.beta.W0), 
                                 beta = log(mean(RECAP1$exp.beta)),
                                 beta.W1 = log(mean(RECAP1$exp.beta.W1)), 
                                 beta.W0 = log(mean(RECAP1$exp.beta.W0)), 
                                 marginal.contrast = mean(RECAP1$marginal.contrast),
                                 W0.contrast = mean(RECAP1$W0.contrast),
                                 W1.contrast = mean(RECAP1$W1.contrast),
                                 log.marginal.contrast = mean(log(RECAP1$marginal.contrast)),
                                 log.W1.contrast = mean(log(RECAP1$W1.contrast)), 
                                 log.W0.contrast = mean(log(RECAP1$W0.contrast))))
  
  details.contrast.marginal <- rbind(details.contrast.marginal, 
                                     c(mean.est.log.marginal.contrast = mean.est.log.marginal.contrast,
                                       mean.est.marginal.contrast = mean.est.marginal.contrast,
                                       relbias.marginal = relbias.marginal,
                                       n = as.numeric(as.character(RECAP1$n[1])), 
                                       n.W1 = mean(RECAP1$n.W1),
                                       n.W0 = mean(RECAP1$n.W0), 
                                       n.sccs = mean(RECAP1$n.sccs), 
                                       n.sccs.W1 = mean(RECAP1$n.sccs.W1),
                                       n.sccs.W0 = mean(RECAP1$n.sccs.W0), 
                                       n.sccs.T1.W1 = mean(RECAP1$n.sccs.T1.W1),
                                       n.sccs.T1.W0 = mean(RECAP1$n.sccs.T1.W0), 
                                       pW = mean(RECAP1$pW), 
                                       p.Z.W1 = mean(RECAP1$p.Z.W1),
                                       p.Z.W0 = mean(RECAP1$p.Z.W0),
                                       beta.T1 = mean(RECAP1$beta.T1), 
                                       beta.T2 = as.numeric(as.character(RECAP1$beta.T2[1])),
                                       alpha.T1 = mean(RECAP1$alpha.T1), 
                                       alpha.T2 = mean(RECAP1$alpha.T2), 
                                       gamma.T1 = mean(RECAP1$gamma.T1), 
                                       gamma.T2 = mean(RECAP1$gamma.T2),
                                       eta.T1 = mean(RECAP1$eta.T1), 
                                       eta.T2 = as.numeric(as.character(RECAP1$eta.T2[1])),
                                       exp.beta = mean(RECAP1$exp.beta), 
                                       exp.beta.W1 = mean(RECAP1$exp.beta.W1), 
                                       exp.beta.W0 = mean(RECAP1$exp.beta.W0), 
                                       beta = log(mean(RECAP1$exp.beta)),
                                       beta.W1 = log(mean(RECAP1$exp.beta.W1)), 
                                       beta.W0 = log(mean(RECAP1$exp.beta.W0)), 
                                       marginal.contrast = mean(RECAP1$marginal.contrast),
                                       W0.contrast = mean(RECAP1$W0.contrast),
                                       W1.contrast = mean(RECAP1$W1.contrast),
                                       log.marginal.contrast = mean(log(RECAP1$marginal.contrast)),
                                       log.W1.contrast = mean(log(RECAP1$W1.contrast)), 
                                       log.W0.contrast = mean(log(RECAP1$W0.contrast))))
}

details.contrast.W1 <- as.data.frame(details.contrast.W1)
details.contrast.W0 <- as.data.frame(details.contrast.W0)
details.contrast.marginal <- as.data.frame(details.contrast.marginal)

details.contrast.W1$exp.beta.T2 <- exp(details.contrast.W1$beta.T2)
details.contrast.W0$exp.beta.T2 <- exp(details.contrast.W0$beta.T2)

### Plot of the results --------------------------------------------------------
details.contrast.W1$eta.T2 <- as.factor(details.contrast.W1$eta.T2)
details.contrast.W1$eta.T2 <- factor(details.contrast.W1$eta.T2, 
                                     labels = c("eta[2]==0", "eta[2]==0.05", 
                                                "eta[2]==0.1"))

details.contrast.W0$eta.T2 <- as.factor(details.contrast.W0$eta.T2)
details.contrast.W0$eta.T2 <- factor(details.contrast.W0$eta.T2, 
                                     labels = c("eta[2]==0", "eta[2]==0.05", 
                                                "eta[2]==0.1"))

results <- rbind(cbind(value = details.contrast.W1$mean.est.log.W1.contrast, 
                       beta.T2 = details.contrast.W1$beta.T2, 
                       eta.T2 = details.contrast.W1$eta.T2, 
                       Contrast = "cond.contrast.W1"),
                 cbind(value = details.contrast.W0$mean.est.log.W0.contrast, 
                       beta.T2 = details.contrast.W0$beta.T2, 
                       eta.T2 = details.contrast.W0$eta.T2, 
                       Contrast = "cond.contrast.W0"))
results           <- as.data.frame(results)
results$value     <- as.numeric(results$value)
results$beta.T2   <- as.numeric(results$beta.T2)
results$eta.T2    <- as.factor(results$eta.T2)
results$eta.T2    <- factor(results$eta.T2, labels = c("eta[2]==0", 
                                                       "eta[2]==0.05", 
                                                       "eta[2]==0.1"))
results$Contrast  <- as.factor(results$Contrast)

pdf(file = paste0("cond.contrasts_example1.pdf"), width = 8, height = 6)
ggplot(data = results, aes(x = beta.T2, y = value)) + 
  geom_line(linewidth = 0.8, aes(linetype = Contrast))  + 
  scale_x_continuous(expression(beta[2]), 
                     breaks = c(0, 0.2, 0.4, 0.6), 
                     labels = c("0", "0.2", "0.4", expression(beta[1]))) +
  scale_y_continuous("mean of logarithm of estimated stratum-specific contrast", 
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), 
                     labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5", 
                                expression(beta[W == 0]), "0.7", 
                                expression(beta[W == 1]))) +
  #scale_y_continuous("mean of estimated log-contrast") +
  facet_wrap(~ eta.T2, labeller = label_parsed) + theme_bw(base_size = 13) +
  scale_linetype_manual(name = "Stratum",
                        values = c(cond.contrast.W1 = "dashed", 
                                   cond.contrast.W0 = "dotted"),
                        labels = c(cond.contrast.W1 = "W = 1", 
                                   cond.contrast.W0 = "W = 0")) 
dev.off()

### Trying to remove the small sample bias by considering population with 
### n = 10^7 individuals instead -----------------------------------------------

n         <- 10^7   
PART      <- 1:20   
nreplic   <- 250    
set.seed(1234)
seed      <- round(abs(rnorm(max(PART)) * 10^5))
param     <- as.data.frame(expand_grid(beta.T2 = BETA.T2, eta.T2 = ETA.T2, 
                                       part = PART))
onerun    <- function(p) {
  beta.T2   <- param[p, 1]
  eta.T2    <- param[p, 2]
  part      <- param[p, 3]
  
  set.seed(seed[part])
  res <- NULL
  
  # Compute other counterfactual probabilities analytically --------------------
  P.T2.doZ1.W1    <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 1 + eta.T2 * 1 * 1) # P(T^{Z=1} = 2 | W = 1)
  P.T2.doZ0.W1    <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 1 + eta.T2 * 0 * 1) # P(T^{Z=0} = 2 | W = 1) 
  P.T2.doZ1.W0    <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 0 + eta.T2 * 1 * 0) # P(T^{Z=1} = 2 | W = 0)
  P.T2.doZ0.W0    <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 0 + eta.T2 * 0 * 0) # P(T^{Z=0} = 2 | W = 0)
  
  P.T2.doZ1       <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 1 + eta.T2 * 1 * 1) * pW + exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 0 + eta.T2 * 1 * 0) * (1 - pW) # P(T^{Z=1} = 2)
  P.T2.doZ0       <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 1 + eta.T2 * 0 * 1) * pW + exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 0 + eta.T2 * 0 * 0) * (1 - pW) # P(T^{Z=0} = 2)
  
  # Focus is on the stratum-specific contrasts
  W1.contrast         <- P.T1.doZ1.W1 / P.T2.doZ1.W1 # stratum-specific quantity estimated in practice from the SCCS data (analytic value)
  W0.contrast         <- P.T1.doZ1.W0 / P.T2.doZ1.W0 
  marginal.contrast   <- (P.T1.doZ1.W1 * p.Z.W1 * pW / 
                            (p.Z.W1 * pW + p.Z.W0 * (1 - pW)) + 
                            P.T1.doZ1.W0 * p.Z.W0 * (1 - pW) / 
                            (p.Z.W1 * pW + p.Z.W0 * (1 - pW))) / 
    (P.T2.doZ1.W1 * p.Z.W1 * pW / (p.Z.W1 * pW + p.Z.W0 * (1 - pW)) + 
       P.T2.doZ1.W0 * p.Z.W0 * (1 - pW) / (p.Z.W1 * pW + p.Z.W0 * (1 - pW))) # marginal quantity estimated in practice from the SCCS data (analytic value)
  
  for (nrep in 1:nreplic){
    
    # Simulate confounder W ----------------------------------------------------
    W         <- rbinom(n, 1, pW) # binary confounder
    n.W1      <- sum(W) 
    n.W0      <- sum(1 - W)
    
    # Simulate exposure Z ------------------------------------------------------
    Z         <- rep(NA, n) # binary exposure
    Z[W == 1] <- rbinom(n.W1, 1, p.Z.W1)
    Z[W == 0] <- rbinom(n.W0, 1, p.Z.W0)
    
    # Simulate the time to event outcome ---------------------------------------
    p1 <- exp(alpha.T1 + beta.T1 * Z + gamma.T1 * W + eta.T1 * Z * W)
    p2 <- exp(alpha.T2 + beta.T2 * Z + gamma.T2 * W + eta.T2 * Z * W)
    
    p1.Z0 <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * W + eta.T1 * 0 * W)
    p1.Z1 <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * W + eta.T1 * 1 * W)
    
    p2.Z0 <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * W + eta.T2 * 0 * W)
    p2.Z1 <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * W + eta.T2 * 1 * W)
    
    eT <- runif(n)
    T <- 2 * (eT <= p2) + 1 * ((eT > p2)&(eT <= (p1 + p2)))
    T.Z1 <- 2 * (eT <= p2.Z1) + 1 * ((eT > p2.Z1)&(eT <= (p1.Z1 + p2.Z1)))
    T.Z0 <- 2 * (eT <= p2.Z0) + 1 * ((eT > p2.Z0)&(eT <= (p1.Z0 + p2.Z0)))
    
    # SCCS data ----------------------------------------------------------------
    n.sccs          <- sum(T != 0 & Z == 1)
    n.sccs.W1       <- sum(T != 0 & Z == 1 & W == 1)
    n.sccs.W0       <- sum(T != 0 & Z == 1 & W == 0)
    n.sccs.T1.W1    <- sum(T == 1 & Z == 1 & W == 1)
    n.sccs.T1.W0    <- sum(T == 1 & Z == 1 & W == 0)
    
    while ((n.sccs.W0 == n.sccs.T1.W0)|
           (n.sccs.W1 == n.sccs.T1.W1)|
           (n.sccs.T1.W0 == 0)|
           (n.sccs.T1.W1 == 0)|
           (sum(T.Z0[W == 1] == 1) == 0)|
           (sum(T.Z1[W == 1] == 1) == 0)|
           (sum(T.Z0[W == 1] == 2) == 0)|
           (sum(T.Z1[W == 1] == 2) == 0)) {
      
      W         <- rbinom(n, 1, pW)
      n.W1      <- sum(W) 
      n.W0      <- sum(1 - W)
      
      Z         <- rep(NA, n) 
      Z[W == 1] <- rbinom(n.W1, 1, p.Z.W1)
      Z[W == 0] <- rbinom(n.W0, 1, p.Z.W0)
      
      p1 <- exp(alpha.T1 + beta.T1 * Z + gamma.T1 * W + eta.T1 * Z * W)
      p2 <- exp(alpha.T2 + beta.T2 * Z + gamma.T2 * W + eta.T2 * Z * W)
      
      p1.Z0 <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * W + eta.T1 * 0 * W)
      p1.Z1 <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * W + eta.T1 * 1 * W)
      
      p2.Z0 <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * W + eta.T2 * 0 * W)
      p2.Z1 <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * W + eta.T2 * 1 * W)
      
      eT <- runif(n)
      T <- 2 * (eT <= p2) + 1 * ((eT > p2)&(eT <= (p1 + p2)))
      T.Z1 <- 2 * (eT <= p2.Z1) + 1 * ((eT > p2.Z1)&(eT <= (p1.Z1 + p2.Z1)))
      T.Z0 <- 2 * (eT <= p2.Z0) + 1 * ((eT > p2.Z0)&(eT <= (p1.Z0 + p2.Z0)))
      
      n.sccs          <- sum(T != 0 & Z == 1)
      n.sccs.W1       <- sum(T != 0 & Z == 1 & W == 1)
      n.sccs.W0       <- sum(T != 0 & Z == 1 & W == 0)
      n.sccs.T1.W1    <- sum(T == 1 & Z == 1 & W == 1)
      n.sccs.T1.W0    <- sum(T == 1 & Z == 1 & W == 0)
    }
    
    pop             <- cbind(W, Z, p1, p2, p1.Z0, p1.Z1, T, T.Z1, T.Z0)       # general population of interest
    sccs            <- as.data.frame(pop[which((Z == 1)&(T != 0)),]) # SCCS population
    
    # Quantity estimated in practice, using the SCCS data ----------------------
    
    W1.contrast.est       <- mean(T[W == 1 & Z == 1 & T != 0] == 1) / mean(T[W == 1 & Z == 1 & T != 0] == 2)
    # mean(T[W == 1 & Z == 1] == 1) / mean(T[W == 1 & Z == 1] == 2)
    
    W0.contrast.est       <- mean(T[W == 0 & Z == 1 & T != 0] == 1) / mean(T[W == 0 & Z == 1 & T != 0] == 2)
    # mean(T[W == 0 & Z == 1] == 1) / mean(T[W == 0 & Z == 1] == 2)
    
    marginal.contrast.est <- mean(T[Z == 1 & T != 0] == 1) / mean(T[Z == 1 & T != 0] == 2)
    
    res <- rbind(res, cbind(marginal.contrast.est = marginal.contrast.est,
                            W1.contrast.est = W1.contrast.est,
                            W0.contrast.est = W0.contrast.est,
                            n = n, n.W1 = n.W1, n.W0 = n.W0,
                            n.sccs = n.sccs, 
                            n.sccs.T1.W1 = n.sccs.T1.W1,
                            n.sccs.T1.W0 = n.sccs.T1.W0, 
                            n.sccs.W1 = n.sccs.W1,
                            n.sccs.W0 = n.sccs.W0,
                            pW = pW,
                            alpha.T1 = alpha.T1, 
                            alpha.T2 = alpha.T2,
                            gamma.T1 = gamma.T1,
                            gamma.T2 = gamma.T2,
                            eta.T1 = eta.T1,
                            eta.T2 = eta.T2,
                            beta.T1 = beta.T1,         
                            beta.T2 = beta.T2,
                            p.Z.W1 = p.Z.W1,
                            p.Z.W0 = p.Z.W0,
                            exp.beta = exp.beta,
                            exp.beta.W1 = exp.beta.W1,
                            exp.beta.W0 = exp.beta.W0,
                            marginal.contrast = marginal.contrast,
                            W0.contrast = W0.contrast,
                            W1.contrast = W1.contrast,
                            part = part))
  }
  myfile  <- paste0("simulres-example1-n", n, "-beta.T2", beta.T2, 
                    "-eta.T2", eta.T2, "-part", part, ".Rdata")
  save(res, file = myfile)
}
param     <- param[1:20,] # only running the scenario when conditions (i) and (ii) hold
resultat  <- mclapply(1:nrow(param), onerun, mc.cores = 20)

RES <- NULL
for (p in 1:nrow(param)) {
  beta.T2   <- param[p, 1]
  eta.T2    <- param[p, 2]
  part      <- param[p, 3]
  
  load(paste0("simulres-example1-n", n, "-beta.T2", beta.T2, 
              "-eta.T2", eta.T2, "-part", part, ".Rdata"))
  RES <- rbind(RES, res)
}
RECAP           <- as.data.frame(RES)
ColNames        <- colnames(RECAP[,c(1:29)])
RECAP[ColNames] <- sapply(RECAP[ColNames], as.numeric)
myfile  <- paste0("SimulationResults-example1-n", n, ".Rdata")
save(RECAP, file = myfile)

n <- 10^7
load(paste0("SimulationResults-example1-n", n, ".Rdata"))
  
# mean of estimated stratum-specific contrast for W = 1
mean.est.log.W1.contrast  <- mean(log(RECAP$W1.contrast.est))
mean.est.W1.contrast      <- mean(RECAP$W1.contrast.est)

beta.W1     <- mean(log(RECAP$exp.beta.W1))
relbias.W1  <- (mean.est.log.W1.contrast - beta.W1) / beta.W1

# mean of estimated stratum-specific contrast for W = 0
mean.est.log.W0.contrast  <- mean(log(RECAP$W0.contrast.est))
mean.est.W0.contrast      <- mean(RECAP$W0.contrast.est)

beta.W0     <- mean(log(RECAP$exp.beta.W0))
relbias.W0  <- (mean.est.log.W0.contrast - beta.W0) / beta.W0

# mean of estimated marginal contrast
mean.est.log.marginal.contrast  <- mean(log(RECAP$marginal.contrast.est))
mean.est.marginal.contrast      <- mean(RECAP$marginal.contrast.est)

beta              <- mean(log(RECAP$exp.beta))
relbias.marginal  <- (mean.est.log.marginal.contrast - beta) / beta

### Trying to remove the small sample bias by considering population with 
### n = 5 * 10^6 individuals instead -----------------------------------------------

n         <- 5 * 10^6
PART      <- 1:20   
nreplic   <- 250    
set.seed(1234)
seed      <- round(abs(rnorm(max(PART)) * 10^5))
param     <- as.data.frame(expand_grid(beta.T2 = BETA.T2, eta.T2 = ETA.T2, 
                                       part = PART))
onerun    <- function(p) {
  beta.T2   <- param[p, 1]
  eta.T2    <- param[p, 2]
  part      <- param[p, 3]
  
  set.seed(seed[part])
  res <- NULL
  
  # Compute other counterfactual probabilities analytically --------------------
  P.T2.doZ1.W1    <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 1 + eta.T2 * 1 * 1) # P(T^{Z=1} = 2 | W = 1)
  P.T2.doZ0.W1    <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 1 + eta.T2 * 0 * 1) # P(T^{Z=0} = 2 | W = 1) 
  P.T2.doZ1.W0    <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 0 + eta.T2 * 1 * 0) # P(T^{Z=1} = 2 | W = 0)
  P.T2.doZ0.W0    <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 0 + eta.T2 * 0 * 0) # P(T^{Z=0} = 2 | W = 0)
  
  P.T2.doZ1       <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 1 + eta.T2 * 1 * 1) * pW + exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 0 + eta.T2 * 1 * 0) * (1 - pW) # P(T^{Z=1} = 2)
  P.T2.doZ0       <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 1 + eta.T2 * 0 * 1) * pW + exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 0 + eta.T2 * 0 * 0) * (1 - pW) # P(T^{Z=0} = 2)
  
  # Focus is on the stratum-specific contrasts
  W1.contrast         <- P.T1.doZ1.W1 / P.T2.doZ1.W1 # stratum-specific quantity estimated in practice from the SCCS data (analytic value)
  W0.contrast         <- P.T1.doZ1.W0 / P.T2.doZ1.W0 
  marginal.contrast   <- (P.T1.doZ1.W1 * p.Z.W1 * pW / 
                            (p.Z.W1 * pW + p.Z.W0 * (1 - pW)) + 
                            P.T1.doZ1.W0 * p.Z.W0 * (1 - pW) / 
                            (p.Z.W1 * pW + p.Z.W0 * (1 - pW))) / 
    (P.T2.doZ1.W1 * p.Z.W1 * pW / (p.Z.W1 * pW + p.Z.W0 * (1 - pW)) + 
       P.T2.doZ1.W0 * p.Z.W0 * (1 - pW) / (p.Z.W1 * pW + p.Z.W0 * (1 - pW))) # marginal quantity estimated in practice from the SCCS data (analytic value)
  
  for (nrep in 1:nreplic){
    
    # Simulate confounder W ----------------------------------------------------
    W         <- rbinom(n, 1, pW) # binary confounder
    n.W1      <- sum(W) 
    n.W0      <- sum(1 - W)
    
    # Simulate exposure Z ------------------------------------------------------
    Z         <- rep(NA, n) # binary exposure
    Z[W == 1] <- rbinom(n.W1, 1, p.Z.W1)
    Z[W == 0] <- rbinom(n.W0, 1, p.Z.W0)
    
    # Simulate the time to event outcome ---------------------------------------
    p1 <- exp(alpha.T1 + beta.T1 * Z + gamma.T1 * W + eta.T1 * Z * W)
    p2 <- exp(alpha.T2 + beta.T2 * Z + gamma.T2 * W + eta.T2 * Z * W)
    
    p1.Z0 <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * W + eta.T1 * 0 * W)
    p1.Z1 <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * W + eta.T1 * 1 * W)
    
    p2.Z0 <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * W + eta.T2 * 0 * W)
    p2.Z1 <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * W + eta.T2 * 1 * W)
    
    eT <- runif(n)
    T <- 2 * (eT <= p2) + 1 * ((eT > p2)&(eT <= (p1 + p2)))
    T.Z1 <- 2 * (eT <= p2.Z1) + 1 * ((eT > p2.Z1)&(eT <= (p1.Z1 + p2.Z1)))
    T.Z0 <- 2 * (eT <= p2.Z0) + 1 * ((eT > p2.Z0)&(eT <= (p1.Z0 + p2.Z0)))
    
    # SCCS data ----------------------------------------------------------------
    n.sccs          <- sum(T != 0 & Z == 1)
    n.sccs.W1       <- sum(T != 0 & Z == 1 & W == 1)
    n.sccs.W0       <- sum(T != 0 & Z == 1 & W == 0)
    n.sccs.T1.W1    <- sum(T == 1 & Z == 1 & W == 1)
    n.sccs.T1.W0    <- sum(T == 1 & Z == 1 & W == 0)
    
    while ((n.sccs.W0 == n.sccs.T1.W0)|
           (n.sccs.W1 == n.sccs.T1.W1)|
           (n.sccs.T1.W0 == 0)|
           (n.sccs.T1.W1 == 0)|
           (sum(T.Z0[W == 1] == 1) == 0)|
           (sum(T.Z1[W == 1] == 1) == 0)|
           (sum(T.Z0[W == 1] == 2) == 0)|
           (sum(T.Z1[W == 1] == 2) == 0)) {
      
      W         <- rbinom(n, 1, pW)
      n.W1      <- sum(W) 
      n.W0      <- sum(1 - W)
      
      Z         <- rep(NA, n) 
      Z[W == 1] <- rbinom(n.W1, 1, p.Z.W1)
      Z[W == 0] <- rbinom(n.W0, 1, p.Z.W0)
      
      p1 <- exp(alpha.T1 + beta.T1 * Z + gamma.T1 * W + eta.T1 * Z * W)
      p2 <- exp(alpha.T2 + beta.T2 * Z + gamma.T2 * W + eta.T2 * Z * W)
      
      p1.Z0 <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * W + eta.T1 * 0 * W)
      p1.Z1 <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * W + eta.T1 * 1 * W)
      
      p2.Z0 <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * W + eta.T2 * 0 * W)
      p2.Z1 <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * W + eta.T2 * 1 * W)
      
      eT <- runif(n)
      T <- 2 * (eT <= p2) + 1 * ((eT > p2)&(eT <= (p1 + p2)))
      T.Z1 <- 2 * (eT <= p2.Z1) + 1 * ((eT > p2.Z1)&(eT <= (p1.Z1 + p2.Z1)))
      T.Z0 <- 2 * (eT <= p2.Z0) + 1 * ((eT > p2.Z0)&(eT <= (p1.Z0 + p2.Z0)))
      
      n.sccs          <- sum(T != 0 & Z == 1)
      n.sccs.W1       <- sum(T != 0 & Z == 1 & W == 1)
      n.sccs.W0       <- sum(T != 0 & Z == 1 & W == 0)
      n.sccs.T1.W1    <- sum(T == 1 & Z == 1 & W == 1)
      n.sccs.T1.W0    <- sum(T == 1 & Z == 1 & W == 0)
    }
    
    pop             <- cbind(W, Z, p1, p2, p1.Z0, p1.Z1, T, T.Z1, T.Z0)       # general population of interest
    sccs            <- as.data.frame(pop[which((Z == 1)&(T != 0)),]) # SCCS population
    
    # Quantity estimated in practice, using the SCCS data ----------------------
    
    W1.contrast.est       <- mean(T[W == 1 & Z == 1 & T != 0] == 1) / mean(T[W == 1 & Z == 1 & T != 0] == 2)
    # mean(T[W == 1 & Z == 1] == 1) / mean(T[W == 1 & Z == 1] == 2)
    
    W0.contrast.est       <- mean(T[W == 0 & Z == 1 & T != 0] == 1) / mean(T[W == 0 & Z == 1 & T != 0] == 2)
    # mean(T[W == 0 & Z == 1] == 1) / mean(T[W == 0 & Z == 1] == 2)
    
    marginal.contrast.est <- mean(T[Z == 1 & T != 0] == 1) / mean(T[Z == 1 & T != 0] == 2)
    
    res <- rbind(res, cbind(marginal.contrast.est = marginal.contrast.est,
                            W1.contrast.est = W1.contrast.est,
                            W0.contrast.est = W0.contrast.est,
                            n = n, n.W1 = n.W1, n.W0 = n.W0,
                            n.sccs = n.sccs, 
                            n.sccs.T1.W1 = n.sccs.T1.W1,
                            n.sccs.T1.W0 = n.sccs.T1.W0, 
                            n.sccs.W1 = n.sccs.W1,
                            n.sccs.W0 = n.sccs.W0,
                            pW = pW,
                            alpha.T1 = alpha.T1, 
                            alpha.T2 = alpha.T2,
                            gamma.T1 = gamma.T1,
                            gamma.T2 = gamma.T2,
                            eta.T1 = eta.T1,
                            eta.T2 = eta.T2,
                            beta.T1 = beta.T1,         
                            beta.T2 = beta.T2,
                            p.Z.W1 = p.Z.W1,
                            p.Z.W0 = p.Z.W0,
                            exp.beta = exp.beta,
                            exp.beta.W1 = exp.beta.W1,
                            exp.beta.W0 = exp.beta.W0,
                            marginal.contrast = marginal.contrast,
                            W0.contrast = W0.contrast,
                            W1.contrast = W1.contrast,
                            part = part))
  }
  myfile  <- paste0("simulres-example1-n", n, "-beta.T2", beta.T2, 
                    "-eta.T2", eta.T2, "-part", part, ".Rdata")
  save(res, file = myfile)
}
param     <- param[1:20,] # only running the scenario when conditions (i) and (ii) hold
resultat  <- mclapply(1:nrow(param), onerun, mc.cores = 20)

RES <- NULL
for (p in 1:nrow(param)) {
  beta.T2   <- param[p, 1]
  eta.T2    <- param[p, 2]
  part      <- param[p, 3]
  
  load(paste0("simulres-example1-n", n, "-beta.T2", beta.T2, 
              "-eta.T2", eta.T2, "-part", part, ".Rdata"))
  RES <- rbind(RES, res)
}
RECAP           <- as.data.frame(RES)
ColNames        <- colnames(RECAP[,c(1:29)])
RECAP[ColNames] <- sapply(RECAP[ColNames], as.numeric)
myfile  <- paste0("SimulationResults-example1-n", n, ".Rdata")
save(RECAP, file = myfile)

5 * 10^6
load(paste0("SimulationResults-example1-n", n, ".Rdata"))

# mean of estimated stratum-specific contrast for W = 1
mean.est.log.W1.contrast  <- mean(log(RECAP$W1.contrast.est))
mean.est.W1.contrast      <- mean(RECAP$W1.contrast.est)
beta.W1     <- mean(log(RECAP$exp.beta.W1))
relbias.W1  <- (mean.est.log.W1.contrast - beta.W1) / beta.W1

# mean of estimated stratum-specific contrast for W = 0
mean.est.log.W0.contrast  <- mean(log(RECAP$W0.contrast.est))
mean.est.W0.contrast      <- mean(RECAP$W0.contrast.est)
beta.W0     <- mean(log(RECAP$exp.beta.W0))
relbias.W0  <- (mean.est.log.W0.contrast - beta.W0) / beta.W0

# mean of estimated marginal contrast
mean.est.log.marginal.contrast  <- mean(log(RECAP$marginal.contrast.est))
mean.est.marginal.contrast      <- mean(RECAP$marginal.contrast.est)
beta              <- mean(log(RECAP$exp.beta))
relbias.marginal  <- (mean.est.log.marginal.contrast - beta) / beta
