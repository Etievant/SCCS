## -----------------------------------------------------------------------------
## Description: This script illustrates Detection bias as presented in Section 8
##              in Etievant, Gail and Follmann (2024)
##
##              The simulation is in Web Appendix D.1 in the Supporting 
##              Information Document
##
##              We consider a scenario where the effect of the vaccine is null 
##              and homogeneous across strata of the confounder
## -----------------------------------------------------------------------------

### load packages --------------------------------------------------------------
library(tidyr)
library(parallel)
library(ggplot2)

### Set fixed parameters -------------------------------------------------------
n             <- 10^6             # (overall) population size
pW            <- 0.6              # P(W = 1)
alpha.T1      <- log(0.00008)     # log baseline risk in period 1 
alpha.T2      <- log(0.00008)     # log baseline risk in period 2 
gamma.T1      <- 0.2              # we use g_1: w -> gamma.T1 * w
gamma.T2      <- 0.2              # we use g_2: w -> gamma.T2 * w
eta.T1        <- 0                # we use h_1: w -> eta.T1 * w
eta.T2        <- 0                # we use h_2: w -> eta.T2 * w
beta.T1       <- 0                # null effect of Z
beta.T2       <- 0   
p.Z.W1        <- 0.2              # P(Z = 1 | W = 1)
p.Z.W0        <- 0.35             # P(Z = 1 | W = 0)
p.detect.T1   <- 0.98             # P(D = 1 | T = 1)

### Compute counterfactual probabilities analytically --------------------------
P.T1.doZ0.W1    <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * 1 + eta.T1 * 0 * 1) # P(T^{Z=0} = 1 | W = 1)
P.T1.doZ0.W0    <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * 0 + eta.T1 * 0 * 0) # P(T^{Z=0} = 1 | W = 0)
P.T1.doZ1.W1    <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 1 + eta.T1 * 1 * 1) # P(T^{Z=1} = 1 | W = 1)
P.T1.doZ1.W0    <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 0 + eta.T1 * 1 * 0) # P(T^{Z=1} = 1 | W = 0)

P.T1.doZ1       <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 1 + eta.T1 * 1 * 1) * pW + exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 0 + eta.T1 * 1 * 0) * (1 - pW) # P(T^{Z=1} = 1)
P.T1.doZ0       <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * 1 + eta.T1 * 0 * 1) * pW + exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * 0 + eta.T1 * 0 * 0) * (1 - pW) # P(T^{Z=0} = 1)

exp.beta        <- P.T1.doZ1 / P.T1.doZ0        # marginal causal risk ratio

### Varying parameter ----------------------------------------------------------
P.DETEC.T2 <- seq(from = p.detect.T1, to = 0.7, length.out = 15) # P(D = 1 | T = 2)

### Set replication parameters -------------------------------------------------
PART      <- 1:20   # splitting the replications
nreplic   <- 250    # total number of replications: Nreplic = nreplic * max(PART) 
set.seed(12345)
seed      <- round(abs(rnorm(max(PART)) * 10^5))
param     <- as.data.frame(expand_grid(p.detect.T2 = P.DETEC.T2, part = PART))

### Function to replicate the example over the replications --------------------
onerun <- function(p) {
  p.detect.T2   <- param[p, 1]
  part          <- param[p, 2]
  
  set.seed(seed[part])
  res <- NULL
  
  # Compute other counterfactual probabilities analytically --------------------
  P.T2.doZ1.W1    <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 1 + eta.T2 * 1 * 1)  # P(T^{Z=1} = 2 | W = 1)
  P.T2.doZ0.W1    <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 1 + eta.T2 * 0 * 1)  # P(T^{Z=0} = 2 | W = 1) 
  P.T2.doZ1.W0    <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * 0 + eta.T2 * 1 * 0)  # P(T^{Z=1} = 2 | W = 0)
  P.T2.doZ0.W0    <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * 0 + eta.T2 * 0 * 0)  # P(T^{Z=0} = 2 | W = 0)

  # Focus is on the marginal contrast as the effect is homogeneous across strata
  marginal.contrast   <- (P.T1.doZ1.W1 * p.Z.W1 * pW / 
                            (p.Z.W1 * pW + p.Z.W0 * (1 - pW)) + 
                            P.T1.doZ1.W0 * p.Z.W0 * (1 - pW) / 
                            (p.Z.W1 * pW + p.Z.W0 * (1 - pW))) * p.detect.T1 / 
    (P.T2.doZ1.W1 * p.Z.W1 * pW / (p.Z.W1 * pW + p.Z.W0 * (1 - pW)) + 
       P.T2.doZ1.W0 * p.Z.W0 * (1 - pW) / (p.Z.W1 * pW + p.Z.W0 * (1 - pW))) / 
    p.detect.T2 # marginal quantity estimated in practice from the SCCS data (analytic value)
  
  for (nrep in 1:nreplic){
    
    # Simulate confounder W ----------------------------------------------------
    W         <- rbinom(n = n, size = 1, prob = pW) # binary confounder
    n.W1      <- sum(W) 
    n.W0      <- sum(1 - W)
    
    # Simulate exposure Z ------------------------------------------------------
    Z.W1      <- sample(c(0,1), n.W1, replace = TRUE, 
                        prob = c(1 - p.Z.W1, p.Z.W1))
    Z.W0      <- sample(c(0,1), n.W0, replace = TRUE, 
                        prob = c(1 - p.Z.W0, p.Z.W0))
    Z         <- rep(0, length.out = n)
    Z[which(W == 1)] <- Z.W1 # binary exposure
    Z[which(W == 0)] <- Z.W0
    
    # Compute probabilities to simulate the time to event outcome --------------
    p1 <- exp(alpha.T1 + beta.T1 * Z + gamma.T1 * W + eta.T1 * Z * W)     # P(T = 1 | Z, W)
    p2 <- exp(alpha.T2 + beta.T2 * Z + gamma.T2 * W + eta.T2 * Z * W)     # P(T = 2 | Z, W) = P(T = 2 | W)
    
    p1.Z1 <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * W + eta.T1 * 1 * W)  # P(T = 1 | do(Z = 1), W)
    p1.Z0 <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * W + eta.T1 * 0 * W)  # P(T = 1 | do(Z = 0), W)
    
    p2.Z1 <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * W + eta.T2 * 1 * W)  # P(T = 2 | do(Z = 1), W)
    p2.Z0 <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * W + eta.T2 * 0 * W)  # P(T = 2 | do(Z = 0), W)
    
    # Simulate the time to event outcome ---------------------------------------
    eT              <- runif(n)
    T               <- 1*(eT <= p1) + 2*((eT > p1)&(eT <= (p1 + p2)))
    T.Z1            <- 1*(eT <= p1.Z1) + 2*((eT > p1.Z1)&(eT <= (p1.Z1 + p2.Z1))) # counterfactual T^{Z = 1}
    T.Z0            <- 1*(eT <= p1.Z0) + 2*((eT > p1.Z0)&(eT <= (p1.Z0 + p2.Z0))) # counterfactual T^{Z = 0}
    
    # Simulate detection indicator D -------------------------------------------
    eD              <- runif(n)
    D               <- rep(0, times = n)
    D[T == 1]       <- 1*(eD[T == 1] <= p.detect.T1)
    D[T == 2]       <- 1*(eD[T == 2] <= p.detect.T2)
    
    D.Z1            <- rep(0, times = n) # counterfactual T^{Z = 1}
    D[T.Z1 == 1]    <- 1*(eD[T.Z1 == 1] <= p.detect.T1)
    D[T.Z1 == 2]    <- 1*(eD[T.Z1 == 2] <= p.detect.T2)
    
    D.Z0            <- rep(0, times = n) # counterfactual T^{Z = 0}
    D[T.Z0 == 1]    <- 1*(eD[T.Z0 == 1] <= p.detect.T1)
    D[T.Z0 == 2]    <- 1*(eD[T.Z0 == 2] <= p.detect.T2)
    
    # SCCS data ----------------------------------------------------------------
    n.sccs        <- sum(T != 0 & Z == 1 & D == 1)
    n.sccs.T1     <- sum(T == 1 & Z == 1 & D == 1) # number of observed vaccinated cases who had their event in period 1
    n.sccs.T2     <- sum(T == 2 & Z == 1 & D == 1)
    
    while (n.sccs == n.sccs.T1) {
      
      W         <- rbinom(n = n, size = 1, prob = pW) 
      n.W1      <- sum(W) 
      n.W0      <- sum(1 - W)
      
      Z.W1      <- sample(c(0,1), n.W1, replace = TRUE, 
                          prob = c(1 - p.Z.W1, p.Z.W1))
      Z.W0      <- sample(c(0,1), n.W0, replace = TRUE, 
                          prob = c(1 - p.Z.W0, p.Z.W0))
      Z         <- rep(0, length.out = n)
      Z[which(W == 1)] <- Z.W1
      Z[which(W == 0)] <- Z.W0
      
      p1 <- exp(alpha.T1 + beta.T1 * Z + gamma.T1 * W + eta.T1 * Z * W)     
      p2 <- exp(alpha.T2 + beta.T2 * Z + gamma.T2 * W + eta.T2 * Z * W)     
      
      p1.Z1 <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * W + eta.T1 * 1 * W)  
      p1.Z0 <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * W + eta.T1 * 0 * W)  
      
      p2.Z1 <- exp(alpha.T2 + beta.T2 * 1 + gamma.T2 * W + eta.T2 * 1 * W)  
      p2.Z0 <- exp(alpha.T2 + beta.T2 * 0 + gamma.T2 * W + eta.T2 * 0 * W)  
      
      eT              <- runif(n)
      T               <- 1*(eT <= p1) + 2*((eT > p1)&(eT <= (p1 + p2)))
      T.Z1            <- 1*(eT <= p1.Z1) + 2*((eT > p1.Z1)&(eT <= (p1.Z1 + p2.Z1))) 
      T.Z0            <- 1*(eT <= p1.Z0) + 2*((eT > p1.Z0)&(eT <= (p1.Z0 + p2.Z0))) 
      
      eD              <- runif(n)
      D               <- rep(0, times = n)
      D[T == 1]       <- 1*(eD[T == 1] <= p.detect.T1)
      D[T == 2]       <- 1*(eD[T == 2] <= p.detect.T2)
      
      D.Z1            <- rep(0, times = n) 
      D[T.Z1 == 1]    <- 1*(eD[T.Z1 == 1] <= p.detect.T1)
      D[T.Z1 == 2]    <- 1*(eD[T.Z1 == 2] <= p.detect.T2)
      
      D.Z0            <- rep(0, times = n) 
      D[T.Z0 == 1]    <- 1*(eD[T.Z0 == 1] <= p.detect.T1)
      D[T.Z0 == 2]    <- 1*(eD[T.Z0 == 2] <= p.detect.T2)
      
      n.sccs        <- sum(T != 0 & Z == 1 & D == 1)
      n.sccs.T1     <- sum(T == 1 & Z == 1 & D == 1)
      n.sccs.T2     <- sum(T == 2 & Z == 1 & D == 1)
    }
    
    pop             <- cbind(W, Z, p1, p2, p1.Z0, p1.Z1, T, T.Z1, T.Z0)       # general population of interest
    sccs            <- as.data.frame(pop[which((Z == 1)&(T != 0)&(D == 1)),]) # SCCS population

    # Quantity estimated in practice, using the SCCS data ----------------------
    marginal.contrast.est <- mean(sccs$T == 1) / mean(sccs$T == 2) # marginal contrast P(T^{Z=1} = 1 | T^{Z=1} != 0) / P(T^{Z=1} = 2 | T^{Z=1} != 0) 
    
    res <- rbind(res, cbind(marginal.contrast.est = marginal.contrast.est,
                            n = n, n.sccs = n.sccs, 
                            n.sccs.T1 = n.sccs.T1, n.sccs.T2 = n.sccs.T2,
                            p.detect.T1 = p.detect.T1, 
                            p.detect.T2 = p.detect.T2,
                            pW = pW, p.Z.W1 = p.Z.W1, p.Z.W0 = p.Z.W0,
                            beta.T1 = beta.T1, beta.T2 = beta.T2, 
                            alpha.T1 = alpha.T1, alpha.T2 = alpha.T2,
                            gamma.T1 = gamma.T1, gamma.T2 = gamma.T2,
                            eta.T1 = eta.T1, eta.T2 = eta.T2,
                            exp.beta = exp.beta,
                            marginal.contrast = marginal.contrast, 
                            part = part))
  }
  myfile  <- paste0("simulres-detectionexample-n", n, "-p.detect.T2", 
                    p.detect.T2, "-part", part, ".Rdata")
  save(res, file = myfile)
}

result <- mclapply(1:nrow(param), onerun, mc.cores = 64)

### Saving the simulation results ----------------------------------------------
RES                 <- NULL
for (p in 1:nrow(param)) {
  p.detect.T2   <- param[p, 1]
  part          <- param[p, 2]
  
  load(paste0("simulres-detectionexample-n", n, "-p.detect.T2", 
              p.detect.T2, "-part", part, ".Rdata"))
  RES <- rbind(RES, res)
}
RECAP               <- as.data.frame(RES)
ColNames            <- colnames(RECAP[,c(1:20)])
RECAP[ColNames]     <- sapply(RECAP[ColNames], as.numeric)
RECAP$p.detect.T2   <- as.factor(RECAP$p.detect.T2)
myfile              <- paste0("SimulationResults-detectionexample.Rdata")
save(RECAP, file = myfile)

### Details of the results  ----------------------------------------------------
load("SimulationResults-detectionexample.Rdata")

Nreplic <- nreplic * max(PART)
param   <- param[,-2]
param   <- param[which(duplicated(param)==FALSE)]

details.marginal.contrast   <- NULL
for (i in 1:length(param)) {
  
  RECAP1 <- RECAP[((i-1) * (Nreplic) + 1):(i * (Nreplic)), ]
  
  # mean of estimated marginal contrast
  mean.est.log.marginal.contrast <- mean(log(RECAP1$marginal.contrast.est))
  mean.est.marginal.contrast <- mean(RECAP1$marginal.contrast.est)
  
  details.marginal.contrast <- rbind(details.marginal.contrast, 
                                     c(mean.est.log.marginal.contrast = mean.est.log.marginal.contrast,
                                       mean.est.marginal.contrast = mean.est.marginal.contrast,
                                       n = as.numeric(as.character(RECAP1$n[1])), 
                                       n.sccs = mean(RECAP1$n.sccs), 
                                       pW = mean(RECAP1$pW),
                                       beta.T1 = mean(RECAP1$beta.T1), 
                                       beta.T2 = mean(RECAP1$beta.T2[1]),
                                       alpha.T1 = mean(RECAP1$alpha.T1), 
                                       alpha.T2 = mean(RECAP1$alpha.T2), 
                                       gamma.T1 = mean(RECAP1$gamma.T1), 
                                       gamma.T2 = mean(RECAP1$gamma.T2),
                                       eta.T1 = mean(RECAP1$eta.T1), 
                                       eta.T2 = mean(RECAP1$eta.T2),
                                       p.detect.T1 = mean(RECAP1$p.detect.T1),
                                       p.detect.T2 = as.numeric(as.character(RECAP1$p.detect.T2[1])), 
                                       exp.beta = mean(RECAP1$exp.beta), 
                                       beta = log(mean(RECAP1$exp.beta)), 
                                       marginal.contrast = mean(RECAP1$marginal.contrast),
                                       log.marginal.contrast = mean(log(RECAP1$marginal.contrast))))
}
details.marginal.contrast <- as.data.frame(details.marginal.contrast)

### Plot of the results --------------------------------------------------------
pdf(file = paste0("marginal.contrast_detectionexample.pdf"), width = 8, height = 6)
ggplot(data = details.marginal.contrast, aes(x = p.detect.T2, 
                                             y = mean.est.log.marginal.contrast)) + 
  geom_line(linewidth = 1)  + 
  scale_x_continuous(expression(P({D == 1}*'|'*{T == 2})), 
                     breaks = c(0.7, 0.8, 0.9, 0.98), 
                     labels = c("0.7", "0.8", "0.9", expression(P({D == 1}*'|'*{T == 1})))) +
  scale_y_continuous("Mean of logarithm of estimated marginal contrast", 
                     breaks = c(0, 0.1, 0.2, 0.3), 
                     labels =c(expression(beta), "0.1", "0.2", "0.3")) +
  theme_bw(base_size = 13) + 
  geom_hline(yintercept = 0, color = "darkgrey")
dev.off()

