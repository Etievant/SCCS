## -----------------------------------------------------------------------------
## Description: This script illustrates Scenario 1 with violation of Condition 
##              (ii) presented in Section 5 in Etievant, Gail and Follmann 
##              (2024)
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
alpha.T1  <- log(0.00004)     # log baseline risk in period 1 
#alpha.T2                     # log baseline risk in period 2 
gamma.T1  <- 0.2              # we use g_1: w -> gamma.T1 * w
#gamma.T2                     # we use g_2: w -> gamma.T2 * w
beta.T1   <- 0.3              # no interaction between Z and W
exp.beta  <- exp(beta.T1)
p.Z.W1    <- 0.2              # P(Z = 1 | W = 1)
p.Z.W0    <- 0.35             # P(Z = 1 | W = 0)

### Compute counterfactual probabilities analytically --------------------------
P.T1.doZ0.W1    <- exp(alpha.T1 + gamma.T1 * 1) # P(T^{Z=0} = 1 | W = 1)
P.T1.doZ0.W0    <- exp(alpha.T1 + gamma.T1 * 0) # P(T^{Z=0} = 1 | W = 0)
P.T1.doZ1.W1    <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 1) # P(T^{Z=1} = 1 | W = 1)
P.T1.doZ1.W0    <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 0) # P(T^{Z=1} = 1 | W = 0)

P.T1.doZ1       <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 1) * pW + exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * 0) * (1 - pW) # P(T^{Z=1} = 1)
P.T1.doZ0       <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * 1) * pW + exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * 0) * (1 - pW) # P(T^{Z=0} = 1)

### Varying parameters ---------------------------------------------------------
GAMMA.T2 <- c(gamma.T1, 0.3, 0.4)
ALPHA.T2 <- seq(from = alpha.T1, to = log(0.00010), length.out = 12)

### Set replication parameters -------------------------------------------------
PART      <- 1:20   # splitting the replications
nreplic   <- 250    # total number of replications: Nreplic = nreplic * max(PART)
set.seed(1234)
seed      <- round(abs(rnorm(max(PART)) * 10^5))
param     <- as.data.frame(expand_grid(gamma.T2 = GAMMA.T2, alpha.T2 = ALPHA.T2, 
                                       part = PART))

### Function to replicate the example over the replications --------------------
onerun <- function(p) {
  alpha.T2  <- param[p, 2]
  gamma.T2  <- param[p, 1]
  part      <- param[p, 3]
  
  set.seed(seed[part])
  res <- NULL
  
  # Compute other counterfactual probabilities analytically --------------------
  P.T2.doZ1.W1    <- exp(alpha.T2 + gamma.T2 * 1) # P(T^{Z=1} = 2 | W = 1)
  P.T2.doZ0.W1    <- exp(alpha.T2 + gamma.T2 * 1) # P(T^{Z=0} = 2 | W = 1) 
  P.T2.doZ1.W0    <- exp(alpha.T2 + gamma.T2 * 0) # P(T^{Z=1} = 2 | W = 0)
  P.T2.doZ0.W0    <- exp(alpha.T2 + gamma.T2 * 0) # P(T^{Z=0} = 2 | W = 0)
  
  P.T2.doZ1       <- exp(alpha.T2 + gamma.T2 * 1) * pW + exp(alpha.T2 + gamma.T2 * 0) * (1 - pW) # P(T^{Z=1} = 2)
  P.T2.doZ0       <- exp(alpha.T2 + gamma.T2 * 1) * pW + exp(alpha.T2 + gamma.T2 * 0) * (1 - pW) # P(T^{Z=0} = 2)

  # Focus is on the marginal contrast as the effect is homogeneous across strata
  marginal.contrast <- (P.T1.doZ1.W1 * p.Z.W1 * pW / 
                          (p.Z.W1 * pW + p.Z.W0 * (1 - pW)) + 
                          P.T1.doZ1.W0 * p.Z.W0 * (1 - pW) / 
                          (p.Z.W1 * pW + p.Z.W0 * (1 - pW))) / 
    (P.T2.doZ1.W1 * p.Z.W1 * pW / (p.Z.W1 * pW + p.Z.W0 * (1 - pW)) + 
       P.T2.doZ1.W0 * p.Z.W0 * (1 - pW) / (p.Z.W1 * pW + p.Z.W0 * (1 - pW))) # marginal quantity estimated in practice from the SCCS data (analytic value)
  
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
    p1 <- exp(alpha.T1 + beta.T1 * Z + gamma.T1 * W)    # P(T = 1 | Z, W)
    p2 <- exp(alpha.T2 + gamma.T2 * W)               # P(T = 2 | Z, W) = P(T = 2 | W)
    
    p1.Z1 <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * W) # P(T = 1 | do(Z = 1), W)
    p1.Z0 <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * W) # P(T = 1 | do(Z = 0), W)
    
    # Simulate the time to event outcome ---------------------------------------
    eT              <- runif(n)
    T               <- 1*(eT <= p1) + 2*((eT > p1)&(eT <= (p1 + p2)))
    T.Z1            <- 1*(eT <= p1.Z1) + 2*((eT > p1.Z1)&(eT <= (p1.Z1 + p2))) # counterfactual T^{Z = 1}
    T.Z0            <- 1*(eT <= p1.Z0) + 2*((eT > p1.Z0)&(eT <= (p1.Z0 + p2))) # counterfactual T^{Z = 0}
    
    # SCCS data ----------------------------------------------------------------
    n.sccs          <- sum(T != 0 & Z == 1)
    n.sccs.W1       <- sum(T != 0 & Z == 1 & W == 1)
    n.sccs.W0       <- sum(T != 0 & Z == 1 & W == 0)
    n.sccs.T1.W1    <- sum(T == 1 & Z == 1 & W == 1)
    n.sccs.T1.W0    <- sum(T == 1 & Z == 1 & W == 0)
    
    while ((n.sccs.W0 == n.sccs.T1.W0)|(n.sccs.W1 == n.sccs.T1.W1)|
           (n.sccs.T1.W0 == 0)|(n.sccs.T1.W1 == 0)) {
      
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
      
      p1 <- exp(alpha.T1 + beta.T1 * Z + gamma.T1 * W)    
      p2 <- exp(alpha.T2 + gamma.T2 * W)              
      
      p1.Z1 <- exp(alpha.T1 + beta.T1 * 1 + gamma.T1 * W) 
      p1.Z0 <- exp(alpha.T1 + beta.T1 * 0 + gamma.T1 * W) 
      
      eT              <- runif(n)
      T               <- 1*(eT <= p1) + 2*((eT > p1)&(eT <= (p1 + p2)))
      T.Z1            <- 1*(eT <= p1.Z1) + 2*((eT > p1.Z1)&(eT <= (p1.Z1 + p2))) 
      T.Z0            <- 1*(eT <= p1.Z0) + 2*((eT > p1.Z0)&(eT <= (p1.Z0 + p2))) 
      
      n.sccs          <- sum(T != 0 & Z == 1)
      n.sccs.W1       <- sum(T != 0 & Z == 1 & W == 1)
      n.sccs.W0       <- sum(T != 0 & Z == 1 & W == 0)
      n.sccs.T1.W1    <- sum(T == 1 & Z == 1 & W == 1)
      n.sccs.T1.W0    <- sum(T == 1 & Z == 1 & W == 0)
    }
    
    pop             <- cbind(W, Z, p1, p2, p1.Z0, p1.Z1, T, T.Z1, T.Z0) # general population of interest
    sccs            <- as.data.frame(pop[which((Z == 1)&(T != 0)),])    # sccs population
    
    # Quantity estimated in practice, using the SCCS data ----------------------
    marginal.contrast.est <- mean(sccs$T == 1) / mean(sccs$T == 2) # marginal contrast P(T^{Z=1} = 1 | T^{Z=1} != 0) / P(T^{Z=1} = 2 | T^{Z=1} != 0) 
    
    res <- rbind(res, cbind(marginal.contrast.est = marginal.contrast.est,
                            n = n, n.sccs = n.sccs, 
                            n.sccs.W1 = n.sccs.W1,
                            n.sccs.W0 = n.sccs.W0, 
                            n.sccs.T1.W1 = n.sccs.T1.W1, 
                            n.sccs.T1.W0 = n.sccs.T1.W0, 
                            pW = pW, p.Z.W1 = p.Z.W1, p.Z.W0 = p.Z.W0,
                            beta.T1 = beta.T1, beta.T2 = 0,
                            alpha.T1 = alpha.T1, alpha.T2 = alpha.T2,
                            gamma.T1 = gamma.T1, gamma.T2 = gamma.T2,
                            exp.beta = exp.beta, 
                            marginal.contrast = marginal.contrast,
                            part = part))
  }
  myfile  <- paste0("simulres-example2-n", n, "-alpha.T2", alpha.T2, 
                    "-gamma.T2", gamma.T2, "-part", part, ".Rdata")
  save(res, file = myfile)
}

result <- mclapply(1:nrow(param), onerun, mc.cores = 64)

### Saving the simulation results ----------------------------------------------
RES             <- NULL
for (p in 1:nrow(param)) {
  alpha.T2  <- param[p, 2]
  gamma.T2  <- param[p, 1]
  part      <- param[p, 3]
  
  load(paste0("simulres-example2-n", n, "-alpha.T2", alpha.T2, 
              "-gamma.T2", gamma.T2, "-part", part, ".Rdata"))
  RES <- rbind(RES, res)
}
RECAP           <- as.data.frame(RES)
ColNames        <- colnames(RECAP[,c(1:19)])
RECAP[ColNames] <- sapply(RECAP[ColNames], as.numeric)
RECAP$n         <- as.factor(RECAP$n)
RECAP$alpha.T2  <- as.factor(RECAP$alpha.T2)
RECAP$gamma.T2  <- as.factor(RECAP$gamma.T2)
myfile          <- paste0("SimulationResults-example2.Rdata")
save(RECAP, file = myfile)

### Details of the results  ----------------------------------------------------
load("SimulationResults-example2.Rdata")

Nreplic <- nreplic * max(PART)
param   <- param[,-3]
param   <- param[which(duplicated.matrix(param)==FALSE),]

details.contrast.marginal   <- NULL
for (i in 1:nrow(param)) {
  
  RECAP1 <- RECAP[((i-1) * (Nreplic) + 1):(i * (Nreplic)), ]
  
  # mean of estimated marginal contrast
  mean.est.log.marginal.contrast  <- mean(log(RECAP1$marginal.contrast.est))
  mean.est.marginal.contrast      <- mean(RECAP1$marginal.contrast.est)
  
  details.contrast.marginal <- rbind(details.contrast.marginal, 
                                     c(mean.est.log.marginal.contrast = mean.est.log.marginal.contrast,
                                       mean.est.marginal.contrast = mean.est.marginal.contrast,
                                       n = as.numeric(as.character(RECAP1$n[1])), 
                                       n.sccs = mean(RECAP1$n.sccs), 
                                       n.sccs.W1 = mean(RECAP1$n.sccs.W1),
                                       n.sccs.W0 = mean(RECAP1$n.sccs.W0), 
                                       n.sccs.T1.W1 = mean(RECAP1$n.sccs.T1.W1),
                                       n.sccs.T1.W0 = mean(RECAP1$n.sccs.T1.W0), 
                                       pW = mean(RECAP1$pW), 
                                       p.Z.W1 = mean(RECAP1$p.Z.W1),
                                       p.Z.W0 = mean(RECAP1$p.Z.W0),
                                       beta.T1 = mean(RECAP1$beta.T1),
                                      alpha.T1 = mean(RECAP1$alpha.T1), 
                                      alpha.T2 = as.numeric(as.character(RECAP1$alpha.T2[1])),
                                      gamma.T1 = mean(RECAP1$gamma.T1), 
                                      gamma.T2 = as.numeric(as.character(RECAP1$gamma.T2[1])),
                                      exp.beta = mean(RECAP1$exp.beta), 
                                      marginal.contrast = mean(RECAP1$marginal.contrast),
                                      log.marginal.contrast = mean(log(RECAP1$marginal.contrast))))
}

details.contrast.marginal <- as.data.frame(details.contrast.marginal)

# Plot of the results ----------------------------------------------------------
details.contrast.marginal$gamma.T2 <- as.factor(details.contrast.marginal$gamma.T2)
details.contrast.marginal$gamma.T2 <- factor(details.contrast.marginal$gamma.T2, 
                                             labels = c("gamma[2]==gamma[1]", 
                                                        "gamma[2]==0.3", 
                                                        "gamma[2]==0.4"))
pdf(file = paste0("marginal.contrast_example2.pdf"), width = 8, height = 6)
ggplot(data = details.contrast.marginal, aes(x = alpha.T2, 
                                             y = mean.est.log.marginal.contrast)) + 
  geom_line(linewidth = 1)  + 
  scale_x_continuous(expression(alpha[2]), 
                     breaks = c(-10.12663, -10, -9.6, -9.2), 
                     labels = c(expression(alpha[1]), "-10", "-9.6", "-9.2")) +
  scale_y_continuous("Mean of logarithm of estimated marginal contrast", 
                     breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.3), 
                     labels =c("-0.6", "-0.4", "-0.2", "0", "0.2", 
                               expression(beta))) +
  facet_wrap(~ gamma.T2, labeller = label_parsed) + theme_bw(base_size = 13)
dev.off()
