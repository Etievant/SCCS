## -----------------------------------------------------------------------------
## Description: This script illustrates the Acceleration Effect presented in
##              Section 5.2 in Etievant, Gail and Follmann (2024)
##
##              The simulation is in Section 8.2 in the Main Document
##
##              We consider a scenario where W is the empty set
## -----------------------------------------------------------------------------

### load packages --------------------------------------------------------------
library(tidyr)
library(parallel)
library(ggplot2)

### Set fixed parameters -------------------------------------------------------
n         <- 10^6             # (overall) population size
alpha.T1  <- log(0.00008)     # log baseline risk in period 1 
alpha.T2  <- log(0.00008)     # log baseline risk in period 2 
beta.T1   <- 0              
beta.T2   <- 0
p.Z       <- 0.3              # P(Z = 1)

### Compute counterfactual probabilities analytically --------------------------
P.T1.doZ0       <- exp(alpha.T1 + beta.T1 * 0)
P.T2.doZ0       <- exp(alpha.T2 + beta.T2 * 0)

### Varying parameter ----------------------------------------------------------
P21       <- seq(from = 0, to = 0.2, length.out = 3) # P(T^{Z=1} = 1 | T^{Z=0} = 2)

### Set replication parameters -------------------------------------------------
PART      <- 1:20   # splitting the replications
nreplic   <- 250     # total number of replications: Nreplic = nreplic * max(PART)
set.seed(12345)
seed      <- round(abs(rnorm(max(PART)) * 10^5))
param     <- as.data.frame(expand_grid(p21 = P21, part = PART))

### Function to replicate the example over the replications --------------------
onerun <- function(p) {
  p21               <- param[p, 1]
  part              <- param[p, 2]
  
  set.seed(seed[part])
  res <- NULL
  
  # Compute other counterfactual probabilities analytically --------------------
  P.T1.doZ1     <- P.T1.doZ0 + p21 * P.T2.doZ0 # P(T^{Z=1} = 1)
  P.T2.doZ1     <- (1 - p21) * P.T2.doZ0       # P(T^{Z=1} = 2)
  
  exp.beta.T1   <- P.T1.doZ1 / P.T1.doZ0 # marginal causal risk ratio over period 1
  exp.beta.T2   <- P.T2.doZ1 / P.T2.doZ0 # marginal causal risk ratio over period 2
  exp.beta.T12  <- (P.T1.doZ1 + P.T2.doZ1) / (P.T1.doZ0 + P.T2.doZ0) # marginal causal risk ratio over period 1 + period 2
  
  # Focus is on the marginal contrast as W = empty set
  marginal.contrast   <- P.T1.doZ1 / P.T2.doZ1
  
  for (nrep in 1:nreplic){
    
    # Simulate exposure Z ------------------------------------------------------
    Z <- rbinom(n, 1, p.Z) # binary exposure
    
    # Simulate the time to event outcome ---------------------------------------
    
    # First simulate the counterfactual outcome under do(Z = 0)
    p1.Z0             <- exp(alpha.T1 + beta.T1 * 0)
    p2.Z0             <- exp(alpha.T2 + beta.T2 * 0)
    eT                <- runif(n)
    T.Z0              <- 1*(eT <= p1.Z0) + 2*((eT > p1.Z0)&(eT <= (p1.Z0 + p2.Z0))) # counterfactual T^{Z = 0}
    
    # Then the counterfactual outcome under do(Z = 1)
    T.Z1              <- T.Z0
    e21               <- runif(sum(T.Z0 == 2))
    T.Z1[(T.Z0 == 2)] <- 1*(e21 <= p21) + 2*((e21 > p21)) # ATS
    
    T                 <- T.Z0
    T[Z == 1]         <- T.Z1[Z == 1] # by consistency
    
    # SCCS data ----------------------------------------------------------------
    n.sccs          <- sum(T != 0 & Z == 1)
    pop             <- cbind(Z, p1.Z0, p2.Z0, T, T.Z1, T.Z0)         # general population of interest
    sccs            <- as.data.frame(pop[which((Z == 1)&(T != 0)),]) # SCCS population
    n.sccs.T1       <- sum(sccs$T == 1) # number vaccinated cases who had their event in period 1
    n.sccs.T2       <- sum(sccs$T == 2)
    
    # Quantity estimated in practice, using the SCCS data ----------------------
    marginal.contrast.est <- mean(sccs$T == 1) / mean(sccs$T == 2) # marginal contrast P(T^{Z=1} = 1 | T^{Z=1} != 0) / P(T^{Z=1} = 2 | T^{Z=1} != 0) 
  
    res <- rbind(res, cbind(marginal.contrast.est = marginal.contrast.est,
                            n = n, n.sccs = n.sccs, 
                            n.sccs.T1 = n.sccs.T1, n.sccs.T2 = n.sccs.T2, 
                            p.Z = p.Z, 
                            p21 = p21,
                            beta.T1 = beta.T1,
                            beta.T2 = beta.T2,
                            alpha.T1 = alpha.T1, alpha.T2 = alpha.T2,
                            sum.T1.Z1 = sum(T.Z1 == 1), 
                            sum.T1.Z0 = sum(T.Z0 == 1),
                            sum.T2.Z1 = sum(T.Z1 == 2), 
                            sum.T2.Z0 = sum(T.Z0 == 2),
                            cases.Z0.Z1 = mean(which(T.Z0 != 0) == which(T.Z1 != 0)),
                            exp.beta.T1 = exp.beta.T1, 
                            exp.beta.T2 = exp.beta.T2,
                            exp.beta.T12 = exp.beta.T12,
                            marginal.contrast = marginal.contrast,
                            part = part))
  }
  myfile  <- paste0("simulres-ats-n", n, "-p21", p21, "-part", part, 
                    ".Rdata")
  save(res, file = myfile)
}

result <- mclapply(1:nrow(param), onerun, mc.cores = 64)

### Saving the simulation results ----------------------------------------------
RES             <- NULL
for (p in 1:nrow(param)) {
  p21       <- param[p, 1]
  part      <- param[p, 2]
  
  load(paste0("simulres-ats-n", n, "-p21", p21, "-part", part, 
              ".Rdata"))
  RES <- rbind(RES, res)
}
RECAP           <- as.data.frame(RES)
ColNames        <- colnames(RECAP[,c(1:21)])
RECAP[ColNames] <- sapply(RECAP[ColNames], as.numeric)
RECAP$n         <- as.factor(RECAP$n)
myfile          <- paste0("SimulationResults-ats.Rdata")
save(RECAP, file = myfile)

### Details of the results  ----------------------------------------------------
load("SimulationResults-ats.Rdata")

Nreplic   <- nreplic * max(PART)
param     <- param[,-2]
param     <- param[which(duplicated(param)==FALSE)]

details.marginal.contrast   <- NULL
for (i in 1:length(param)) {
  
  RECAP1 <- RECAP[((i-1) * (Nreplic) + 1):(i * (Nreplic)), ]
  
  exp.beta.T1 = mean(RECAP1$exp.beta.T1)
  exp.beta.T2 = mean(RECAP1$exp.beta.T2)
  exp.beta.T12 = mean(RECAP1$exp.beta.T12)

  # mean of estimated marginal contrast
  mean.est.log.marginal.contrast  <- mean(log(RECAP1$marginal.contrast))
  mean.est.marginal.contrast      <- mean(RECAP1$marginal.contrast)
  
  details.marginal.contrast <- rbind(details.marginal.contrast, 
                                     c(mean.est.log.marginal.contrast = mean.est.log.marginal.contrast,
                                       mean.est.marginal.contrast = mean.est.marginal.contrast,
                                       n = as.numeric(as.character(RECAP1$n[1])), 
                                       n.sccs = mean(RECAP1$n.sccs), 
                                       n.sccs.T1 = mean(RECAP1$n.sccs.T1),
                                       n.sccs.T2 = mean(RECAP1$n.sccs.T2), 
                                       p.Z = mean(RECAP1$p.Z), 
                                       beta.T1 = mean(RECAP1$beta.T1), 
                                       beta.T2 = mean(RECAP1$beta.T2[1]),
                                       alpha.T1 = mean(RECAP1$alpha.T1), 
                                       alpha.T2 = mean(RECAP1$alpha.T2), 
                                       p21 = as.numeric(as.character(RECAP1$p21[1])),
                                       exp.beta.period1 = mean(exp.beta.T1),
                                       exp.beta.period2 = mean(exp.beta.T2),
                                       exp.beta.period12 = mean(exp.beta.T12),
                                       beta.period1 = mean(log(exp.beta.T1)),
                                       beta.period2 = mean(log(exp.beta.T2)),
                                       beta.period12 = mean(log(exp.beta.T12)),
                                       marginal.contrast = mean(RECAP1$marginal.contrast),
                                       log.marginal.contrast = mean(log(RECAP1$marginal.contrast))))
}
details.marginal.contrast <- as.data.frame(details.marginal.contrast)

### Plot of the results --------------------------------------------------------
results <- rbind(cbind(value = details.marginal.contrast$mean.est.log.marginal.contrast, 
                       p21 = details.marginal.contrast$p21, 
                       Contrast = "marginal.contrast"),
                 cbind(value = details.marginal.contrast$beta.period1, 
                       p21 = details.marginal.contrast$p21, 
                       Contrast = "effect.period1"),
                 cbind(value = details.marginal.contrast$beta.period2, 
                       p21 = details.marginal.contrast$p21, 
                       Contrast = "effect.period2"),
                 cbind(value = details.marginal.contrast$beta.period12, 
                       p21 = details.marginal.contrast$p21, 
                       Contrast = "effect.period12"))

results           <- as.data.frame(results)
results$value     <- as.numeric(results$value)
results$p21       <- as.numeric(results$p21)
results$Contrast  <- as.factor(results$Contrast)

legend_labels     <- c(marginal.contrast = "Mean of logarithm of estimated marginal contrast",
                       effect.period1 = expression(log(P(T^{Z==1} == 1)/P(T^{Z==0} == 1))), 
                       effect.period2 = expression(log(P(T^{Z==1} == 2) / P(T^{Z==0} == 2))),
                       effect.period12 = expression(log((P(T^{Z==1} == 1) + P(T^{Z==1} == 2)) / (P(T^{Z==0} == 1) + P(T^{Z==0} == 2)))))
                  
pdf(file = paste0("marginal.contrast_ats.pdf"), width = 8, height = 6)
ggplot(data = results, aes(x = p21, y = value, linetype = Contrast, 
                           colour = Contrast)) + 
  geom_line(linewidth = 0.8) + 
  xlab(expression(P({T^{Z==1} == 1}*'|'*{T^{Z==0} == 2}))) + 
  ylab("") + 
  theme_bw(base_size = 13) + 
  scale_colour_manual(name = "", 
                      values = c("orange", "pink", "lightblue", "black"), 
                      labels = legend_labels) +
  scale_linetype_manual(name = "", 
                        values = c("dotdash", "dashed", "dotted", "solid"), 
                        labels = legend_labels)
dev.off()

