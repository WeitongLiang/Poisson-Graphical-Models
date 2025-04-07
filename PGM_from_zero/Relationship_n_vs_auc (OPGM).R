library(pheatmap)
library(XMRF)
library(pROC)
library(verification)
# library(here)

setwd("./Poisson-Graphical-Models")
purrr::walk(list.files("./PGM_from_zero/helper_fns/", pattern = "*.R$", full.names=TRUE), source,.GlobalEnv)

# range of n
N_range <- 5 * c(2,4,6,8)
AUC_range <- numeric(length = length(N_range))
Norm_range <- numeric(length = length(N_range))
p <- 4
M <- 3

for (N in 1:length(N_range)) {
  message(paste0("Current n = ", N_range[N]))
  n <- N_range[N]
  tmp_AUC <- numeric(length = M)
  tmp_norm <- numeric(length = M)
  for (m in 1:M){
    sim <- XMRF.Sim(n = n, p = p, model = "LPGM", graph.type = "scale-free")
    lmax <- lambdaMax(t(sim$X))
    lambda <- 0.01 * sqrt(log(p) / n) * lmax
    fit <- XMRF(sim$X, method = "PGM", lambda, parallel=FALSE)
    
    # attach result to store
    THETA <- fit$network[fit$opt.index][[1]]
    tmp_AUC[m] <- auc_roc_analysis(sim$B, THETA)
    tmp_norm[m] <- norm(sim$B-THETA, type = "F")
  }
  AUC_range[N] <- mean(tmp_AUC)
  Norm_range[N] <- mean(tmp_norm)
}


