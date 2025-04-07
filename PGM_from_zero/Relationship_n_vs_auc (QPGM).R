library(pheatmap)
library(XMRF)
library(pROC)
library(verification)
library(here)

# load function
# here::set_here("~/Documents/Duke University/Courses/STA 841 Categotrical Data Analysis/Poisson-Graphical-Models")
setwd("~/Documents/Duke University/Courses/STA 841 Categotrical Data Analysis/Poisson-Graphical-Models")
purrr::walk(list.files("./PGM_from_zero/helper_fns/", pattern = "*.R$", full.names=TRUE), source,.GlobalEnv)

# Node-conditional distribution via GLM: P(Z) = exp(theta B(Z) - C(Z) - D(theta))
B <- function(x) {
  return(x) # Canonical form of sufficient statistic
}

C <- function(x) {
  return(x^2) # Quadratic base measure
}

# range of n
N_range <- c(5, 10, 30, 100, 200)
AUC_range <- numeric(length = length(N_range))
Norm_range <- numeric(length = length(N_range))
p <- 5
sim <- XMRF.Sim(n = max(N_range), p = p, model = "TPGM", graph.type = "hub")
for (N in 1:length(N_range)) {
  # learn edges
  THETA = learn_edge(
    method = "QPGM",
    data = sim$X[, 1:N_range[N]],
    sufficient_statistic = B,
    base_measure = C,
    lambda = sqrt(log(p)/N_range[N]),
    M = 50
  )
  # attach result to store
  AUC_range[N] <- auc_roc_analysis(sim$B, THETA)
  Norm_range[N] <- norm(sim$B-THETA, type = "F")
  print(AUC_range)
}

# plot(N_range, AUC_range, type = "l")
