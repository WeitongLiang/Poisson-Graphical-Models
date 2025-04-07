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
B <- function(X,R,R0=0) {
  Bx = X
  Bx[X>R] = (R+R0)/2
  ind = X>R0 & X<=R
  Bx[ind] =(-X[ind]^2 +2*R*X[ind]-R0^2)/(2*(R-R0))
  return(Bx)
}

C <- function(x) {
  return(log(factorial(x)))
}

# range of n
N_range <- c(10, 20, 50, 100)
AUC_range <- numeric(length = length(N_range))
Norm_range <- numeric(length = length(N_range))
p <- 5
sim <- XMRF.Sim(n = max(N_range), p = p, model = "TPGM", graph.type = "hub")
for (N in 1:length(N_range)) {
  # learn edges
  THETA = learn_edge(
    method = "Sub-linear PGM",
    data = sim$X[, 1:N_range[N]],
    sufficient_statistic = function(x) {B(x, R = max(sim$X))},
    base_measure = C,
    lambda = sqrt(log(p)/N_range[N]),
    M = 50
  )
  # attach result to store
  AUC_range[N] <- auc_roc_analysis(sim$B, THETA)
  Norm_range[N] <- norm(sim$B-THETA, type = "F")
  print(Norm_range[N])
}

# plot(N_range, AUC_range, type = "l")
