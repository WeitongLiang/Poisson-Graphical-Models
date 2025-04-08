library(XMRF)
library(pROC)
library(verification)
library(glmnet)       # For penalized GLMs
library(igraph)       # Network generation and analysis

setwd("~/Documents/Duke University/Courses/STA 841 Categotrical Data Analysis/Poisson-Graphical-Models")
purrr::walk(list.files("./PGM_from_zero/helper_fns/", pattern = "*.R$", full.names=TRUE), source,.GlobalEnv)

# Function to simulate Poisson PGM data via Gibbs sampling
simulate_pgm_data <- function(adj_matrix, n, omega, burn_in = 2*1e5) {
  p <- nrow(adj_matrix)
  X <- matrix(0, nrow = n, ncol = p)  # Samples in rows, variables in columns
  theta <- adj_matrix * omega
  
  x_current <- rpois(p, lambda = 1)
  for (iter in 1:(burn_in + n)) {
    for (s in 1:p) {
      neighbors <- which(theta[s, ] != 0)
      eta <- sum(theta[s, neighbors] * x_current[neighbors])
      x_current[s] <- rpois(1, lambda = exp(eta))
    }
    if (iter > burn_in) X[iter - burn_in, ] <- x_current
  }
  return(X)  # No transposition
}

# fit Original Poisson Graphical Model
fit_Original_Poisson_Graphical_Model <- function(X, lambda) {  # Changed parameter name
  p <- ncol(X)  # Variables in columns
  adjacency <- matrix(0, p, p)
  
  for (s in 1:p) {
    fit <- glmnet(
      x = X[, -s, drop = FALSE],  # All samples, other variables
      y = X[, s],                 # Current variable
      family = "poisson",
      alpha = 1,
      lambda = lambda,
      standardize = TRUE
    )
    coefs <- as.vector(coef(fit))[-1]  # Remove intercept
    adjacency[s, -s] <- (abs(coefs) > 1e-6) * 1
  }
  return((adjacency + t(adjacency)) / 2)  # Symmetrize
}

# calculate_success_rate <- function(true_adj, pred_adj) {
#   # Only compare upper triangular parts to avoid double-counting undirected edges
#   upper_true <- true_adj[upper.tri(true_adj)]
#   upper_pred <- pred_adj[upper.tri(pred_adj)]
# 
#   # Calculate proportion of matching edges
#   sum(upper_true == upper_pred) / length(upper_true)
# }

calculate_success_rate <- function(true_adj, pred_adj) {
  all(true_adj[upper.tri(true_adj)] == pred_adj[upper.tri(pred_adj)]) * 1
}

results_list <- list() # storage results
# N_range <- c(20, 50, 70, 80, 100 * seq(from = 1, to = 15, by = 2))
# N_range <- 1e5 * seq(from = 10, to = 15, by = 2)
N_range <- c(50, 100*(1:10), 2000 * (1: 5), 1e4 * (1:5))
Norm_range <- numeric(length = length(N_range))
Success_range <- numeric(length = length(N_range))
# p_values <- c(64, 100, 169, 225)
p_values <- c(4, 9, 16)
M <- 50
omega <- -0.1         # Edge weight (negative for Poisson PGM)

for (p in p_values) {
  message(paste("\nRunning experiments for p =", p))
  # Generate 4-nearest neighbor lattice graph
  adj_matrix <- make_lattice(length = sqrt(p), dim = 2, nei = 1) |> 
    as_adjacency_matrix(sparse = FALSE)
  diag(adj_matrix) <- 0 # Remove self-loops
  
  for (N in 1:length(N_range)) {
    # message(paste0("Current n = ", N_range[N]))
    n <- N_range[N]
    tmp_AUC <- numeric(length = M)
    tmp_norm <- numeric(length = M)
    tmp_success <- numeric(length = M)
    for (m in 1:M){
      cat(paste("Calculating success rate in progress (n = ", n, "): ", floor(100*(m/M)), "%", collapse=""),"\r")
      flush.console()
      
      X <- simulate_pgm_data(adj_matrix, n, omega)
      THETA <- fit_Original_Poisson_Graphical_Model(X = X, lambda = sqrt(log(p) / n))
      
      # attach result to store
      tmp_norm[m] <- norm(adj_matrix-THETA, type = "F")
      tmp_success[m] <- calculate_success_rate(adj_matrix, THETA)
    }
    Success_range[N] <- mean(tmp_success)
    Norm_range[N] <- mean(tmp_norm)
  }
  
  # Store results for this p
  results_list[[as.character(p)]] <- list(
    p = p,
    N_range = N_range,
    Success_range = Success_range
  )
}

plot(N_range[1:20], results_list[["4"]]$Success_range[1:20], type = "l", col = 1, 
     ylim = c(0, 1), xlab = "Sample Size (n)", ylab = "Success Rate",
     main = "Edge Recovery Performance")
lines(N_range, results_list[["9"]]$Success_range, type = "l", col = 2)
lines(N_range, results_list[["16"]]$Success_range, type = "l", col = 3)
lines(N_range, results_list[["225"]]$Success_range, type = "l", col = 4)
legend("bottomright", legend = paste("p =", p_values), col = 1:4, lty = 1, pch = 1)

# saveRDS(results_list, file = "./Edge_recovery_performance.RData")
# tmp <- readRDS(file = "Edge_recovery_performance.RData")
