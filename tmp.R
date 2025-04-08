# ---------------------------
# Simulation Experiment Code
# ---------------------------

# Load required libraries
library(glmnet)       # For penalized GLMs
library(igraph)       # Network generation and analysis
library(pROC)         # ROC curve calculation

# Simulation parameters
p <- 64               # Number of nodes (64, 100, 169, 225 as in paper)
n <- 500              # Sample size (vary to reproduce Figure 1)
omega <- -0.1         # Edge weight (negative for Poisson PGM)
n_replicates <- 50    # Number of simulation replicates

# Generate 4-nearest neighbor lattice graph
adj_matrix <- make_lattice(length = sqrt(p), dim = 2, nei = 2) |> 
  as_adjacency_matrix(sparse = FALSE)
diag(adj_matrix) <- 0 # Remove self-loops

# Function to simulate Poisson PGM data via Gibbs sampling
# Corrected function to simulate Poisson PGM data via Gibbs sampling
simulate_pgm_data <- function(adj_matrix, n, omega, burn_in = 100) {
  p <- nrow(adj_matrix)
  X <- matrix(0, nrow = n, ncol = p)
  theta <- adj_matrix * omega  # Apply edge weight
  
  # Initialize current state vector
  x_current <- rep(0, p)
  
  # Gibbs sampling
  for (iter in 1:(burn_in + n)) {
    for (s in 1:p) {
      # Compute the conditional mean for variable s
      neighbors <- which(theta[s, ] != 0)
      if (length(neighbors) == 0) {
        lambda <- exp(0)  # No neighbors
      } else {
        eta <- sum(theta[s, neighbors] * x_current[neighbors])
        lambda <- exp(eta)
      }
      # Sample new value for variable s
      x_current[s] <- rpois(1, lambda = lambda)
    }
    
    # Store after burn-in
    if (iter > burn_in) {
      X[iter - burn_in, ] <- x_current  # Correct assignment
    }
  }
  return(t(X))
}

# Function to fit node-wise GLM and evaluate edges
fit_pgm <- function(X, lambda) {
  p <- ncol(X)
  est_adj <- matrix(0, p, p)
  
  for (s in 1:p) {
    fit <- cv.glmnet(X[, -s], X[, s], family = "poisson", alpha = 1)
    coefs <- as.vector(coef(fit, s = "lambda.min"))[-1] # Remove intercept
    est_adj[s, -s] <- coefs != 0 # Binary adjacency
  }
  return(est_adj)
}

# Simulation loop
success_rates <- numeric(n_replicates)
for (rep in 1:n_replicates) {
  # Generate data
  X <- simulate_pgm_data(adj_matrix, n, omega)
  
  # Fit model (lambda chosen as per Theorem 1)
  lambda_n <- sqrt(log(p) / n)
  est_adj <- fit_pgm(X, lambda = lambda_n)
  
  # Evaluate edge recovery
  tp <- sum(est_adj == 1 & adj_matrix == 1)
  fp <- sum(est_adj == 1 & adj_matrix == 0)
  tn <- sum(est_adj == 0 & adj_matrix == 0)
  fn <- sum(est_adj == 0 & adj_matrix == 1)
  
  success_rates[rep] <- tp / (tp + fn) # True positive rate
}

# Plot results (as in Figure 1)
plot(success_rates, type = "b", xlab = "Replicate", ylab = "Success Rate",
     main = "Edge Recovery Performance")