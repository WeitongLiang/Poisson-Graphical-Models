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