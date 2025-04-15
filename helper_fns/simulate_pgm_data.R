simulate_pgm_data <- function(adj_matrix, n, omega, burn_in = 2e5, rate_cap = 50) {
  p <- nrow(adj_matrix)
  X <- matrix(0, nrow = n, ncol = p)
  theta <- adj_matrix * omega
  x_current <- rpois(p, lambda = 1)
  for (iter in 1:(burn_in + n)) {
    for (s in 1:p) {
      neighbors <- which(theta[s, ] != 0)
      eta <- sum(theta[s, neighbors] * x_current[neighbors])
      lambda <- min(exp(eta), rate_cap)
      x_current[s] <- rpois(1, lambda = lambda)
    }
    if (iter > burn_in) X[iter - burn_in, ] <- x_current
  }
  return(X)
}
