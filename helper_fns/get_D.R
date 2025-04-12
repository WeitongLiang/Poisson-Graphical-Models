get_D <- function(sufficient_statistic, base_measure, tol = 1e-7) {
  # approximate the log-normalization function with interpolation
  
  # function to find the summation of the kernel
  log_series_sum <- function(t, max_z = 100, tol = 1e-12) {
    z_vals <- 0:max_z
    log_terms <- t * z_vals - z_vals^2
    
    # Early stopping if tail terms are negligible
    max_term <- max(log_terms)
    log_terms <- log_terms[log_terms > max_term + log(tol)]  # keep only significant terms
    A <- max(log_terms)
    
    sum_exp <- sum(exp(log_terms - A))
    return(A + log(sum_exp))
  }
  
  # points
  t <- seq(from = -100, to = 100, length.out = 1e3)
  value <- sapply(t, log_series_sum)
  
  return(splinefun(t, value))
}