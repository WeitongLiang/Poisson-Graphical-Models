fit_Sublinear_Poisson_Graphical_Model_optim <- function(X, lambda, threshold = 1e-1){
  #' Fit Sublinear Poisson Graphical Model (SPGM)
  #' 
  #' @param X n-by-p count matrix (rows = samples, cols = variables)
  #' @param lambda L1-penalty
  #' @return p-by-p symmetric parameter matrix (or optionally binary adjacency matrix)
  
  cat("Start training SPGM...\n")
  
  # ===== Sublinear sufficient statistic =====
  B <- function(x, R0, R) {
    out <- numeric(length(x))
    
    not_na <- !is.na(x)
    
    # Case 1: x <= R0
    out[not_na & x <= R0] <- x[not_na & x <= R0]
    
    # Case 2: R0 < x < R
    idx_middle <- not_na & x > R0 & x < R
    out[idx_middle] <- -1 / (2 * (R - R0)) * x[idx_middle]^2 +
      R / (R - R0) * x[idx_middle] -
      R0^2 / (2 * (R - R0))
    
    # Case 3: x >= R
    out[not_na & x >= R] <- (R + R0) / 2
    
    return(out)
  }
  
  
  # ===== Log-normalizing constant (numerically approximated) =====
  LogNormalizing_spgm <- function(th, R0, R) {
    tryCatch({
      n_vals <- 0:R
      
      B_vec <- numeric(length(n_vals))
      B_vec[n_vals <= R0] <- n_vals[n_vals <= R0]
      middle_idx <- n_vals > R0 & n_vals < R
      B_vec[middle_idx] <- -1/(2*(R - R0)) * n_vals[middle_idx]^2 +
        R / (R - R0) * n_vals[middle_idx] -
        R0^2 / (2 * (R - R0))
      B_vec[n_vals >= R] <- (R + R0) / 2
      
      log_terms <- th * B_vec - lgamma(n_vals + 1)
      m <- max(log_terms)
      log_sum <- m + log(sum(exp(log_terms - m)))
      
      return(log_sum)
    }, error = function(e) {
      return(log(1e10))
    })
  }
  
  # ===== Wrapped versions for batch use =====
  R0 <- 0
  R <- max(X)
  
  B_sublinear <- function(x) B(x, R0, R)
  
  D_sublinear <- function(theta_vec) {
    sapply(theta_vec, function(th) LogNormalizing_spgm(th, R0, R))
  }
  
  # ===== Objective function for nodewise regression =====
  neg_loglik_l1_qpgm <- function(beta, x, y, lambda) {
    eta <- as.vector(x %*% beta)
    d_val <- D_sublinear(eta)
    b_val <- B_sublinear(y)
    
    # Fallback for numerical problems
    if (any(is.nan(d_val)) || any(is.infinite(d_val)) || any(is.na(d_val)) ||
        any(is.nan(b_val)) || any(is.infinite(b_val)) || any(is.na(b_val))) {
      return(1e10)
    }
    
    loglik <- sum(b_val * eta - log(factorial(y)) - d_val)
    penalty <- lambda * sum(abs(beta))
    out <- -loglik + penalty
    
    if (!is.finite(out)) return(1e10)
    return(out)
  }
  
  # ===== Main optimization loop =====
  p <- ncol(X)
  n <- nrow(X)
  THETA <- matrix(0, nrow = p, ncol = p)
  
  for (s in 1:p) {
    cat(sprintf("Training node %d out of %d nodes... Progress = %.2f%%\n", s, p, s/p * 100))
    flush.console()
    
    X_others <- X[, -s]
    Y <- X[, s]
    
    converged <- FALSE
    maxit <- 200
    while (!converged) {
      maxit <- maxit * 2
      fit <- optim(
        par = rep(0, p - 1),
        fn = neg_loglik_l1_qpgm,
        x = X_others,
        y = Y,
        lambda = lambda,
        method = "BFGS",
        control = list(maxit = maxit)
      )
      converged <- (fit$convergence == 0)
      if (!converged) {
        message(sprintf("Not converging, increasing maxit to %d", maxit))
      }
    }
    
    THETA[s, ] <- THETA[s, ] + append(fit$par, 0, after = s - 1)
    THETA[, s] <- THETA[, s] + append(fit$par, 0, after = s - 1)  # symmetric
  }
  
  # Optional binarization of edges
  # THETA_bin <- (abs(THETA) > threshold) * 1
  
  return(THETA)
}
