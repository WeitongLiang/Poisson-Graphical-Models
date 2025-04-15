# ===== Function: fit QPGM with optim() =====
fit_Quadratic_Poisson_Graphical_Model_optim <- function(X, lambda, threshold = 1e-1){
  #' Fit Quadratic Poisson Graphical Model
  #' 
  #' @param X a n-by-p matrix, where n is the number of samples, p is the number of nodes (variables)
  #' @param lambda L1-regularization coefficient
  #' @return p-by-p adjacency matrix indicating the edeges (1 means edge, 0 means no edge)
  cat("Start training QPGM...")
  
  # optimization target 
  neg_loglik_l1_qpgm <- function(beta, x, y, lambda) {
    eta <- x %*% beta
    loglik <- sum(y * eta - y^2 - D(eta))
    penalty <- lambda * sum(abs(beta))
    return(-loglik + penalty)
  }
  # prepare params
  p <- ncol(X); n <- nrow(X)
  THETA <- matrix(0, nrow = p, ncol = p)
  for (s in 1:p) {
    cat(paste0("Training node ", s, " out of ", p," nodes... Progress = ", round(s/p*100,2), "%\n"))
    flush.console()
    # prepare data
    X_others <- X[, -s]
    Y <- X[, s]
    # fit node-conditional distribution
    converged <- FALSE; maxit <- 200
    while(!converged){
      maxit <- maxit * 2
      fit <- optim(
        par = rep(0, p-1),
        fn = neg_loglik_l1_qpgm,
        x = X_others,
        y = Y,
        lambda = lambda,
        method = "BFGS",  # or "L-BFGS-B" for box constraints
        control = list(maxit = maxit)
      )
      converged <- ifelse(fit$convergence==0, TRUE, FALSE)
      if (!converged){
        message(paste0("Not converging, training at maxit=",maxit*2))
      }
    }
    # append value to THETA
    THETA[s, ] <- THETA[s, ] + append(fit$par, 0, after = s-1)
    THETA[, s] <- THETA[, s] + append(fit$par, 0, after = s-1) # symmetry
  }
  # print(THETA)
  # THETA <- (THETA>threshold) + 0
  return(THETA)
}
