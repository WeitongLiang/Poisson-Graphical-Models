# ===== Function: fit SPGM with optim() =====
fit_Sublinear_Poisson_Graphical_Model_optim <- function(X, lambda, threshold = 1e-1){
  #' Fit Sublinear Poisson Graphical Model
  #' 
  #' @param X a n-by-p matrix, where n is the number of samples, p is the number of nodes (variables)
  #' @param lambda L1-regularization coefficient
  #' @return p-by-p adjacency matrix indicating the edeges (1 means edge, 0 means no edge)
  cat("Start training QPGM...")
  
  Sublinear <- function(x, R0, R){
    if (x==0) {return(0)}
    else if (x>=R) {
      return((R+R0)/2)
    } else {
      return(
        -1/(2*(R-R0))*x^2 + R/(R-R0)*x - R0^2/(2*(R-R0))
      )
    }
  }
  
  LogNormalizing_spgm <- function(th, R0, R){
    part_1 <- sum(sapply(0:R0, function(n){exp(th*n)/factorial(n)}))
    part_2 <- sum(sapply(R0+1:R, function(n){
      exp((-1/(2*(R-R0))*n^2 + R/(R-R0)*n - R0^2/(2*(R-R0)))*th)/factorial(n)
    }))
    part_3 <- exp(th * ((R+R0)/2)) * (1-sum(sapply(0:R, function(n){1/factorial(n)})))
    return(log(part_1 + part_2 + part_3))
  }
  
  B_sublinear <- function(x) {
    return(B(x, R0=0, R = max(X)))  # R0 is set as 0
  }
  D_sublinear <- function(theta) {
    return(LogNormalizing_spgm(theta, R0=0, R = max(X)))
  }
  
  # optimization target 
  neg_loglik_l1_qpgm <- function(beta, x, y, lambda) {
    eta <- x %*% beta
    loglik <- sum(B_sublinear(y) * eta - log(factorial(y)) - D_sublinear(eta))
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
        fn = neg_loglik_l1,
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
