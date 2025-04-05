QPGM.path.neighborhood <-
  function(X, Y, nlams = 10, startb = 0, lambda = NULL) {
    n <- nrow(X)
    p <- ncol(X)
    
    # Add quadratic terms: X² and X_j X_k (optional)
    X_sq <- X^2
    X_interact <- matrix(0, n, p*(p-1)/2)
    colnames(X_interact) <- rep("", p*(p-1)/2)  # Avoid naming issues
    idx <- 1
    for (j in 1:(p-1)) {
      for (k in (j+1):p) {
        X_interact[, idx] <- X[, j] * X[, k]
        idx <- idx + 1
      }
    }
    Xq <- cbind(X, X_sq, X_interact)  # Combined design matrix
    
    if (is.null(lambda)) {
      lmax <- lambdaMax(Xq)
      lambda <- exp(seq(log(lmax), log(0.001 * lmax), length.out = nlams))
    }

    thr <- 1e-6
    maxit <- 1e6
    Xt <- cbind(1, Xq)  # Add intercept
    L <- max(eigen(t(Xq) %*% Xq, only.values = TRUE)$values)  # Lipschitz constant
    
    alphas <- 0
    Bmat <- matrix(0, ncol(Xq), nlams)  # Now includes quadratic terms
    
    if (sum(startb) == 0) {
      Bhat <- matrix(runif(ncol(Xt)), ncol(Xt), 1)
      Bhat[1] <- 0  # Intercept
    } else {
      Bhat <- startb
    }
    
    for (i in 1:nlams) {
      iter <- 1
      ind <- 1
      while (thr < ind & iter < maxit) {
        oldb <- Bhat
        # Gradient: Xt'(Y - exp(Xt * Bhat))
        grad <- t(Xt) %*% (Y - exp(Xt %*% Bhat))
        tmp <- Bhat - grad / L
        # Soft-thresholding (L1 penalty)
        Bhat <- sign(tmp) * pmax(abs(tmp) - lambda[i] / L, 0)
        Bhat[1] <- tmp[1]  # No penalty on intercept
        ind <- sum((Bhat - oldb)^2)
        iter <- iter + 1
      }
      alphas[i] <- -Bhat[1]
      Bmat[, i] <- -Bhat[2:ncol(Xt), drop = FALSE]  # Store all coefficients
    }
    
    return(list(alpha = alphas, Bmat = Bmat, lambda = lambda))
  }