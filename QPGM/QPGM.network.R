QPGM.network <-
  function(X, R = max(X), nlams = 10, lambda = NULL, parallel = TRUE, ncores = ncores) {
    ghat <- array(0, dim = c(nrow(X), nrow(X), nlams))

    if (is.null(lambda)) {
      lmax <- lambdaMax(t(X)) # Recompute for QPGM if X includes quadratic terms
      lambda <- exp(seq(log(lmax), log(0.01 * lmax), length.out = nlams))
    }

    wrapper <- function(i) {
      fit <- QPGM.path.neighborhood(t(X[-i, ]), X[i, ], nlams = nlams, lambda = lambda)
      fit$beta <- as.matrix(fit$Bmat)
      # (Rest of the indexing logic remains the same)
      return(ghat[i, , ])
    }

    # (Parallel/posting processing unchanged)
    return(ghat)
  }
