## ===== Function to fit Original Poisson Graphical Model =====
fit_Original_Poisson_Graphical_Model <- function(X, lambda) {
  #' Fit Original Poisson Graphical Model
  #' 
  #' @param X a n-by-p matrix, where n is the number of samples, p is the number of nodes (variables)
  #' @param lambda L1-regularization coefficient
  #' @return p-by-p adjacency matrix indicating the edeges (1 means edge, 0 means no edge)
  p <- ncol(X)
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
    adjacency[s, -s] <- (abs(coefs) > 1e-6*p) * 1
  }
  # make adjacency symmetric
  adjacency <- adjacency + t(adjacency)
  adjacency[adjacency!=0] = 1
  return(adjacency)
}