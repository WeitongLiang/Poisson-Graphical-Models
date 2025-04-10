quadratic_poisson_family <- function(link = "log", y_max = 100) {
  # Set up link function
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  stats <- make.link(linktemp)
  
  # Function to compute normalizing constant Z(λ)
  compute_z <- function(lambda) {
    ys <- 0:y_max
    sum(lambda^ys / exp(ys^2))
  }
  
  # Function to compute mean μ(λ) = E[y]
  compute_mean <- function(lambda) {
    ys <- 0:y_max
    probs <- lambda^ys / exp(ys^2)
    sum(ys * probs) / compute_z(lambda)
  }
  
  # Function to compute variance Var(y)
  compute_var <- function(lambda) {
    ys <- 0:y_max
    probs <- lambda^ys / exp(ys^2)
    z <- compute_z(lambda)
    mu <- sum(ys * probs) / z
    sum((ys - mu)^2 * probs) / z
  }
  
  # Create lookup tables for faster computation
  lambda_grid <- exp(seq(-10, 10, length.out = 1001))
  mean_grid <- sapply(lambda_grid, compute_mean)
  var_grid <- sapply(lambda_grid, compute_var)
  
  # Create approximating functions
  approx_mean <- approxfun(lambda_grid, mean_grid)
  approx_var <- approxfun(lambda_grid, var_grid)
  
  # Variance function for the family
  variance <- function(mu) {
    # Find lambda that gives this mean
    lambda_est <- suppressWarnings(
      uniroot(function(x) approx_mean(x) - mu, 
              range(lambda_grid))$root
    )
    approx_var(lambda_est)
  }
  
  # Deviance residuals
  dev.resids <- function(y, mu, wt) {
    # Find lambda for each mu
    lambdas <- sapply(mu, function(m) {
      suppressWarnings(
        uniroot(function(x) approx_mean(x) - m, 
                range(lambda_grid))$root
      )
    })
    
    # Compute log-probabilities
    log_p_mu <- function(y, lambda) {
      y * log(lambda) - y^2 - log(compute_z(lambda))
    }
    
    # Saturated model (mu = y)
    log_p_sat <- sapply(seq_along(y), function(i) {
      if (y[i] == 0) 0 else log_p_mu(y[i], lambdas[i])
    })
    
    2 * wt * (log_p_sat - log_p_mu(y, lambdas))
  }
  
  # Initialize
  initialize <- expression({
    if (any(y < 0)) stop("y must be non-negative")
    n <- rep(1, nobs)
    mustart <- y + 0.1
  })
  
  # AIC calculation
  aic <- function(y, n, mu, wt, dev) {
    -2 * sum(log(dcustom(y, mu))) + 2 * length(n)
  }
  
  # Custom density function (for AIC)
  dcustom <- function(y, mu) {
    lambdas <- sapply(mu, function(m) {
      suppressWarnings(
        uniroot(function(x) approx_mean(x) - m, 
                range(lambda_grid))$root
      )
    })
    sapply(seq_along(y), function(i) {
      (lambdas[i]^y[i] / exp(y[i]^2)) / compute_z(lambdas[i])
    })
  }
  
  structure(list(
    family = "custom_discrete_exp",
    link = linktemp,
    linkfun = stats$linkfun,
    linkinv = stats$linkinv,
    variance = variance,
    dev.resids = dev.resids,
    mu.eta = stats$mu.eta,
    initialize = initialize,
    validmu = function(mu) all(mu > 0),
    valideta = stats$valideta,
    aic = aic,
    dcustom = dcustom
  ), class = "family")
}