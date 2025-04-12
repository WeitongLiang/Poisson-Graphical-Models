message("Loading Quatratic Poisson Graphical Model Family...")

# ===== Helper functions for Quadratic PGM =====
N = 100 # boundary
theta <- seq(from = -19, to = 19, length.out = 1001) # theta (or eta in traditional GLM)

## ===== log normalizing constant =====
D <- function(th) {
  xs <- 0:N
  log_terms <- th*xs - xs^2
  max_log <- max(log_terms)
  log(sum(pmax(exp(log_terms - max_log), .Machine$double.eps))) + max_log
}
D <- Vectorize(D)

## ===== p.d.f. of node conditional distribution =====
pdf <- function(x, theta) {exp(theta*x - x^2 - D(theta))}

## ===== g^{-1}: inverse of link function =====
# g_inv <- function(th) {
#   sum(sapply(0:N, function(x){
#     sapply(0:N, function(y){
#       (x)*exp((x+y)*th - (x^2+y^2))
#     })
#   }))
# }
g_inv <- function(eta) {
  sapply(eta, function(t) {
    t <- pmin(pmax(t, -19), 19)
    xs <- 0:N
    terms <- xs * exp(t * xs - xs^2 - D(t))
    val <- sum(terms, na.rm = TRUE)
    ifelse(is.finite(val), val, 1e-8)
  })
}

## ===== g: link function =====
# mu_vals <- sapply(theta, g_inv)
mu_vals <- g_inv(theta)
ord <- order(mu_vals)
if(any(duplicated(mu_vals))) {stop("g duplicates")}
g <- approxfun(x = mu_vals[ord], y = theta[ord], rule = 2)

## ===== V: Variance function V(mu) =====
V_vals <- sapply(theta, function(th) {
  v <- sum(sapply(0:N, function(x) x^2 * exp(th * x - x^2 - D(th)))) - (g_inv(th))^2
  ifelse(is.finite(v) && v > 1e-6, v, 1e-6)
})
V <- approxfun(x = mu_vals[ord], y = V_vals[ord], rule = 2)
# V <- function(mu) {
#   vapply(mu, function(m) {
#     th <- g(m)
#     xs <- 0:N
#     log_terms <- th*xs - xs^2 - D(th)
#     max_log <- max(log_terms)
#     probs <- exp(log_terms - max_log)/sum(exp(log_terms - max_log))
#     var <- sum(xs^2 * probs) - (sum(xs * probs))^2
#     ifelse(is.finite(var) & var > 0, var, 1e-8)
#   }, numeric(1))
# }

## ===== Deviance.residual =====
dev.resid <- function(y, mu, wt) {
  # print(paste("mu:", paste(mu[1:5], collapse = ",")))
  # print(paste("g(y):", paste(g(y)[1:5], collapse = ",")))
  # print(paste("g(mu):", paste(g(mu)[1:5], collapse = ",")))
  # print(paste("D(g(y)):", paste(D(g(y)[1:5]), collapse = ",")))
  # print(paste("D(g(mu)):", paste(D(g(mu)[1:5]), collapse = ",")))
  if (any(is.null(mu)) || any(is.na(mu))) {
    stop("mu is NULL or NA in dev.resid")
  }
  2 * wt * ((g(y)-g(mu))*y-D(g(y))+D(g(mu)))
}

## ===== AIC function =====
# aic <- function(y, n, mu, wt, dev) {
#   -2 * sum(
#     sapply(1:n, function(i){log(pdf(y[i], g(mu[i])))})
#   )
# }
aic <- function(y, n, mu, wt, dev) {
  -2 * sum(
    sapply(seq_along(y), function(i) {
      log(pdf(y[i], g(mu[i])))
    })
  )
}

## ===== d(mu)/d(eta) =====
# mu.eta <- function(th) {
#   sum(sapply(0:N, function(x){
#     sapply(0:N, function(y){
#       (x*(x+y))*exp((x+y)*th - (x^2+y^2))
#     })
#   }))
# }

# mu.eta <- function(eta) {
#   vapply(eta, function(t) {
#     t <- pmin(pmax(t, -20), 20)
#     xs <- 0:N
#     
#     # Calculate in log space
#     log_terms <- t*xs - xs^2 - D(t)
#     if (any(!is.finite(log_terms))) return(1e-8)
#     
#     max_log <- max(log_terms)
#     rel_probs <- exp(log_terms - max_log)
#     probs <- rel_probs/sum(rel_probs)
#     
#     # Calculate variance
#     mean_val <- sum(xs * probs)
#     second_moment <- sum(xs^2 * probs)
#     variance <- second_moment - mean_val^2
#     
#     ifelse(is.finite(variance) & variance > 0, variance, 1e-8)
#   }, numeric(1))
# }
mu.eta <- function(eta) {
  h <- 1e-4
  val <- (g_inv(eta + h) - g_inv(eta - h)) / (2 * h)
  ifelse(is.finite(val) & val > 1e-8, val, 1e-8)
}

mu.eta <- Vectorize(mu.eta)

## ===== valid mu =====
validmu <- function(mu) {
  all(is.finite(mu)) && all(mu >= g_inv(-19)) && all(mu <= g_inv(19))
}

## ===== valid eta =====
valideta <- function(eta) {
  all(is.finite(eta)) && all(eta >= -19) && all(eta <= 19)
}

## ===== family for Quadratic PGM =====
family_QPGM <- function() {
  # initialize <- expression({
  #   if (any(y < 0)) stop("Negative values not allowed")
  #   n <- rep.int(1, nobs)
  #   mustart <- pmin(pmax(y + 0.1, 0.1), g_inv(20))
  # })
  # initialize <- expression({
  #   if (any(y < 0)) stop("Negative values not allowed")
  #   # n <- rep.int(1, nobs)
  #   # y_adj <- pmin(pmax(y, 0.1), N - 0.1)
  #   # mustart <- y_adj
  #   mustart <- pmax(y, 0.1)
  # })
  initialize <- expression({
    if (any(y < 0)) stop("Negative values not allowed")
    n <- rep.int(1, nobs)
    mustart <- pmax(y, 0.1)  # Use mustart instead of mu here
  })
  
  structure(list(
    family = "QPGM",
    link = "custom",
    linkfun = Vectorize(g),
    linkinv = Vectorize(g_inv),
    variance = V,
    dev.resids = dev.resid,
    mu.eta = mu.eta,  # Already properly vectorized
    initialize = initialize,
    validmu = validmu,
    valideta = valideta,
    aic = aic
  ), class = "family")
}

message("Success!")

## ===== Function to fit Quadratic Poisson Graphical Model =====
fit_Quadratic_Poisson_Graphical_Model <- function(X, lambda) {
  #' Fit Quadratic Poisson Graphical Model
  #' 
  #' @param X a n-by-p matrix, where n is the number of samples, p is the number of nodes (variables)
  #' @param lambda L1-regularization coefficient
  #' @return p-by-p adjacency matrix indicating the edeges (1 means edge, 0 means no edge)
  p <- ncol(X)
  adjacency <- matrix(0, p, p)
  
  for (s in 1:p) {
    if (length(unique(X[, s])) < 2) {
      stop("Colinearity occors!")
    }
    fit <- glmnet(
      x = X[, -s, drop = FALSE],  # All samples, other variables
      y = X[, s],                 # Current variable
      family = family_QPGM(),
      alpha = 1,
      lambda = lambda
      # standardize = TRUE
    )
    coefs <- as.vector(coef(fit))[-1]  # Remove intercept
    adjacency[s, -s] <- (abs(coefs) > 1e-6) * 1
  }
  # make adjacency symmetric
  adjacency <- adjacency + t(adjacency)
  adjacency[adjacency!=0] = 1
  return(adjacency)
}

