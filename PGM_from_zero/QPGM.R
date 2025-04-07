library(pheatmap)
library(XMRF)
library(pROC)
library(verification)

set.seed(1231)
# simulate data
sim <- XMRF.Sim(n = 200, p = 50, model = "TPGM", graph.type = "hub")
simDat <- sim$X
p <- nrow(simDat)
n <- ncol(simDat)


# Node-conditional distribution via GLM: P(Z) = exp(theta B(Z) - C(Z) - D(theta))
B <- function(x) {
  return(x) # Canonical form of sufficient statistic
}

C <- function(x) {
  return(x^2) # Quadratic base measure
}

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

# Perform the M-estimation
learn_edge <- function(data, sufficient_statistic, log_normalization, lambda, threshold = 0.1) {
  p = nrow(data)
  THETA_prob <- matrix(0, nrow=p, ncol=p)
  for (s in 1:100){
    cat(paste("QPGM: learning edges ... in progress: ", floor(100*(s/100)), "%", collapse=""),"\r")
    flush.console()
    THETA <- NULL
    for (i in 1:p){
      X_minus_s <- simDat[-i,]
      X_s <- simDat[i,]
      
      optimization_target <- function(
    par, X_s, X_minus_s, sufficient_statistic, log_normalization, lambda
      ){
        n <- ncol(X_minus_s)
        log_likelihood <- mean(sapply(
          1:n, 
          function(m) {
            k <- t(par) %*% sufficient_statistic(X_minus_s[,m])
            return(-sufficient_statistic(X_s[m])*k + log_normalization(k))
          }
        ))
        l1_regularization <- lambda * sum(abs(par))
        return(log_likelihood + l1_regularization)
      }
      
      fit <- optim(par = rnorm(p-1, 0, 5),
                   optimization_target,
                   X_s = X_s,
                   X_minus_s = X_minus_s,
                   sufficient_statistic = B,
                   log_normalization = D,
                   lambda = lambda,
                   method = "Nelder-Mead"
      )
      theta <- fit$par
      theta <- append(theta, 1, after = i-1)
      THETA <- rbind(THETA, theta)
    }
    THETA_prob <- THETA_prob + (THETA>0.01)
  }
  return(THETA_prob/100)
}

D <- get_D(sufficient_statistic = B, base_measure = C)
THETA = learn_edge(simDat, B, C, lambda = 1)

roc_curve <- roc(as.numeric(sim$B), as.numeric(THETA))
with(
  roc_curve, 
  plot(
    1-specificities, 
    sensitivities, 
    type = 'l',
    main = "ROC for QPGM",
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    xlim = c(0,1), ylim = c(0,1),
    lwd = 2
  )
)
axis(1, axTicks(1), labels=F)
lines(c(0, 1), c(0, 1), col="grey")
