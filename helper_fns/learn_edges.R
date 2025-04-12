learn_edge <- function(method, data, sufficient_statistic, base_measure, lambda, M, threshold = 0.1) {
  D <- get_D(
    sufficient_statistic = sufficient_statistic,
    base_measure = base_measure
  )
  
  p = nrow(data)
  THETA_prob <- matrix(0, nrow=p, ncol=p)
  for (s in 1:M) {
    cat(paste(method, ": learning edges with data ( p =",p, "n =",ncol(data), ") in progress: ", round(100*(s/M), 2), "%", collapse=""),"\r")
    flush.console()
  
    THETA <- NULL
    for (i in 1:p){
      X_minus_s <- data[-i,]
      X_s <- data[i,]
      optimization_target <- function(
    par, X_s, X_minus_s, sufficient_statistic, log_normalization, lambda
      ){
        n <- ncol(X_minus_s)
        log_likelihood <- sapply(
          1:n, 
          function(m) {
            k <- t(par) %*% sufficient_statistic(X_minus_s[,m])
            return(-sufficient_statistic(X_s[m])*k + D(k))
          }
        )
        l1_regularization <- lambda * sum(abs(par))
        return(mean(log_likelihood) + l1_regularization)
      }
      
      fit <- optim(par = rnorm(p-1, 0, 100),
                   optimization_target,
                   X_s = X_s,
                   X_minus_s = X_minus_s,
                   sufficient_statistic = sufficient_statistic,
                   log_normalization = D,
                   lambda = lambda,
                   method = "Nelder-Mead"
      )
      theta <- fit$par
      theta <- append(theta, 1, after = i-1)
      THETA <- rbind(THETA, theta)
    }
    edge <- abs(THETA)>0.01
    for (row in 1:nrow(THETA)) {
      for (col in 1:row) {
        if (any(c(edge[row,col], edge[col,row]))) {
          edge[row, col] <- TRUE
          edge[col, row] <- TRUE
        }
      }
    }
    THETA_prob <- THETA_prob + edge
  }
  return(THETA_prob/M)
}