QPGM <-
function(X, method = "QPGM", R = max(X), stability = "bootstrap", N = 100, beta = 0.05, 
         lmin = 0.00001, lambda.path = NULL, nlams = 20, ncores = 4, parallel = FALSE, sth = 0.8) {
  
  if (is.null(lambda.path)) {
    lmax = lambdaMax(t(X))  # Update if X includes quadratic terms
    lambda.path = exp(seq(log(lmax), log(lmin), length.out = nlams))
  }
  
  # (Rest of stability selection logic remains identical)
  
  wrapper <- function(i) {
    index = sample(1:ncol(X), b, replace = replaceF)
    ghat.path$raw = QPGM.network(X[, index], R, nlams = length(lambda.path), 
                                lambda = lambda.path, parallel = parallel, ncores = ncores)
    return(ghat.path)
  }
  
  # (Post-processing unchanged)
  class(ghat) = "GMS"
  return(ghat)
}