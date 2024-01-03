inv.pois <- function(s, mu.hat) {
  qres <- qpois(s, lambda = mu.hat) - 1
  pres <- ppois(qres, lambda = mu.hat)
  return(pres)
}


resid.pois <- function(model) {
  # fitted.values
  y <- model$model[, 1]
  lambda1f <- fitted.values(model)
  n <- length(y)
  res <- ppois(y, lambda = lambda1f)

  # residuals
  pres <- sapply(res, inv.pois, lambda1f)
  diag(pres) <- 0
  empcdf <- apply(pres, 2, sum) / (n - 1)
  return(empcdf)
}
