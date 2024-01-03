inv.nb <- function(s, mu.hat, sizef) {
  qres <- qnbinom(s, size = sizef, mu = mu.hat) - 1
  pres <- pnbinom(qres, size = sizef, mu = mu.hat)
  return(pres)
}


resid.nb <- function(model) {
  # fitted.values
  y <- model$model[, 1]
  lambda1f <- fitted.values(model)
  size1f <- summary(model)$theta
  n <- length(y)
  res <- pnbinom(y, mu = lambda1f, size = size1f)

  # residuals
  pres <- sapply(res, inv.nb, mu.hat = lambda1f, sizef = size1f)
  diag(pres) <- 0
  empcdf <- apply(pres, 2, sum) / (n - 1)
  return(empcdf)
}
