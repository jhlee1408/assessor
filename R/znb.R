#' @keywords internal
inv.znb <- function(s, pzero, mu.hat, size1f) {
  qres <- ifelse(s < (pzero + (1 - pzero) * (pnbinom(0, size = size1f, mu = mu.hat))), 0,
    qnbinom(pmax((s - pzero) / (1 - pzero), 0), mu = mu.hat, size = size1f)
  ) - 1
  pres <- ifelse(qres == -1, 0, (pzero + (1 - pzero) * (pnbinom(qres, size = size1f, mu = mu.hat))))
  return(pres)
}


#' @keywords internal
resid.znb <- function(model) {
  # fitted.values
  y <- model$model[, 1]
  mu.hat <- predict(model, type = "count")
  pzero <- predict(model, type = "zero")
  size1f <- model$theta
  n <- length(y)
  res <- (pzero + (1 - pzero) * (pnbinom(y, size = size1f, mu = mu.hat)))

  # residuals
  empcdf <- rep(NA,0)
  for(i in 1:n){
    pres <- inv.znb(res[i], pzero, mu.hat, size1f)
    pres[i] <-0
    empcdf[i] <- sum(pres)/(n-1)
  }
  return(empcdf)
}
