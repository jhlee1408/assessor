#' @keywords internal
inv.zpois <- function(s, pzero, meanpoisson) {
  qres <- ifelse(s < (pzero + (1 - pzero) * (ppois(0, lambda = meanpoisson))), 0,
    qpois(pmax((s - pzero)/(1 - pzero), 0), lambda = meanpoisson)
  ) - 1
  pres <- ifelse(qres == -1, 0, (pzero + (1 - pzero) * (ppois(qres, meanpoisson))))
  return(pres)
}

#' @keywords internal
resid.zpois <- function(model) {
  # fitted.values
  y <- model$model[, 1]
  meanpoisson <- predict(model, type = "count")
  pzero <- predict(model, type = "zero")
  n <- length(y)
  res <- (pzero + (1 - pzero) * (ppois(y, meanpoisson)))

  # resid
  empcdf <- rep(NA,n)
  for(i in 1:n){
    pres <- inv.zpois(res[i], pzero, meanpoisson)
    pres[i] <- 0
    empcdf[i] <- sum(pres)/(n-1)
  }

  return(empcdf)
}
