inv.zpois <- function(s, pzero, meanpoisson){
  qres <- ifelse(s < (pzero + (1 - pzero) * (ppois(0, lambda = meanpoisson))), 0,
                 qpois(pmax((s - pzero) / (1 - pzero), 0), lambda = meanpoisson))-1
  pres <- ifelse(qres==-1,0,(pzero + (1 - pzero) * (ppois(qres, meanpoisson))))
  return(pres)
}


resid.zpois <- function(model){
  # fitted.values
  y <- model$model[,1]
  meanpoisson <- predict(model, type="count")
  pzero <- predict(model, type="zero")
  n <- length(y)
  res <-(pzero + (1 - pzero) * (ppois(y, meanpoisson)))

  # residuals
  pres <- sapply(res, inv.zpois, pzero, meanpoisson)
  diag(pres) <-0
  empcdf <- apply(pres, 2,sum)/(n-1)
  return(empcdf)
}
