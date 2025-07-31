#' @keywords internal
resid.pois <- function(model) {
  # fitted.values
  y <- model$model[, 1]
  lambda1f <- fitted.values(model)
  n <- length(y)
  res <- ppois(y, lambda = lambda1f)

  empcdf <- rep(NA,n)
  for(i in 1:n){
    qres <- qpois(res[i], lambda=lambda1f)-1
    pres <- ppois(qres,lambda = lambda1f)
    pres[i] <- 0
    empcdf[i] <-sum(pres)/(n-1)
  }
  return(empcdf)
}
