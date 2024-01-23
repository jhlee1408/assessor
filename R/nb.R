
resid.nb <- function(model){
  # fitted.values
  y <- model$model[, 1]
  lambda1f <- fitted.values(model)
  size1f <- summary(model)$theta
  n <- length(y)
  res <- pnbinom(y, mu = lambda1f, size = size1f)

  empcdf <- rep(NA,n)
  for(i in 1:n){
    qres <- qnbinom(res[i], mu=lambda1f, size=size1f)-1
    pres <- pnbinom(qres,mu=lambda1f,size=size1f)
    pres[i] <- 0
    empcdf[i] <-sum(pres)/(n-1)
  }
  return(empcdf)
}
