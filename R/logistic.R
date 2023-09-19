inv.bin0 <- function(s, q10){
  qres <- 1*(s>=q10)*(s<1)*0 + 1*(s<q10)*(s<1)*(-1) + 1*(s==1)
  pres <- ifelse(qres==0, q10, ifelse(qres==1,1,0))
  return(pres)
}


inv.ordi1 <- function(s,q10){
  qses <- (1*(s==0)*(-2)+1*(s>=q10)*(s<1)*0 + 1*(s<q10)*(s<1)*(-1)+1*(s==1))+1
  pses <- ifelse(qses==0,q10,ifelse(qses==1,1,0))
  return(pses)
}

resid.logi <- function(model){
  # fitted.values
  k <- length(model$lev)
  out <- as.numeric(factor(model$model[,1],ordered = TRUE))
  n <- length(out)

  fitprob <- fitted(model)
  q <- t(apply(fitprob,1, cumsum))
  inde <- cbind(1:n,out)
  res <- q[inde]

  res <- matrix(rep(res,n),n,n,byrow=TRUE)
  note <- list()

  for(i in 1:k){
    note[[i]] <- fitprob[,i]*(res > q[,i])
  }
  pres <- Reduce("+", note)

  diag(pres) <- 0
  empcdf <- colSums(pres)/(n-1)

  ses <- ifelse(out==k,q[,1],0)
  ses <- matrix(rep(ses,n),n,n,byrow=TRUE)
  pses <- ifelse(ses==0,0,ifelse(ses<q[,1],q[,k-1],1))

  diag(pses) <- 0
  rempcdf <- colSums(pses)/(n-1)
  finalecdf <- c(empcdf[out<k],rempcdf[out==k])

  return(finalecdf)
}
