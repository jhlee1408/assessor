inv.bin0 <- function(s, q10){
  qres <- 1*(s>=q10)*(s<1)*0 + 1*(s<q10)*(s<1)*(-1) + 1*(s==1)
  pres <- ifelse(qres==0, q10, ifelse(qres==1,1,0))
  return(pres)
}


inv.bin1 <- function(s,q10){
  qses <- (1*(s==0)*(-2) + 1*(s>=q10)*(s<1)*0 + 1*(s<q10)*(s<1)*(-1)+1*(s==1))+1
  pses <- ifelse(qses==0,q10,ifelse(qses==1,1,0))
  return(pses)
}


resid.bin <- function(model){
  # fitted.values
  y <- model$y
  n <- length(y)
  q10 <- 1- model$fitted.values
  res <- ifelse(y==0, q10, 1)
  ses <- ifelse(y==1, q10, 0)
  # residuals
  pres <- sapply(res, inv.bin0, q10)
  pses <- sapply(ses, inv.bin1, q10)

  diag(pres) <-0; diag(pses) <-0
  empcdf <- apply(pres, 2,sum)/(n-1)
  rempcdf <- apply(pses ,2 ,sum)/(n-1)

  fin.empcdf <- rep(NA, n)
  fin.empcdf[y==0] <- empcdf[y==0]
  fin.empcdf[y==1] <- empcdf[y==1]
  return(fin.empcdf)
}
