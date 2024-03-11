#' @keywords internal
resid.bin <- function(model) {
  # fitted.values
  y <- model$y
  n <- length(y)
  q10 <- 1 - model$fitted.values
  res <- ifelse(y == 0, q10, 1)
  ses <- ifelse(y == 1, q10, 0)

  # residuals
  empcdf <- rep(NA,n)
  for(i in 1:n){
    if(i %in% which(y==1)) next
    qres <- 1*(res[i] >= q10)*(res[i] < 1)*0 + 1 * (res[i]<q10) * (res[i]<1)*(-1)+1*(res[i]==1)
    pres <- ifelse(qres ==0, q10, ifelse(qres == 1, 1, 0))
    pres[i] <- 0
    empcdf[i] <- sum(pres)/ (n-1)
  }

  for(i in 1:n){
    if(i %in% which(y==0)) next
    qses <- (1*(ses[i]==0)*(-2)+1*(ses[i]>=q10)*(ses[i]<1)*0+1*(ses[i] < q10)*(ses[i]<1)*(-1)+1*(ses[i]==1))+1
    pses <- ifelse(qses ==0, q10, ifelse(qses == 1, 1, 0))
    pses[i] <- 0
    empcdf[i] <- sum(pses)/ (n-1)
  }
  return(empcdf)
}

