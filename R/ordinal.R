#' @keywords internal
resid.ordi <- function(model) {
  # fitted.values
  k <- length(model$lev)
  out <- as.numeric(factor(model$model[, 1], ordered = TRUE))
  n <- length(out)

  fitprob <- fitted(model)
  q <- t(apply(fitprob, 1, cumsum))
  inde <- cbind(1:n, out)
  res <- q[inde]


   # for loop without max values
  empcdf <- rep(NA, n)
  for(i in 1:n){
    if(i %in% which(out==k)) next
    note <- matrix(NA, ncol=k, nrow=n)
    for(p in 1:k){
        note[,p] <- fitprob[, p] * (res[i] > q[, p])
    }
    note.sum <- rowSums(note)
    note.sum[i] <- 0
    empcdf[i] <- sum(note.sum)/(n-1)
  }

  # for loop with max values
  ses <- ifelse(out == k, q[, 1], 0)
  for(i in 1:n){
    if(i %in% which(out != k)) next
    pses <- (ses[i] < q[,1])*q[,k-1]
    pses[pses==0] <- 1
    pses[i] <- 0
    empcdf[i] <- sum(pses)/(n-1)
  }
  return(empcdf)
}
