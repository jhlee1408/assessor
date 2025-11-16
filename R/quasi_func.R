#' @keywords internal
listvec <- function(x) {
  x[1]:x[2]
}

## Binary
## \hat{U}. y is the outcome, and q0 is the fitted probability of zero
margin01 <- function(u, y, q0, h) {
  wei <- 1 * ((q0 - u)^2 < 5 * h^2) * (1 - ((q0 - u)^2) / h^2 / 5)
  l <- sum(wei * 1 * (y == 0)) / sum(wei)
  l
}
margin01 <- Vectorize(margin01, "u")

## bandwidth selector
bandwidth01 <- function(y, q0) {
  bw <- npregbw(ydat = 1 * (y == 0), xdat = q0, ckertype = "epanechnikov")
  return(bw$bw)
}

#' @keywords internal
#' @importFrom utils modifyList
resid.bin_quasi <- function(model, line_args, ...){
  y <- model$y
  q10 <- 1 - model$fitted.values
  h <- bandwidth01(y = y, q0 = q10)
  x.input <- seq(0,1, length.out=101)

  qq_defaults <- list(
    type='l',
    main = "Quasi, Binary",
    ylab = expression(hat(U) * "(s)"),
    xlab = "s",
    cex.lab = 2, cex.axis = 2, cex.main = 2,
    lwd = 2,
    xlim = c(min(pbinom(0, size=1,prob=q10)),
             max(pbinom(0, size=1,prob=q10)))
  )
  qq_args <- utils::modifyList(qq_defaults, list(...))
  do.call(base::plot,
          c(list(x = x.input,
                 y = margin01(x.input, y=y, q0=q10, h=h)
            ),qq_args))
  default_abline <- list(
    a   = 0,
    b   = 1,
    col = "red",
    lty = 5,
    lwd = 2
  )
  ab_args <- modifyList(default_abline, line_args)
  do.call(abline, ab_args)
}


## Poisson
marginesti <- function(u, y, lambdaf, h) {
  n <- length(y)
  # F^{-1}(u)
  inv1 <- qpois(u, lambdaf)
  inv1m <- ifelse(inv1 > 0, inv1 - 1, 0)
  n1 <- cbind(inv1, inv1m)
  # H^+ and H^_
  p1 <- ppois(n1, lambdaf)
  # find out which one is closer to u
  ind1 <- apply(abs(p1 - u), 1, which.min)
  # weight function
  wei <- 1 * ((p1[cbind(1:n, ind1)] - u)^2 < 5 * h^2) *
    (1 - ((p1[cbind(1:n, ind1)] - u)^2) / h^2 / 5)
  l <- sum(wei * 1 * (y <= n1[cbind(1:n, ind1)])) / sum(wei)
  l
}
marginesti <- Vectorize(marginesti, "u")

####### bandwidth selector
bandwidthp <- function(y, lambdaf){
  n <- length(y)
  newout <- unlist(sapply(split(cbind(qpois(0.1, lambdaf), qpois(0.9, lambdaf)), 1:n), listvec))
  newlambda <- rep(lambdaf, times = qpois(0.9, lambdaf) - qpois(0.1, lambdaf) + 1)
  newx <- ppois(newout, newlambda)
  newy <- 1 * (rep(y, times = qpois(0.9, lambdaf) - qpois(0.1, lambdaf) + 1) <= newout)
  bws <- npregbw(ydat = newy[which(newx <= 0.9)], xdat = newx[which(newx <= 0.9)], ckertype = "epanechnikov")
  return(bws$bw)
}
#' @keywords internal
resid.pois_quasi <- function(model, line_args, ...){
  y <- model$y
  lambda1f <- model$fitted.values
  h <- bandwidthp(y = y, lambdaf = lambda1f)
  x.input <- seq(0,1, length.out= 101)

  defaults <- list(
    type='l',
    main = "Quasi, Poisson",
    ylab = expression(hat(U) * "(s)"),
    xlab = "s",
    cex.lab = 2, cex.axis = 2, cex.main = 2,
    lwd = 2,
    xlim = c(min(ppois(0, lambda = lambda1f)), 1)
  )
  qq_args <- utils::modifyList(defaults, list(...))
  do.call(base::plot,
          c(list(x = x.input,
                 y = marginesti(u=x.input, y=y, lambdaf= lambda1f, h=h)),
            qq_args))
  default_abline <- list(
    a   = 0,
    b   = 1,
    col = "red",
    lty = 5,
    lwd = 2
  )
  ab_args <- modifyList(default_abline, line_args)
  do.call(abline, ab_args)
}

############
##   NB   ##
############
marginnb <- function(u, y, lambdaf, sizef,h) {
  n <- length(y)
  inv1 <- qnbinom(u, mu = lambdaf, size = sizef)
  inv1m <- ifelse(inv1 > 0, inv1 - 1, 0)
  n1 <- cbind(inv1, inv1m)
  p1 <- pnbinom(n1, mu = lambdaf, size = sizef)
  ind1 <- apply(abs(p1 - u), 1, which.min)
  wei <- 1 * ((p1[cbind(1:n, ind1)] - u)^2 < 5 * h^2) *
    (1 - ((p1[cbind(1:n, ind1)] - u)^2) / h^2 / 5)
  l <- sum(wei * 1 * (y <= n1[cbind(1:n, ind1)])) / sum(wei)
  l
}

marginnb <- Vectorize(marginnb, "u")

## bandwidth selector
bandwidthnb <- function(y, lambdaf, sizef) {
  n <- length(y)
  newout <- unlist(sapply(split(cbind(qnbinom(0.1, mu = lambdaf, size = sizef), qnbinom(0.9, mu = lambdaf, size = sizef)), 1:n), listvec))
  newlambda <- rep(lambdaf, times = qnbinom(0.9, mu = lambdaf, size = sizef) - qnbinom(0.1, mu = lambdaf, size = sizef) + 1)
  newx <- pnbinom(newout, mu = newlambda, size = sizef)
  newy <- 1 * (rep(y, times = qnbinom(0.9, mu = lambdaf, size = sizef) - qnbinom(0.1, mu = lambdaf, size = sizef) + 1) <= newout)

  newy1 <- newy[which(newx <= 0.9)]
  newx1 <- newx[which(newx <= 0.9)]
  bws <- npregbw(ydat = newy1, xdat = newx1, ckertype = "epanechnikov")
  return(bws$bw)
}

#' @keywords internal
resid.nb_quasi <- function(model, line_args, ...){
  y <- model$y
  lambda1f <- model$fitted.values
  size1f <- summary(model)$theta
  h <- bandwidthnb(y = y, lambdaf = lambda1f, sizef = size1f)
  x.input <- seq(0,1, length.out=101)
  defaults <- list(
    type= 'l',
    main = "Quasi, NB",
    ylab = expression(hat(U) * "(s)"),
    xlab = "s",
    cex.lab = 2, cex.axis = 2, cex.main = 2,
    lwd = 2,
    xlim = c(min(pnbinom(0, size=size1f, mu=lambda1f)), 1)
    )
  qq_args <- utils::modifyList(defaults, list(...))
  do.call(base::plot,
          c(list(
            x = x.input,
            y = marginnb(u=x.input, y=y, lambdaf=lambda1f, sizef=size1f, h=h)
            ), qq_args))
  default_abline <- list(
    a   = 0,
    b   = 1,
    col = "red",
    lty = 5,
    lwd = 2
  )
  ab_args <- modifyList(default_abline, line_args)
  do.call(abline, ab_args)
}

###########
## Ord  ###
###########
### \hat{U} y is the outcome, q0=P(y<=1), q1=P(y<=2)
marginm <- function(x, y, p1, h) { # p1 = t(apply(model$fitted.values,1 ,cumsum))[,-n]
  n <- length(y)
  ind1 <- apply(abs(p1 - x), 1, which.min)
  wei <- 1 * ((p1[cbind(1:n, ind1)] - x)^2 < 5 * h^2) *
    (1 - ((p1[cbind(1:n, ind1)] - x)^2) / h^2 / 5)
  l <- sum(wei * 1 * (y <= (ind1 - 1))) / sum(wei)
  l
}

marginm <- Vectorize(marginm, "x")

### Bandwidth selection
bandwidthord <- function(y, p1) {
  k <- max(y)
  n <- length(y)
  xdat. <- rep(NA, n*k)
  for(i in 0:(k-1)) xdat.[(i*n+1):((i+1)*n)] <- p1[, (i+1)]
  ydat. <- rep(NA,n*k)
  for(i in 0:(k-1)) ydat.[(i*n+1):((i+1)*n)] <- 1*(y <= i)
  bw <- npregbw(ydat = ydat., xdat = xdat., ckertype = "epanechnikov")
  return(bw$bw)
}

#' @keywords internal
resid.ordi_quasi <- function(model, line_args, ...){
  y <- as.numeric(model$model[,1])-1
  p1 <- t(apply(model$fitted.values,1 ,cumsum))
  h <- bandwidthord(y = y, p1=p1)
  x.input <- seq(0,1, length.out=101)
  defaults <- list(
    type='l',
    main = "Quasi, Ordinal",
    ylab = expression(hat(U) * "(s)"),
    xlab = "s",
    cex.lab = 2, cex.axis = 2, cex.main = 2,
    lwd = 2,
    xlim = c(0,1))
  qq_args <- utils::modifyList(defaults, list(...))
  do.call(base::plot,
          c(list(x = x.input,
                 y = marginm(x=x.input, y=y, p1=p1, h=h),
                 qq_args)))
  default_abline <- list(
            a   = 0,
            b   = 1,
            col = "red",
            lty = 5,
            lwd = 2
          )
  ab_args <- modifyList(default_abline, line_args)
  do.call(abline, ab_args)
}

###########
## zpois ##
###########
marginzerop <- function(u, y, pzero, meanpoisson,h) {
  n <- length(y)
  inv1 <- ifelse(u < (pzero + (1 - pzero) * (ppois(0, lambda = meanpoisson))), 0,
                 qpois(pmax((u - pzero) / (1 - pzero), 0), lambda = meanpoisson)
  )
  inv1m <- ifelse(inv1 > 0, inv1 - 1, 0)
  n1 <- cbind(inv1, inv1m)
  p1 <- cbind(
    (pzero + (1 - pzero) * (ppois(inv1, meanpoisson))),
    (pzero + (1 - pzero) * (ppois(inv1m, meanpoisson)))
  )
  ind1 <- apply(abs(p1 - u), 1, which.min)
  wei <- 1 * ((p1[cbind(1:n, ind1)] - u)^2 < 5 * h^2) *
    (1 - ((p1[cbind(1:n, ind1)] - u)^2) / h^2 / 5)
  l <- sum(wei * 1 * (y <= n1[cbind(1:n, ind1)])) / sum(wei)
  l
}
marginzerop <- Vectorize(marginzerop, "u")

### Bandwidth selection
bandwidth0p <- function(y, pzero, meanpoisson) {
  n <- length(y)
  newout <- unlist(sapply(split(cbind(qpois(0.1, meanpoisson), qpois(0.9, meanpoisson)), 1:n), listvec))
  newlambda <- rep(meanpoisson, times = qpois(0.9, meanpoisson) - qpois(0.1, meanpoisson) + 1)
  newzero <- rep(pzero, times = qpois(0.9, meanpoisson) - qpois(0.1, meanpoisson) + 1)
  newx <- newzero + (1 - newzero) * ppois(newout, newlambda)
  newy <- 1 * (rep(y, times = qpois(0.9, meanpoisson) - qpois(0.1, meanpoisson) + 1) <= newout)
  return(npregbw(ydat = newy[which(newx <= 0.9 & newx >= 0.1)], xdat = newx[which(newx <= 0.9 & newx > 0.1)], ckertype = "epanechnikov")$bw)
}

#' @keywords internal
resid.zpois_quasi <- function(model, line_args,...){
  y <- model$model[, 1]
  meanpoisson <- predict(model, type = "count")
  pzero <- predict(model, type = "zero")
  h <- bandwidth0p(y = y, pzero = pzero, meanpoisson = meanpoisson)

  x.input <- seq(0,1, length.out=101)

  defaults <- list(
    type='l',
    main = "Quasi, 0-Inflated Poisson",
    ylab = expression(hat(U) * "(s)"),
    xlab = "s",
    cex.lab = 2, cex.axis = 2, cex.main = 2,
    lwd = 2,
    xlim = c(min(pzero + (1 - pzero) * (ppois(0, meanpoisson))), 1)
  )
  qq_args <- utils::modifyList(defaults, list(...))
  do.call(base::plot,
          c(list(x = x.input,
                 y = marginzerop(u=x.input, y=y, pzero = pzero, meanpoisson = meanpoisson, h=h)
                 ),qq_args))
  default_abline <- list(
            a   = 0,
            b   = 1,
            col = "red",
            lty = 5,
            lwd = 2
          )
  ab_args <- modifyList(default_abline, line_args)
  do.call(abline, ab_args)
}


###########
##  znb  ##
###########
marginzerop.nb <- function(u, y, pzero, size1f, mu.hat, h) {
  n <- length(y)
  inv1 <- ifelse(u < (pzero + (1 - pzero) * (pnbinom(0, size = size1f, mu = mu.hat))), 0,
                 qnbinom(pmax((u - pzero) / (1 - pzero), 0), mu = mu.hat, size = size1f)
  )
  inv1m <- ifelse(inv1 > 0, inv1 - 1, 0)
  n1 <- cbind(inv1, inv1m)
  p1 <- cbind(
    (pzero + (1 - pzero) * (pnbinom(inv1, size=size1f, mu=mu.hat))),
    (pzero + (1 - pzero) * (pnbinom(inv1m, size=size1f, mu=mu.hat)))
  )
  ind1 <- apply(abs(p1 - u), 1, which.min)
  wei <- 1 * ((p1[cbind(1:n, ind1)] - u)^2 < 5 * h^2) *
    (1 - ((p1[cbind(1:n, ind1)] - u)^2) / h^2 / 5)
  l <- sum(wei * 1 * (y <= n1[cbind(1:n, ind1)])) / sum(wei)
  l
}

marginzerop.nb <- Vectorize(marginzerop.nb, "u")

### Bandwidth selection
bandwidth0p.nb <- function(y, pzero, mu.hat, size1f) {
  n <- length(y)
  newout <- unlist(sapply(split(cbind(qnbinom(0.1, mu=mu.hat, size=size1f), qnbinom(0.9, mu=mu.hat, size=size1f)), 1:n), listvec))
  newlambda <- rep(mu.hat, times = qnbinom(0.9, mu=mu.hat, size=size1f) - qnbinom(0.1, mu=mu.hat, size=size1f) + 1)
  newzero <- rep(pzero, times = qnbinom(0.9, mu=mu.hat, size=size1f) - qnbinom(0.1, mu=mu.hat, size=size1f) + 1)
  newx <- newzero + (1 - newzero) * ppois(newout, newlambda)
  newy <- 1 * (rep(y, times = qnbinom(0.9, mu=mu.hat, size=size1f) - qnbinom(0.1, mu=mu.hat, size=size1f) + 1) <= newout)

  return(npregbw(ydat = newy[which(newx <= 0.9 & newx >= 0.1)], xdat = newx[which(newx <= 0.9 & newx > 0.1)], ckertype = "epanechnikov")$bw)
}

#' @keywords internal
resid.znb_quasi <- function(model, line_args, ...){
  y <- model$model[, 1]
  mu.hat <- predict(model, type = "count")
  size1f <- model$theta
  pzero <- predict(model, type = "zero")
  h <- bandwidth0p.nb(y = y, pzero = pzero, mu.hat = mu.hat, size1f= size1f)
  x.input <- seq(0,1, length.out=101)

  defaults <- list(
    type='l',
    main = "Quasi, 0-Inflated NB",
    ylab = expression(hat(U) * "(s)"),
    xlab = "s",
    cex.lab = 2, cex.axis = 2, cex.main = 2,
    lwd = 2,
    xlim = c(min(pzero + (1 - pzero) * (pnbinom(0, size=size1f,mu=mu.hat))), 1)
  )
  qq_args <- utils::modifyList(defaults, list(...))
  do.call(base::plot,
          c(list(x = x.input,
                 y = marginzerop.nb(u=x.input, y=y, mu.hat =mu.hat,pzero=pzero, size1f = size1f, h=h)
                 ),qq_args))
  default_abline <- list(
            a   = 0,
            b   = 1,
            col = "red",
            lty = 5,
            lwd = 2
          )
  ab_args <- modifyList(default_abline, line_args)
  do.call(abline, ab_args)
}



