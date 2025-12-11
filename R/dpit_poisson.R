#' Residuals for regression models with Poisson outcomes
#'
#' Calculates Double Probability Integral Transform (DPIT) residuals for
#' regression model whose response follows a Poisson distribution.
#' Unlike `resid_disc`, which infers the fitted CDF from a model object,
#' `dpit_pois` requires only the vector of fitted means, so it can be applied
#' to GLMs, GAMs, or other custom estimators—as long as
#' Poisson fitted values are available.
#'
#' @usage dpit_pois(fitted, y =NULL, plot = TRUE, scale="uniform", line_args=list(), ...)
#'
#' @param fitted Numeric vector of fitted means.
#' @param y Outcome variable.
#' @param plot Logical; if `TRUE` (default) a QQ-plot of the residuals is displayed.
#' @param scale You can choose the scale of the residuals among normal and uniform scales. The sample quantiles of the residuals are plotted against the theoretical quantiles of a standard normal distribution under the normal scale, and against the theoretical quantiles of a uniform (0,1) distribution under the uniform scale. The default scale is normal.
#' @param line_args A named list of graphical parameters passed to
#'   \code{graphics::abline()} to modify the reference (red) 45° line
#'   in the QQ plot. If left empty, a default red dashed line is drawn.
#' @param ... Additional graphical arguments passed to
#'   \code{stats::qqplot()} for customizing the QQ plot (e.g., \code{pch},
#'   \code{col}, \code{cex}, \code{xlab}, \code{ylab}).
#'
#' @returns DPIT residuals. If `plot=TRUE`, also produces a QQ plot.
#'
#' @import stats
#' @import graphics
#' @export
#'
#' @references Yang, Lu. "Double Probability Integral Transform Residuals for Regression Models with Discrete Outcomes." arXiv preprint arXiv:2308.15596 (2023).
#' @examples
#' library(assessor)
#' library(mgcv)
#' set.seed(1234)
#' n  <- 500
#' x1 <- runif(n, 0, 1)
#' x2 <- runif(n, 0, 1)
#' f1 <- function(x)  2 * sin(2*pi*x)
#' f2 <- function(x) -1.5 * (x - .5)^2
#' eta <- -0.5 + f1(x1) + f2(x2)
#' mu  <- exp(eta)
#' theta_true <- 1.3
#' y <- rnbinom(n, size = theta_true, mu = mu)
#' dat <- data.frame(y, x1, x2)
#' gam_nb <- gam(y ~ s(x1) + s(x2),
#'               family = nb(),
#'               data   = dat,
#'               method = "REML")
#' gam_pois <- gam(y ~ s(x1) + s(x2),
#'                 family = poisson(),
#'                 data = dat,
#'                 method = "REML")
#' dpit_nb(gam_nb$fitted.values, y=y, size=gam_nb$family$getTheta(TRUE))
#' dpit_pois(gam_pois$fitted.values, y=y)
dpit_pois <- function(fitted, y=NULL, plot=TRUE, scale="uniform", line_args=list(), ...) {
  lambda1f <- fitted
  n <- length(y)
  res <- ppois(y, lambda = lambda1f)
  empcdf <- rep(NA, n)
  for (i in 1:n) {
    qres <- qpois(res[i], lambda = lambda1f) - 1
    pres <- ppois(qres, lambda = lambda1f)
    pres[i] <- 0
    empcdf[i] <- sum(pres)/(n - 1)
  }
  if (plot == TRUE) {
    qqplot.resid(empcdf, scale, line_args, ...)
  } else {
    if (scale == "normal") empcdf <- qnorm(empcdf)
    if (scale == "uniform") empcdf <- empcdf
  }
  return(empcdf)
}
