#' Residuals for regression models with negative binomial outcomes
#'
#' Calculates Double Probability Integral Transform (DPIT) residuals for
#' regression model whose response follows a negative binomial distribution.
#' Unlike `resid_disc`, which infers the fitted CDF from a model object,
#' `dpit_nb` requires only the vector of fitted means and a dispersion parameter, so it can be applied
#' to GLMs, GAMs, or other custom estimators—as long as fitted values are available.
#'
#' @usage dpit_nb(fitted, y=NULL, size=NULL, plot = TRUE, scale="uniform")
#'
#' @param fitted Numeric vector of fitted means.
#' @param y Outcome variable.
#' @param size Dispersion (size) parameter for the negative binomial distribution; e.g., a `glm.nb` fitted with a negative‐binomial family, you can extract it with `summary(model)$theta`.
#' @param plot Logical; if `TRUE` (default) a QQ-plot of the residuals is displayed.
#' @param scale You can choose the scale of the residuals among normal and uniform scales. The sample quantiles of the residuals are plotted against the theoretical quantiles of a standard normal distribution under the normal scale, and against the theoretical quantiles of a uniform (0,1) distribution under the uniform scale. The default scale is normal
#'
#' @returns DPIT residuals and their mean CRPS. If `plot=TRUE`, also produces a QQ plot.
#'
#' @import stats
#' @import graphics
#' @importFrom scoringRules crps_unif
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
dpit_nb <- function(fitted, y = NULL, size = NULL, plot=TRUE, scale="uniform") {
  if(is.null(size)) stop("size has to be specified")
  if(is.null(y)) stop("y is not given.")
  lambda1f <- fitted
  size1f <- size
  n <- length(fitted)
  res <- pnbinom(y, mu = lambda1f, size = size1f)

  empcdf <- rep(NA,n)
  for(i in 1:n){
    qres <- qnbinom(res[i], mu=lambda1f, size=size1f)-1
    pres <- pnbinom(qres,mu=lambda1f,size=size1f)
    pres[i] <- 0
    empcdf[i] <-sum(pres)/(n-1)
  }
  res <- empcdf
  if (plot == TRUE) {
    if (scale == "normal") {
      empcdf <- qnorm(empcdf)
      n <- length(empcdf)
      qqplot(qnorm(ppoints(n)), (empcdf),
             main = "QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
             cex.lab = 1, cex.axis = 1, cex.main = 1.5, lwd = 1.5
      )
      abline(0, 1, col = "red", lty = 5, cex.lab = 2, cex.axis = 2, cex.main = 2, lwd = 1.5)
    }
    if (scale == "uniform") {
      n <- length(empcdf)
      qqplot(ppoints(n), empcdf,
             main = "QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
             cex.lab = 1, cex.axis = 1, cex.main = 1.5, lwd = 1.5
      )
      abline(0, 1, col = "red", lty = 5, cex.lab = 2, cex.axis = 2, cex.main = 2, lwd = 1.5)
    }
  } else {
    if (scale == "normal") empcdf <- qnorm(empcdf)
    if (scale == "uniform") empcdf <- empcdf
  }
  return(list("DPIT"= empcdf,
              "CRPS"= mean(crps_unif(res))))
}
