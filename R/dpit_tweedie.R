#' Residuals for regression models with tweedie outcomes
#'
#' Computes DPIT residuals for Tweedie-distributed outcomes using the observed responses (\code{y}),
#' their fitted mean values (\code{mu}), the variance power parameter
#' (\eqn{\xi}), and the dispersion parameter (\eqn{\phi}).
#'
#' @usage dpit_tweedie(y, mu, xi, phi, plot=TRUE, scale="normal", line_args=list(), ...)
#' @param y Observed outcome vector.
#' @param mu Vector of fitted mean values of each outcomes.
#' @param xi Value of \eqn{\xi} such that the variance is \eqn{Var[Y] = \phi\mu^\xi}
#' @param phi Dispersion parameter \eqn{\phi}.
#' @param plot A logical value indicating whether or not to return QQ-plot
#' The sample quantiles of the residuals are plotted against
#' @param scale You can choose the scale of the residuals among `normal` and `uniform`.
#' the theoretical quantiles of a standard normal distribution under the normal scale,
#' and against the theoretical quantiles of a uniform (0,1) distribution under the uniform scale.
#'  The default scale is `normal`.
#' @param line_args A named list of graphical parameters passed to
#'   \code{graphics::abline()} to modify the reference (red) 45° line
#'   in the QQ plot. If left empty, a default red dashed line is drawn.
#' @param ... Additional graphical arguments passed to
#'   \code{stats::qqplot()} for customizing the QQ plot (e.g., \code{pch},
#'   \code{col}, \code{cex}, \code{xlab}, \code{ylab}).
#'
#' @returns DPIT residuals.
#'
#' @details
#' For formulation details on semicontinuous outcomes, see \code{\link{dpit}}.
#' @import tweedie
#'
#'
#' @examples
#' ## Tweedie model
#' library(tweedie)
#' library(statmod)
#' n <- 500
#' x11 <- rnorm(n)
#' x12 <- rnorm(n)
#' beta0 <- 5
#' beta1 <- 1
#' beta2 <- 1
#' lambda1 <- exp(beta0 + beta1 * x11 + beta2 * x12)
#' y1 <- rtweedie(n, mu = lambda1, xi = 1.6, phi = 10)
#' # Choose parameter p
#' # True model
#' model1 <-
#'   glm(y1 ~ x11 + x12,
#'     family = tweedie(var.power = 1.6, link.power = 0)
#'   )
#' y1 <- model1$y
#' p.max <- get("p", envir = environment(model1$family$variance))
#' lambda1f <- model1$fitted.values
#' phi1f <- summary(model1)$dis
#' resid.tweedie <- dpit_tweedie(y= y1, mu=lambda1f, xi=p.max, phi=phi1f)
#' @export
dpit_tweedie <- function(y, mu, xi, phi, plot=TRUE, scale="normal", line_args=list(), ...){
  n <- length(y)
  p.max <- xi
  lambda1f <- mu
  phi1f <- phi
  p1f <- dtweedie(rep(0, n), mu = lambda1f, xi = p.max, phi = phi1f)
  cdf1 <- ptweedie(y, mu = lambda1f, xi = p.max, phi = phi1f)
  func <- ecdf(p1f)
  newp <- cdf1 * func(cdf1)
  .dpit_finalize(newp, plot=plot, scale=scale, line_args=line_args, ...)
}
