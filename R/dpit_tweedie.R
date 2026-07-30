#' Residuals for regression models with tweedie outcomes
#'
#' Computes DPIT residuals for Tweedie-distributed outcomes using the observed responses (\code{y}),
#' their fitted mean values (\code{mu}), the variance power parameter
#' (\eqn{\xi}), and the dispersion parameter (\eqn{\phi}).
#'
#' @usage dpit_tweedie(y, mu, xi, phi)
#' @param y Observed outcome vector.
#' @param mu Vector of fitted mean values of each outcomes.
#' @param xi Value of \eqn{\xi} such that the variance is \eqn{Var[Y] = \phi\mu^\xi}
#' @param phi Dispersion parameter \eqn{\phi}.
#'
#' @returns A `dpit` object containing DPIT residuals.
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
#' n <- 300
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
#' dpit.tweedie <- dpit_tweedie(y= y1, mu=lambda1f, xi=p.max, phi=phi1f)
#' resid.tweedie <- residuals(dpit.tweedie)
#' plot(dpit.tweedie)
#' @export
dpit_tweedie <- function(y, mu, xi, phi) {
  n <- length(y)
  p.max <- xi
  lambda1f <- mu
  phi1f <- phi
  p1f <- dtweedie(rep(0, n), mu = lambda1f, xi = p.max, phi = phi1f)
  cdf1 <- ptweedie(y, mu = lambda1f, xi = p.max, phi = phi1f)
  func <- ecdf(p1f)
  newp <- cdf1 * func(cdf1)
  .new_dpit(newp, method = "Tweedie")
}
