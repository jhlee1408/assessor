#' Residuals for a tobit model
#'
#' Computes DPIT residuals for tobit regression models using the observed
#' responses (`y`) and their corresponding fitted distributional parameters (`mu`, `sd`).
#'
#' @usage dpit_tobit(y, mu, sd, plot=TRUE, scale="normal", line_args=list(), ...)
#' @param y An observed outcome vector.
#' @param mu A vector of fitted mean values of latent variables.
#' @param sd A standard deviation of latent variables.
#' @param plot A logical value indicating whether or not to return QQ-plot
#' @param scale You can choose the scale of the residuals among `normal` and `uniform`.
#' The sample quantiles of the residuals are plotted against
#' the theoretical quantiles of a standard normal distribution under the normal scale,
#' and against the theoretical quantiles of a uniform (0,1) distribution under the uniform scale.
#'  The default scale is `normal`.
#' @param line_args A named list of graphical parameters passed to
#'   \code{graphics::abline()} to modify the reference (red) 45° line
#'   in the QQ plot. If left empty, a default red dashed line is drawn.
#' @param ... Additional graphical arguments passed to
#'   \code{stats::qqplot()} for customizing the QQ plot (e.g., \code{pch},
#'   \code{col}, \code{cex}, \code{xlab}, \code{ylab}).
#' @returns DPIT residuals.
#'
#' @details
#' For formulation details on semicontinuous outcomes, see \code{\link{dpit}}.
#'
#' @examples
#' ## Tobit regression model
#' library(VGAM)
#' n <- 500
#' beta13 <- 1
#' beta14 <- -3
#' beta15 <- 3
#'
#' set.seed(1234)
#' x11 <- runif(n)
#' x12 <- runif(n)
#' lambda1 <- beta13 + beta14 * x11 + beta15 * x12
#' sd0 <- 0.3
#' yun <- rnorm(n, mean = lambda1, sd = sd0)
#' y <- ifelse(yun >= 0, yun, 0)
#'
#' # Using VGAM package
#' # True model
#' fit1 <- vglm(formula = y ~ x11 + x12,
#'              tobit(Upper = Inf, Lower = 0, lmu = "identitylink"))
#' # Missing covariate
#' fit1miss <- vglm(formula = y ~ x11,
#'                  tobit(Upper = Inf, Lower = 0, lmu = "identitylink"))
#'
#' resid.tobit1 <- dpit_tobit(y = y, mu = VGAM::fitted(fit1), sd = sd0)
#' resid.tobit2 <- dpit_tobit(y = y, mu = VGAM::fitted(fit1miss), sd = sd0)
#'
#' # Using AER package
#' library(AER)
#' # True model
#' fit2 <- tobit(y ~ x11 + x12, left = 0, right = Inf, dist = "gaussian")
#' # Missing covariate
#' fit2miss <- tobit(y ~ x11, left = 0, right = Inf, dist = "gaussian")
#'
#' resid.aer1 <- dpit_tobit(y = y, mu = fitted(fit2), sd = sd0)
#' resid.aer2 <- dpit_tobit(y = y, mu = fitted(fit2miss), sd = sd0)
#' @export
dpit_tobit <- function(y, mu, sd,
                       plot=TRUE, scale="normal", line_args=list(), ...){
  if(!is.numeric(y)) y <- as.numeric(y)
  p1f <- pnorm(0, mean = mu, sd = sd)
  cdf1 <- pnorm(y, mean = mu, sd = sd)
  Fhat <- stats::ecdf(p1f)
  newp <- as.vector(cdf1 * Fhat(cdf1))
  return(newp)
}

#' @rawNamespace S3method(dpit,vglm)
dpit.vglm <- function(model,
                      plot = TRUE,
                      scale = "normal",
                      line_args = list(),
                      ...) {
  y <- model@y
  fitted <- VGAM::fitted(model)
  sd <- do.call(model@misc$link[2],args=list(theta=coef(model)[2],inverse = TRUE))
  newp <- dpit_tobit(y=y, mu=fitted, sd=sd)
  res_u <- newp
  .dpit_finalize(res_u, plot = plot, scale = scale, line_args = line_args, ...)
}


#' @rawNamespace S3method(dpit,tobit)
dpit.tobit <- function(model,
                       plot = TRUE,
                       scale = "normal",
                       line_args = list(),
                       ...) {
  y <- as.numeric(model$y)
  fitted <- VGAM::fitted(model)
  sd <- summary(model)$scale
  newp <- dpit_tobit(y=y, mu = fitted, sd=sd)
  res_u <- newp
  .dpit_finalize(res_u, plot = plot, scale = scale, line_args = line_args, ...)
}
