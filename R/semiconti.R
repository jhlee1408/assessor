#' Residuals for regression models with semicontinuous outcomes
#'
#' Calculates the DPIT residuals for regression models with semi-continuous outcomes.
#' The semi-continuous regression model such as
#' a Tweedie regression model from `tweedie` package or a Tobit regression model
#' from `VGAM`, `AER` packages is used in this function.
#'
#' @usage resid_semiconti(model, plot=TRUE, scale = "normal", line_args=list(), ...)
#' @seealso [resid_2pm()]
#'
#' @param model Model object (e.g., `tweedie`, `vglm`, and `tobit`)
#' @param plot A logical value indicating whether or not to return QQ-plot
#' @param scale You can choose the scale of the residuals between `normal` and `uniform` scales. The default scale is `normal`.
#' @param line_args A named list of graphical parameters passed to
#'   \code{graphics::abline()} to modify the reference (red) 45° line
#'   in the QQ plot. If left empty, a default red dashed line is drawn.
#' @param ... Additional graphical arguments passed to
#'   \code{stats::qqplot()} for customizing the QQ plot (e.g., \code{pch},
#'   \code{col}, \code{cex}, \code{xlab}, \code{ylab}).

#'
#' @returns Residuals. If plot=TRUE, also produces a QQ plot.
#'
#' @details
#' The DPIT residual for the \eqn{i}th semicontinuous observation is defined as follows:
#' \deqn{\hat{r}_i = \frac{\hat{F}(Y_i|X_i)}{n}\sum_{j=1}^{n}I\bigg(\hat{p}_0(X_j) \leq \hat{F}(Y_i|X_i)\bigg),}
#' which has a null distribution of uniformity.
#' \eqn{\hat{F}} refers to the fitted cumulative distribution function,
#' and \eqn{\hat{p}_0} refers to the fitted probability of being zero.
#'
#'
#' @references Lu Yang (2024). Diagnostics for Regression Models with Semicontinuous Outcomes, Biometrics, https://arxiv.org/abs/2401.06347
#'
#' @importFrom stats ecdf
#' @import tweedie
#' @importFrom tweedie dtweedie
#' @importFrom tweedie ptweedie
#' @importFrom VGAM fitted
#' @export
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
#' resid.tweedie <- resid_semiconti(model1)
#'
#' ## Tobit regression model
#' library(VGAM)
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
#' fit1 <- vglm(formula = y ~ x11 + x12, tobit(Upper = Inf, Lower = 0, lmu = "identitylink"))
#' # Missing covariate
#' fit1miss <- vglm(formula = y ~ x11, tobit(Upper = Inf, Lower = 0, lmu = "identitylink"))
#'
#' resid.tobit1 <- resid_semiconti(fit1, plot = TRUE)
#' resid.tobit2 <- resid_semiconti(fit1miss, plot = TRUE)
#'
#' # Using AER package
#' library(AER)
#' # True model
#' fit2 <- tobit(y ~ x11 + x12, left = 0, right = Inf, dist = "gaussian")
#' # Missing covariate
#' fit2miss <- tobit(y ~ x11, left = 0, right = Inf, dist = "gaussian")
#' resid.aer1 <- resid_semiconti(fit2, plot = TRUE)
#' resid.aer2 <- resid_semiconti(fit2miss, plot = TRUE)
resid_semiconti <- function(model, plot = TRUE, scale = "normal", line_args=list(), ...) {
  if (!(scale %in% c("normal", "uniform"))) stop("scale has to be either normal or uniform")
  is.vglm <- isS4(model)
  if (is.vglm) {
    if (!(paste(model@call)[1] %in% c("vglm"))) stop("model has to be tweedie, vglm or tobit")
  }
  if (!is.vglm) {
    if (!is.null(model$family)) {
      model.family <- model$family[[1]]
    } # Tweedie
    else {
      model.family <- "AER"
    } # AER
    if (!(paste(model$call)[1] %in% c("tobit")) & model.family != "Tweedie") stop("model has to be tweedie, vglm or tobit")
  }

  if (!is.vglm && model.family == "Tweedie") {
    y1 <- model$y
    p.max <- get("p", envir = environment(model$family$variance))
    n <- length(y1)
    lambda1f <- model$fitted.values
    phi1f <- summary(model)$dis
    p1f <- dtweedie(rep(0, n), mu = lambda1f, xi = p.max, phi = phi1f)
    cdf1 <- ptweedie(y1, mu = lambda1f, xi = p.max, phi = phi1f)

    func <- ecdf(p1f)

    newp <- cdf1 * func(cdf1)
  }

  if (!is.vglm && model.family == "AER") {
    p1f <- pnorm(0, mean = VGAM::fitted(model), sd = summary(model)$scale)
    cdf1 <- pnorm(y, mean = VGAM::fitted(model), sd = summary(model)$scale)
    newp <- cdf1 * ecdf(p1f)(cdf1)
  }

  if (is.vglm) {
    y <- model@y
    p1f <- pnorm(0, mean = fitted(model), sd = exp(coef(model)[2]))
    cdf1 <- pnorm(y, mean = fitted(model), sd = exp(coef(model)[2]))
    newp <- cdf1 * ecdf(p1f)(cdf1)
    newp <- as.vector(newp)
  }

  if (plot == T) {
    qqplot.resid(newp, scale, line_args, ...)
  } else {
    if (scale == "normal") newp <- qnorm(newp)
    if (scale == "uniform") newp <- newp
  }
  return(newp)
}
