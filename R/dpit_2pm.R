#' Residuals for regression models with two-part outcomes
#'
#' Calculates DPIT residuals with model for semi-continuous outcomes.
#' `dpit_2pm` can be used either with `model0` and `model1` or with `part0` and `part1` as arguments.
#'
#' @usage dpit_2pm(model0, model1, y, part0, part1)
#'
#'
#' @param model0 Model object for 0 outcomes (e.g., logistic regression)
#' @param model1 Model object for the continuous part (gamma regression)
#' @param y Semicontinuous outcomes.
#' @param part0 Alternative argument to `model0`. One can supply the sequence of probabilities \eqn{P(Y_i=0),~i=1,\ldots,n}.
#' @param part1 Alternative argument to `model1`. One can fit a regression model on the positive data and supply their probability integral transform. Note that the length of `part1` is the number of positive values in `y` and can be shorter than `part0`.
#'
#'
#' @details
#' For formulation details on semicontinuous outcomes, see \code{\link{dpit}}.
#'
#' In two-part models, the probability of zero can be modeled using a logistic regression, `model0`,
#' while the positive observations can be modeled using a gamma regression, `model1.`
#' Users can choose to use different models and supply the resulting probabilities of zero and probability integral transforms.
#'  `part0` should be the sequence of fitted probabilities of zeros \eqn{\hat{p}_0(\mathbf{X}_i) ,~i=1,\ldots,n}.
#'  `part1` should be the probability integral transform of the positive part \eqn{\hat{G}(Y_i|\mathbf{X}_i)}.
#'  Note that the length of `part1` is the number of positive values in `y` and can be shorter than `part0`.
#' Use `residuals()`, `summary()`, and `plot()` on the returned object to select
#' the residual scale, summarize the values, and draw the QQ plot.
#'
#'
#' @returns A `dpit` object containing DPIT residuals.
#'
#' @importFrom stats ecdf
#' @importFrom MASS gamma.dispersion
#' @export
#'
#' @examples
#' library(MASS)
#' n <- 500
#' beta10 <- 1
#' beta11 <- -2
#' beta12 <- -1
#' beta13 <- -1
#' beta14 <- -1
#' beta15 <- -2
#' x11 <- rnorm(n)
#' x12 <- rbinom(n, size = 1, prob = 0.4)
#'
#' p1 <- 1 / (1 + exp(-(beta10 + x11 * beta11 + x12 * beta12)))
#' lambda1 <- exp(beta13 + beta14 * x11 + beta15 * x12)
#' y2 <- rgamma(n, scale = lambda1 / 2, shape = 2)
#' y <- rep(0, n)
#' u <- runif(n, 0, 1)
#' ind1 <- which(u >= p1)
#' y[ind1] <- y2[ind1]
#'
#' # models as input
#' mgamma <- glm(y[ind1] ~ x11[ind1] + x12[ind1], family = Gamma(link = "log"))
#' m10 <- glm(y == 0 ~ x12 + x11, family = binomial(link = "logit"))
#' dpit.model <- dpit_2pm(model0 = m10, model1 = mgamma, y = y)
#' resid.model <- residuals(dpit.model, scale = "normal")
#' summary(dpit.model, scale = "normal")
#' plot(dpit.model, scale = "normal")
#'
#' # PIT as input
#' cdfgamma <- pgamma(y[ind1],
#'   scale = mgamma$fitted.values * gamma.dispersion(mgamma),
#'   shape = 1 / gamma.dispersion(mgamma)
#' )
#' p1f <- m10$fitted.values
#' dpit.pit <- dpit_2pm(y = y, part0 = p1f, part1 = cdfgamma)
#' resid.pit <- residuals(dpit.pit, scale = "uniform")
#' summary(dpit.pit, scale = "uniform")
#' plot(dpit.pit, scale = "uniform")

dpit_2pm <- function(model0, model1, y, part0, part1) {
  if (missing(y)) stop("argument y is missing, with no default")
  if (sum(!(y >= 0)) != 0) stop("y has to be nonnegative")


  if (!missing(model0) && !missing(model1) && !missing(y)) {
    if (model0$family[[1]] != "binomial") stop("model0 has to be a logistic regression")
    if (model1$family[[1]] != "Gamma") stop("model1 has to be a gamma regression")
    if (length(model1$fitted.values) != sum(y > 0)) stop("Length of the fitted values of model1 has to be the same as the number of positive values in y")
    if (length(model0$fitted.values) != length(y)) {
      stop("Length of the fitted values of model0 has to be the same as the length of y")
    } else {
      n <- length(y)
      cdfgamma <- pgamma(y[y > 0],
        scale = model1$fitted.values * gamma.dispersion(model1),
        shape = 1 / gamma.dispersion(model1)
      )
      p1f <- model0$fitted.values
      cdf1 <- rep(0, n)
      cdf1[y == 0] <- model0$fitted.values[y == 0]
      cdf1[y > 0] <- model0$fitted.values[which(y > 0)] + (1 - model0$fitted.values[which(y > 0)]) * cdfgamma
      newp <- cdf1 * ecdf(p1f)(cdf1)
    }
  }

  if (!missing(part0) && !missing(part1) && !missing(y)) {
    if (length(part0) != length(y)) stop("Length of part0 has to be the same as the length of y")
    if (sum(y != 0) != length(part1)) stop("Length of part1 has to be the same as the number of positive values in y")
    if (sum((part1 < 0) + (part1 > 1)) != 0) stop("Values of part0 and part1 have to be between 0 and 1")
    if (sum((part0 < 0) + (part0 > 1)) != 0) stop("Values of part0 and part1 have to be between 0 and 1")
    n <- length(y)
    cdf1 <- rep(0, n)
    cdf1[y == 0] <- part0[y == 0]
    cdf1[y > 0] <- part0[y > 0] + (1 - part0[y > 0]) * part1
    newp <- cdf1 * ecdf(part0)(cdf1)
  }

  if (!missing(model0) && !missing(part1) && !missing(y)) {
    if (length(model0$fitted.values) != length(y)) stop("Length of the fitted values of model0 has to be the same as the length of y")
    if (sum(y != 0) != length(part1)) stop("Length of part1 has to be the same as the number of positive values in y")
    if (model0$family[[1]] != "binomial") stop("model0 has to be a logistic regression")
    if (sum((part1 < 0) + (part1 > 1)) != 0) stop("Values of part0 and part1 have to be between 0 and 1")
    n <- length(y)
    cdf1 <- rep(0, n)
    part0 <- model0$fitted.values
    cdf1[y == 0] <- part0[y == 0]
    cdf1[y > 0] <- part0[y > 0] + (1 - part0[y > 0]) * part1
    newp <- cdf1 * ecdf(part0)(cdf1)
  }

  if (!missing(part0) && !missing(model1) && !missing(y)) {
    if (length(model1$fitted.values) != sum(y > 0)) stop("Length of the fitted values of model1 has to be the same as the number of positive values in y")
    if (length(part0) != length(y)) stop("Length of part0 has to be the same as the length of y")
    if (sum((part0 < 0) + (part0 > 1)) != 0) stop("Values of part0 and part1 have to be between 0 and 1")
    if (model1$family[[1]] != "Gamma") stop("model1 has to be a gamma regression")
    n <- length(y)
    cdf1 <- rep(0, n)
    cdfgamma <- pgamma(y[y > 0],
      scale = model1$fitted.values * gamma.dispersion(model1),
      shape = 1 / gamma.dispersion(model1)
    )
    cdf1[y == 0] <- part0[y == 0]
    cdf1[y > 0] <- part0[y > 0] + (1 - part0[y > 0]) * cdfgamma
    newp <- cdf1 * ecdf(part0)(cdf1)
  }

  out <- .new_dpit(newp, method = "Two-part")
  model_calls <- list()
  if (!missing(model0)) {
    model_calls$model0 <- stats::getCall(model0)
  }
  if (!missing(model1)) {
    model_calls$model1 <- stats::getCall(model1)
  }
  if (length(model_calls) > 0L) {
    out$call <- model_calls
  }
  out
}
