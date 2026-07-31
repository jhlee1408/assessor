#' Residuals for regression models with two-part outcomes
#'
#' Calculates DPIT residuals with model for semi-continuous outcomes.
#' For each component, `dpit_2pm` accepts either a fitted model or supplied
#' probabilities: exactly one of `model0` and `part0`, and exactly one of
#' `model1` and `part1`.
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
#' Exactly one of `model0` and `part0`, and exactly one of `model1` and
#' `part1`, must be supplied. Model and probability inputs may be mixed.
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
  if (!is.numeric(y) || any(!is.finite(y)) || any(y < 0)) {
    stop("y must contain finite nonnegative values.", call. = FALSE)
  }

  has_model0 <- !missing(model0)
  has_model1 <- !missing(model1)
  has_part0 <- !missing(part0)
  has_part1 <- !missing(part1)

  if (has_model0 == has_part0 || has_model1 == has_part1) {
    stop(
      paste(
        "Supply exactly one of model0 or part0, and exactly one of",
        "model1 or part1."
      ),
      call. = FALSE
    )
  }

  n <- length(y)
  positive <- y > 0

  if (has_model0) {
    if (stats::family(model0)$family != "binomial") {
      stop("model0 has to be a binomial regression.", call. = FALSE)
    }
    part0 <- stats::fitted.values(model0)
  }

  if (length(part0) != n) {
    stop("Length of part0 has to be the same as the length of y.", call. = FALSE)
  }
  if (!is.numeric(part0) || any(!is.finite(part0)) || any(part0 < 0) || any(part0 > 1)) {
    stop("Values of part0 have to be finite and between 0 and 1.", call. = FALSE)
  }

  if (has_model1) {
    if (stats::family(model1)$family != "Gamma") {
      stop("model1 has to be a gamma regression.", call. = FALSE)
    }
    if (length(stats::fitted.values(model1)) != sum(positive)) {
      stop(
        paste(
          "Length of the fitted values of model1 has to be the same as the",
          "number of positive values in y."
        ),
        call. = FALSE
      )
    }
    dispersion <- MASS::gamma.dispersion(model1)
    part1 <- stats::pgamma(
      y[positive],
      scale = stats::fitted.values(model1) * dispersion,
      shape = 1 / dispersion
    )
  }

  if (length(part1) != sum(positive)) {
    stop(
      "Length of part1 has to be the same as the number of positive values in y.",
      call. = FALSE
    )
  }
  if (!is.numeric(part1) || any(!is.finite(part1)) || any(part1 < 0) || any(part1 > 1)) {
    stop("Values of part1 have to be finite and between 0 and 1.", call. = FALSE)
  }

  cdf1 <- numeric(n)
  cdf1[!positive] <- part0[!positive]
  cdf1[positive] <- part0[positive] + (1 - part0[positive]) * part1
  newp <- cdf1 * stats::ecdf(part0)(cdf1)

  out <- .new_dpit(newp, method = "Two-part")
  model_calls <- list()
  if (has_model0) {
    model_calls$model0 <- stats::getCall(model0)
  }
  if (has_model1) {
    model_calls$model1 <- stats::getCall(model1)
  }
  if (length(model_calls) > 0L) {
    out$call <- model_calls
  }
  out
}
