#' Residuals for regression models with two-part outcomes
#'
#' Calculates DPIT proposed residuals for model for semi-continuous outcomes.
#' `resid_2pm` can be used either with `model0` and `model1` or with `part0` and `part1` as arguments.
#'
#' @usage resid_2pm(model0, model1, y, part0, part1, plot=TRUE, scale = "normal")
#'
#' @seealso [resid_semiconti()]
#'
#' @param model0 Model object for 0 outcomes (e.g., logistic regression)
#' @param model1 Model object for the continuous part (gamma regression)
#' @param y Semicontinuous outcome variables
#' @param part0 Alternative argument to `model0`. One can supply the sequence of probabilities \eqn{P(Y_i=0),~i=1,\ldots,n}.
#' @param part1 Alternative argument to `model1`. One can fit a regression model on the positive data and supply their probability integral transform. Note that the length of `part1` is the number of positive values in `y` and can be shorter than `part0`.
#' @param plot A logical value indicating whether or not to return QQ-plot
#' @param scale You can choose the scale of the residuals among `normal` and `uniform` scales. The default scale is `normal`.
#'
#' @importFrom scoringRules crps_unif
#'
#' @details
#' The DPIT residuals for regression models with semi-continuous outcomes are \deqn{\hat{r}_i=\frac{\hat{F}(Y_i|\mathbf{X}_i)}{n}\sum_{j=1}^n1\left(\hat{p}_0(\mathbf{X}_j)\leq \hat{F}(Y_i|\mathbf{X}_i)\right), i=1,\ldots,n,}
#' where \eqn{\hat{p}_0(\mathbf{X}_i)} is the fitted probability of zero, and \eqn{\hat{F}(\cdot|\mathbf{X}_i)} is the  fitted cumulative distribution function for the \eqn{i}th observation. Furthermore, \deqn{\hat{F}(y|\mathbf{x})=\hat{p}_0(\mathbf{x})+\left(1-\hat{p}_0(\mathbf{x})\right)\hat{G}(y|\mathbf{x})}
#' where \eqn{\hat{G}} is the fitted cumulative distribution for the positive data.
#'
#' In two-part models, the probability of zero can be modeled using a logistic regression, `model0`,
#' while the positive observations can be modeled using a gamma regression, `model1.`
#' Users can choose to use different models and supply the resulting probability transforms.
#'  `part0` should be the sequence of fitted probabilities of zeros \eqn{\hat{p}_0(\mathbf{X}_i) ,~i=1,\ldots,n}.
#'  `part1` should be the probability integral transform of the positive part \eqn{\hat{G}(Y_i|\mathbf{X}_i)}.
#'  Note that the length of `part1` is the number of positive values in `y` and can be shorter than `part0`.
#'
#'
#' @returns residuals and their mean CRPS. If plot=TRUE, also produces a QQ plot.
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
#' resid.model <- resid_2pm(model0 = m10, model1 = mgamma, y = y)
#'
#' # PIT as input
#' cdfgamma <- pgamma(y[ind1],
#'   scale = mgamma$fitted.values * gamma.dispersion(mgamma),
#'   shape = 1 / gamma.dispersion(mgamma)
#' )
#' p1f <- m10$fitted.values
#' resid.pit <- resid_2pm(y = y, part0 = p1f, part1 = cdfgamma)
resid_2pm <- function(model0, model1, y, part0, part1, plot = TRUE, scale = "normal") {
  if (!(scale %in% c("normal", "uniform"))) stop("scale has to be either normal or uniform")
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

  if (plot == T) {
    if (scale == "normal") {
      newp <- qnorm(newp)
      n <- length(newp)
      qqplot(qnorm(ppoints(n)), newp[is.finite(newp)],
        main = "QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
        cex.lab = 1, cex.axis = 1, cex.main = 1.5, lwd = 1.5
      )
      abline(0, 1, col = "red", lty = 5, cex.lab = 2, cex.axis = 2, cex.main = 2, lwd = 1.5)
    }
    if (scale == "uniform") {
      n <- length(newp)
      qqplot(ppoints(n), newp[is.finite(newp)],
        main = "QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
        cex.lab = 1, cex.axis = 1, cex.main = 1.5, lwd = 1.5
      )
      abline(0, 1, col = "red", lty = 5, cex.lab = 2, cex.axis = 2, cex.main = 2, lwd = 1.5)
    }
  } else {
    if (scale == "normal") newp <- qnorm(newp)
    if (scale == "uniform") newp <- newp
  }
  return(list("residuals" = newp,
               "CRPS"= crps_unif(newp)))
}
