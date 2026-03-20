#' @keywords internal
inv.zpois <- function(s, pzero, meanpoisson) {
  qres <- ifelse(s < (pzero + (1 - pzero) * (ppois(0, lambda = meanpoisson))), 0,
    qpois(pmax((s - pzero)/(1 - pzero), 0), lambda = meanpoisson)
  ) - 1
  pres <- ifelse(qres == -1, 0, (pzero + (1 - pzero) * (ppois(qres, meanpoisson))))
  return(pres)
}


#' Residuals for regression models with zero-inflated Poisson outcomes
#'
#' Computes DPIT residuals for regression models with zero-inflated Poisson
#' outcomes using the observed counts(`y`) and their fitted distributional
#' parameters(`mu`, `pzero`).
#'
#'
#' @usage dpit_zpois(y, mu, pzero, plot=TRUE, scale="normal", line_args=list(), ...)
#' @param y An observed outcome vector.
#' @param mu A vector of fitted mean values for the count (non-zero) component.
#' @param pzero A vector of fitted probabilities for the zero-inflation component.
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
#' @importFrom stats family
#' @export
#'
#' @details
#' For formulation details on discrete outcomes, see \code{\link{dpit}}.
#'
#'
#' @examples
#' ## Zero-Inflated Poisson
#' library(pscl)
#' n <- 500
#' set.seed(1234)
#' # Covariates
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.7)
#' # Coefficients
#' beta0 <- -2
#' beta1 <- 2
#' beta2 <- 1
#' beta00 <- -2
#' beta10 <- 2
#'
#' # Mean of Poisson part
#' lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
#' # Excess zero probability
#' p0 <- 1 / (1 + exp(-(beta00 + beta10 * x1)))
#' ## simulate outcomes
#' y0 <- rbinom(n, size = 1, prob = 1 - p0)
#' y1 <- rpois(n, lambda1)
#' y <- ifelse(y0 == 0, 0, y1)
#' ## True model
#' modelzero1 <- zeroinfl(y ~ x1 + x2 | x1, dist = "poisson", link = "logit")
#' y1 <- modelzero1$y
#' mu1    <- stats::predict(modelzero1, type = "count")
#' pzero1 <- stats::predict(modelzero1, type = "zero")
#' resid.zero1 <- dpit_zpois(y= y1, pzero=pzero1, mu=mu1)
#'
#' ## Zero inflation
#' modelzero2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' y2 <- modelzero2$y
#' mu2    <- fitted(modelzero2)
#' resid.zero2 <- dpit_pois(y= y2, mu=mu2)
#'
#'
#' @export
dpit_zpois <- function(y, mu, pzero,
                       plot=TRUE, scale="normal", line_args=list(), ...) {
  if (missing(y) || missing(pzero) || missing(mu)) {
    stop("y, pzero, and mu are required.", call. = FALSE)
  }
  y <- as.numeric(y)
  n <- length(y)

  if (length(pzero) != n) stop("length(pzero) must equal length(y).", call. = FALSE)
  if (length(mu)    != n) stop("length(mu) must equal length(y).", call. = FALSE)
  if (any(y < 0)) stop("y must be nonnegative.", call. = FALSE)

  res <- pzero + (1 - pzero) * stats::ppois(y, lambda = mu)

  empcdf <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    pres <- inv.zpois(res[i], pzero, mu)  # inv.zpois is assumed to exist
    pres[i] <- 0
    empcdf[i] <- sum(pres) / (n - 1)
  }
  .dpit_finalize(empcdf, plot=plot, scale=scale, line_args=line_args,...)
}
