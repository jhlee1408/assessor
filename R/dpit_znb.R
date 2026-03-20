#' @keywords internal
inv.znb <- function(s, pzero, mu.hat, size1f) {
  qres <- ifelse(s < (pzero + (1 - pzero) * (pnbinom(0, size = size1f, mu = mu.hat))), 0,
                 qnbinom(pmax((s - pzero) / (1 - pzero), 0), mu = mu.hat, size = size1f)
  ) - 1
  pres <- ifelse(qres == -1, 0, (pzero + (1 - pzero) * (pnbinom(qres, size = size1f, mu = mu.hat))))
  return(pres)
}

#' Residuals for regression models with zero-inflated negative binomial outcomes
#'
#' Computes DPIT residuals for regression models with zero-inflated negative
#' binomial outcomes using the observed counts (`y`) and their fitted distributional
#' parameters (`mu`, `pzero`, `size`).
#'
#' @usage dpit_znb(y, mu, pzero, size, plot=TRUE, scale="normal", line_args=list(), ...)
#' @param y An observed outcome vector.
#' @param mu A vector of fitted mean values for the count (non-zero) component.
#' @param pzero A vector of fitted probabilities for the zero-inflation component.
#' @param size A dispersion parameter of the negative binomial distribution.
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
#'
#' @details
#' For formulation details on discrete outcomes, see \code{\link{dpit}}.
#'
#' @examples
#' ## Zero-Inflated Negative Binomial
#' library(pscl)
#' n <- 500
#' set.seed(1234)
#'
#' # Covariates
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.7)
#'
#' # Coefficients
#' beta0 <- -2
#' beta1 <-  2
#' beta2 <-  1
#' beta00 <- -2
#' beta10 <-  2
#'
#' # NB dispersion (size = theta; larger => closer to Poisson)
#' theta_true <- 1.2
#'
#' # Mean of NB count part
#' mu_true <- exp(beta0 + beta1 * x1 + beta2 * x2)
#'
#' # Excess zero probability (logit)
#' p0 <- 1 / (1 + exp(-(beta00 + beta10 * x1)))
#'
#' ## simulate outcomes
#' z  <- rbinom(n, size = 1, prob = 1 - p0)              # 1 => from NB, 0 => structural zero
#' y1 <- rnbinom(n, size = theta_true, mu = mu_true)     # NB count draw
#' y  <- ifelse(z == 0, 0, y1)
#'
#' ## True model
#' modelzero1 <- zeroinfl(y ~ x1 + x2 | x1, dist = "negbin", link = "logit")
#' y1 <- modelzero1$y
#' mu1    <- stats::predict(modelzero1, type = "count")
#' pzero1 <- stats::predict(modelzero1, type = "zero")
#' theta1 <- modelzero1$theta
#' resid.zero1 <- dpit_znb(y = y1, pzero = pzero1, mu = mu1, size = theta1)
#'
#' ## Ignoring zero-inflation: NB only
#' modelzero2 <- MASS::glm.nb(y ~ x1 + x2)
#' y2 <- modelzero2$y
#' mu2    <- fitted(modelzero2)
#' theta2 <- modelzero2$theta
#' resid.zero2 <- dpit_nb(y = y2, mu = mu2, size = theta2)
#' @export
dpit_znb <- function(y, mu, pzero, size, plot=TRUE, scale="normal", line_args=list(), ...) {
  mu.hat <- mu
  pzero <- pzero
  size1f <- size
  n <- length(y)
  res <- (pzero + (1 - pzero) * (pnbinom(y, size = size1f, mu = mu.hat)))

  # residuals
  empcdf <- rep(NA,0)
  for(i in 1:n){
    pres <- inv.znb(res[i], pzero, mu.hat, size1f)
    pres[i] <-0
    empcdf[i] <- sum(pres)/(n-1)
  }
  .dpit_finalize(empcdf, plot=plot, scale=scale, line_args=line_args, ...)
}
