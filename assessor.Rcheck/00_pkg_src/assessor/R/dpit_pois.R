#' Residuals for regression models with poisson outcomes
#'
#' Computes DPIT residuals for Poisson outcomes regression using the observed counts (`y`) and their
#' corresponding fitted mean values (`mu`).
#'
#' @usage dpit_pois(y, mu, plot=TRUE, scale="normal", line_args=list(), ...)
#' @param y An observed outcome vector.
#' @param mu A vector of fitted mean values.
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
#' For formulation details on discrete outcomes, see \code{\link{dpit}}.
#'
#' @examples
#' ## Poisson example
#' n <- 500
#' set.seed(1234)
#' # Covariates
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.7)
#' # Coefficients
#' beta0 <- -2
#' beta1 <- 2
#' beta2 <- 1
#' lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
#' y <- rpois(n, lambda1)
#'
#' # True model
#' poismodel <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' y1 <- poismodel$y
#' p1f <- fitted(poismodel)
#' resid.poi <- dpit_pois(y=y1, mu=p1f)
#'
#' @export
dpit_pois <- function(y, mu, plot=TRUE, scale="normal", line_args=list(), ...) {
  n <- length(y)
  lambda1f <- mu
  res <- ppois(y, lambda = lambda1f)
  empcdf <- rep(NA,n)

  for(i in 1:n){
    qres <- qpois(res[i], lambda=lambda1f)-1
    pres <- ppois(qres,lambda = lambda1f)
    pres[i] <- 0
    empcdf[i] <-sum(pres)/(n-1)
  }
  .dpit_finalize(empcdf, plot=plot, scale=scale, line_args = line_args, ...)
}
