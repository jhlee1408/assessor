#' Residuals for regression models with binary outcomes
#'
#'
#' Computes DPIT residuals for regression models with binary outcomes
#' using the observed responses (`y`) and their fitted distributional parameters(`prob`).
#'
#' @usage dpit_bin(y, prob, plot=TRUE, scale="normal", line_args=list(), ...)
#' @param y An observed outcome vector.
#' @param prob A vector of fitted probabilities of one.
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
#'
#' @returns DPIT residuals.
#'
#' @details
#' For formulation details on discrete outcomes, see \code{\link{dpit_pois}}.
#'
#' @examples
#' ## Binary example
#' n <- 500
#' set.seed(1234)
#' # Covariates
#' x1 <- rnorm(n, 1, 1)
#' x2 <- rbinom(n, 1, 0.7)
#' # Coefficients
#' beta0 <- -5
#' beta1 <- 2
#' beta2 <- 1
#' beta3 <- 3
#' q1 <- 1 / (1 + exp(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x1 * x2))
#' y1 <- rbinom(n, size = 1, prob = 1 - q1)
#'
#' # True model
#' model01 <- glm(y1 ~ x1 * x2, family = binomial(link = "logit"))
#' fitted1 <- fitted(model01)
#' y1 <- model01$y
#' resid.bin1 <- dpit_bin(y=y1, prob=fitted1)
#'
#' # Missing covariates
#' model02 <- glm(y1 ~ x1, family = binomial(link = "logit"))
#' y2 <- model02$y
#' fitted2 <- fitted(model02)
#' resid.bin2 <- dpit_bin(y=y2, prob=fitted2)
#' @export
dpit_bin <- function(y, prob,
                     plot=TRUE,
                     scale="normal",
                     line_args=list(), ...) {
  n <- length(y)
  prob <- 1-prob
  # fitted.values
  res <- ifelse(y == 0, prob, 1)
  ses <- ifelse(y == 1, prob, 0)
  # residuals
  empcdf <- rep(NA,n)
  for(i in 1:n){
    if(i %in% which(y==1)) next
    qres <- 1*(res[i] >= prob)*(res[i] < 1)*0 + 1 * (res[i]<prob) * (res[i]<1)*(-1)+1*(res[i]==1)
    pres <- ifelse(qres ==0, prob, ifelse(qres == 1, 1, 0))
    pres[i] <- 0
    empcdf[i] <- sum(pres)/ (n-1)
  }

  for(i in 1:n){
    if(i %in% which(y==0)) next
    qses <- (1*(ses[i]==0)*(-2)+1*(ses[i]>=prob)*(ses[i]<1)*0+1*(ses[i] < prob)*(ses[i]<1)*(-1)+1*(ses[i]==1))+1
    pses <- ifelse(qses ==0, prob, ifelse(qses == 1, 1, 0))
    pses[i] <- 0
    empcdf[i] <- sum(pses)/ (n-1)
  }
  .dpit_finalize(empcdf, plot=plot, scale=scale, line_args = line_args,...)
}
