#' Residuals for regression models with zero-inflated outcomes
#'
#' Caluates the DPIT residuals for a regression model with zero-inflated discrete outcome.
#' A zero-inflated model from `pscl` is used in this function.
#'
#' @usage resid_zeroinfl(model, plot=TRUE, scale='normal')
#' @param model Model object, which is the output of `pscl::zeroinfl`.
#' @param plot  A logical value indicating whether or not to return QQ-plot.
#' @param scale You can choose the scale of the residuals among `normal` and `uniform` scales. The default scale is `normal`.
#'
#' @returns DPIT residuals. If `plot=TRUE`, also produces a QQ plot.
#' @importFrom stats family
#' @export
#'
#' @references Yang, Lu. "Double Probability Integral Transform Residuals for Regression Models with Discrete Outcomes." arXiv preprint arXiv:2308.15596 (2023).
#'
#' @examples
#' ## Zero-Inflated Poisson
#' library(pscl)
#' n <- 500
#' set.seed(1234)
#' # Covariates
#' x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.7)
#' # Coefficients
#' beta0 <- -2; beta1 <- 2; beta2 <- 1
#' beta00 <- -2; beta10 <- 2
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
#' resid_zeroinfl(modelzero1,plot=TRUE, scale="uniform")
#'
#' ## Zero inflation
#' modelzero2 <- glm(y~x1+x2, family=poisson(link="log"))
#' resid_disc(modelzero2, plot = TRUE,scale="normal")


resid_zeroinfl <- function(model, plot=TRUE, scale="normal"){

  # Model checking
  zero.test <- (paste(model$call)[1] %in% c("zeroinfl"))
  if(zero.test) model.family <- model$dist
  else stop("model has to be pscl::zeroinfl")

  #
  if(model.family == "poisson" && zero.test){
    empcdf <- resid.zpois(model)
  }
  if(model.family == "negbin" && zero.test){
    empcdf <- resid.znb(model)
  }
  if(plot==T){
    empcdf2 <- empcdf[empcdf!=1]
    n <- length(empcdf2)
    qqplot(qnorm(ppoints(n)),qnorm(empcdf), main="QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
           cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5)
    abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5)
  }
  if(scale=="normal") empcdf <- qnorm(empcdf)
  return(empcdf)
}
