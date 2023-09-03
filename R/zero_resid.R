#' Residuals for zeroinflated regression model
#'
#' `resid_zeroinfl` is used for calculating residuals based on probability integral transform.
#' A zeroinflated model from `pscl` is recommended in this package.
#'
#' @usage resid_zeroinfl(model, plot=TRUE)
#' @param model glm model object (eg. `zeroinfl()` from `pscl`)
#' @param plot  A logical value indicating whether or not to return QQ-plot
#'
#' @returns proposed residuals
#' @importFrom stats family
#' @export
#'
#' @examples
#' ## zeroinflated poisson model
#' library(pscl)
#' # simulation
#' n <- 1e3
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.7)
#'
#' b1 <- 2;b2 <- 1;b0 <- -1
#' b00 <- -2;b10 <- 2
#' lambda1 <- exp(b0 + b1 * x1 + b2 * x2)
#' p0 <- 1 / (1 + exp(-(b00 + b10*x1)))
#' y0 <- rbinom(n, size = 1, prob = 1 - p0)
#' y1 <- rpois(n, lambda1)
#' y <- ifelse(y0 == 0, 0, y1)
#'
#' m1 <- zeroinfl(y ~ x1 + x2 | x1, dist = "poisson", link = "logit")
#' zpois.resid <- resid_zeroinfl(m1) ## with qqplot
#' zpois.resid <- resid_zeroinfl(m1, plot=FALSE) ## without qqplot

resid_zeroinfl <- function(model = stop("model must be specified"), plot=TRUE){

  # Model checking
  zero.test <- (paste(model$call)[1] %in% c("zeroinfl"))
  if(zero.test) model.family <- model$dist
  else stop("model should be zeroinfl model")

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
    qqplot(qnorm(ppoints(n)),qnorm(empcdf2), main="QQ-plot", xlab = "Theoretical Quantiles", ylab = "Empirical Quantiles",
           cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5,xlim=c(-3,3),ylim=c(-3,3))
    abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5,xlim=c(-3,3),ylim=c(-3,3))
  }
  return(empcdf)
}
