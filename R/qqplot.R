#' qqplot with DPIT residuals
#'
#' `qqresid` is fucntion for plotting QQ-plot with residuals based on the probability integral transform.
#' Note that this function does not return residuals. To get both residuals and QQ-plot,
#' use [resid_disc()], [resid_semiconti()] and [resid_zeroinfl()].
#'
#' @usage qqresid(model, scale="normal")
#'
#' @param model glm model object (eg. `glm()`, `glm.nb()`, `zeroinfl()`, and `polr()`)
#' @param scale You can choose the scale of qqplot among `normal` and `uniform` scales. The defalut scale is `normal`.
#'
#' @details
#' The calculated residuals follows the uniform distribution under the adequate model assumption.
#' Also, taking the normal quantile transformation to DPIT residuals,
#' \deqn{\Phi^{-1}[\hat{r}(Y_i|X_i)], i=1,2,3,...,n,}
#' the standard normal distribution serves as the null pattern. \cr
#'
#' Users can choose `scale` argument between `normal` and `uniform`. The default value is `normal` as an usual QQ plot.
#'
#'
#' @importFrom stats qqplot
#' @importFrom stats ppoints
#' @importFrom stats qnorm
#' @importFrom graphics abline
#' @export
#' @seealso [resid_disc()], [resid_semiconti()], [resid_zeroinfl()]
#'
#' @examples
#' library(glmnet)
#' n <- 1e2
#' b <- c(2, 1, -2)
#' x1 <- rnorm(n); x2 <- rbinom(n,1,0.7)
#' y <-  rpois(n, exp(b[1]+b[2]*x1+b[3]*x2))
#'
#' m1 <- glm(y~x1+x2, family=poisson)
#' qqresid(m1) ## qqplot of poisson regression


qqresid <- function(model, scale="normal"){
  if(is.null(resid)){

    glm.test <- (paste(model$call)[1] %in% c("glm", "glm.nb"))
    zero.test <- (paste(model$call)[1] %in% c("zeroinfl"))
    tweedie.test <- ifelse(model$family[[1]]=='Tweedie', T,F)


    if(glm.test) model.family <- paste(family(model)[[1]])
    if(zero.test) model.family <- model$dist

    title <- paste(paste(model$call)[1], model.family, "Model QQ-plot" , sep=" ")

    if(glm.test && !tweedie.test) empcdf <- resid_disc(model)
    if(zero.test) empcdf <- resid_zeroinfl(model)
    if(glm.test && tweedie.test) empcdf <- resid_semiconti(model)

    n <- length(empcdf)

    if(scale == "normal"){
      qqplot(qnorm(ppoints(n)),qnorm(empcdf), main=title, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
             cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5)
      abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5)
    }
    if(scale == "uniform"){
      qqplot(ppoints(n),empcdf, main=title, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
             cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5)
      abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5)
    }
  }
  if(is.null(model)){
    empcdf <- resid
    n <- length(empcdf)
    qqplot(qnorm(ppoints(n)),qnorm(empcdf), main=title, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
           cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5)
    abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5)
  }
}

