#' qqplot with DPIT residuals
#'
#' Makes a QQ-plot of the DPIT residuals calculated from `resid_disc`, `resid_semiconti` or `resid_zeroinfl`.
#' The plot should be close to the diagonal if the model is correctly specified.
#' Note that this function does not return residuals. To get both residuals and QQ-plot,
#' use [resid_disc()], [resid_semiconti()] and [resid_zeroinfl()].
#'
#' @usage qqresid(model, scale="normal")
#'
#' @param model Fitted model object (e.g., `glm()`, `glm.nb()`, `zeroinfl()`, and `polr()`)
#' @param scale You can choose the scale of qqplot among `normal` and `uniform` scales.
#' The sample quantiles of the residuals are plotted against
#' the theoretical quantiles of a standard normal distribution under the normal scale,
#' and against the theoretical quantiles of a uniform (0,1) distribution under the uniform scale.
#'  The defalut scale is `normal`.
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
#' n <- 1e2
#' b <- c(2, 1, -2)
#' x1 <- rnorm(n); x2 <- rbinom(n,1,0.7)
#' y <-  rpois(n, exp(b[1]+b[2]*x1+b[3]*x2))
#'
#' m1 <- glm(y~x1+x2, family=poisson)
#' qqresid(m1, scale="normal") ## qqplot of poisson regression
#' qqresid(m1, scale="uniform")


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

