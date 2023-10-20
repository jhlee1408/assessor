#' Residuals for regression models with semicontinuous outcomes
#'
#' Calculates residuals for regression models with semi-continuous outcomes.
#' Specifically, a tweedie regression model from `tweedie` package or a tobit regression model
#' from `VGAM`, `AER` packages is used in this fucntion.
#'
#' @usage resid_semiconti(model, plot=TRUE, scale = "normal")
#'
#'
#' @param model model object(e.g., `tweedie`, `vglm`, and `tobit`)
#' @param plot A logical value indicating whether or not to return QQ-plot
#' @param scale You can choose the scale of the residuals among `normal` and `uniform` scales. The default scale is `normal`.
#'
#' @returns residuals. If plot=TRUE, also produces a QQ plot.
#'
#' @details
#' The proposed residual for the \eqn{i}th observation is defined as follows:
#' \deqn{\hat{r}_i = \frac{\hat{F}(Y_i|X_i)}{n}\sum_{j=1}^{n}I\bigg(\hat{p}_0(X_j) \leq \hat{F}(Y_i|X_i)\bigg)}
#' , which has a null distribution of uniformity. \eqn{\hat{F}} refers to the fitted cumulative distribution function.
#'
#'
#'
#' @importFrom stats ecdf
#' @import tweedie
#' @importFrom tweedie dtweedie
#' @importFrom tweedie ptweedie
#' @importFrom VGAM fitted
#' @export
#'

resid_semiconti <- function(model, plot=TRUE, scale = "normal"){
  is.vglm <- isS4(model)
  if(!is.vglm){
    if(!is.null(model$family)) model.family <- model$family[[1]]
    else model.family <- "AER"
  }


  if(!is.vglm && model.family == "Tweedie"){
    y1 <- model$y
    p.max <- get("p",envir=environment(model$family$variance))
    n <- length(y1)
    lambda1f <- model$fitted.values
    phi1f <- summary(model)$dis
    p1f <- dtweedie(rep(0,n),mu=lambda1f, xi=p.max,phi=phi1f )
    cdf1 <- ptweedie(y1,mu=lambda1f, xi=p.max, phi=phi1f )

    func <- ecdf(p1f)

    newp <- cdf1*func(cdf1)
  }

  if(!is.vglm && model.family == "AER"){
    p1f <- pnorm(0,mean=VGAM::fitted(model),sd=summary(model)$scale)
    cdf1 <- pnorm(y,mean=VGAM::fitted(model),sd=summary(model)$scale)
    newp <- cdf1*ecdf(p1f)(cdf1)
  }

  if(is.vglm){
    y <- model@y
    p1f <- pnorm(0,mean=fitted(model),sd=exp(coef(model)[2]))
    cdf1 <- pnorm(y,mean=fitted(model),sd=exp(coef(model)[2]))
    newp <- cdf1*ecdf(p1f)(cdf1)
  }


  if(plot==T){
    if(scale=="normal"){
      newp <- qnorm(newp)
      n <- length(newp)
      qqplot(qnorm(ppoints(n)),(newp),main="QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
             cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5)
      abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5)
    }
    if(scale=="uniform"){
      newp <- newp
      n <- length(newp)
      qqplot(ppoints(n),newp, main="QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
             cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5)
      abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5)
    }
  }
  else{
    if(scale=="normal") newp <- qnorm(newp)
    if(scale=="uniform") newp <- newp
  }
  return(newp)
}


