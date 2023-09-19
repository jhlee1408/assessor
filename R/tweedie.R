#' Residuals for regression models with semicontinuous outcomes
#'
#' `resid.semiconti` is used to calculate newly proposed residuals for semi-continuous outcomes regression such as tweedie model.
#' A model object of semicontinuous regression from `tweedie` package is recommended.
#'
#' @usage resid_semiconti(model, plot=TRUE, scale = "normal")
#'
#'
#' @param model model object(using tweedie family)
#' @param plot A logical value indicating whether or not to return QQ-plot
#' @param scale You can choose the scale of the residuals among `normal` and `uniform` scales. The default scale is `normal`.
#'
#' @returns The double probability integral transform residuals(DPIT residuals).
#'
#' @details
#' The proposed residuals are defined as
#' \deqn{\hat{r}_i = \frac{\hat{F}(Y_i|X_i)}{n}\sum_{j=1}^{n}I\bigg(\hat{p}_0(X_j) \leq \hat{F}(Y_i|X_i)\bigg)}
#' , which has a null distribution of uniformity.
#'
#'
#'
#' @importFrom stats ecdf
#' @import tweedie
#' @importFrom tweedie dtweedie
#' @importFrom tweedie ptweedie
#' @export
#'

resid_semiconti <- function(model, plot=TRUE, scale = "normal"){
  model.family <- model$family[[1]]
  if(model.family == "Tweedie"){
    y1 <- model$y.
    p.max <- get("p",envir=environment(model$family$variance))
    n <- length(y1)
    lambda1f <- model$fitted.values
    phi1f <- summary(model)$dis
    p1f <- dtweedie(rep(0,n),mu=lambda1f, xi=p.max,phi=phi1f )
    cdf1 <- ptweedie(y1,mu=lambda1f, xi=p.max,phi=phi1f )

    func <- ecdf(p1f)

    newp <- cdf1*func(cdf1)
    if(plot==T){
      empcdf2 <- newp[newp!=1]
      n <- length(empcdf2)
      qqplot(qnorm(ppoints(n)),qnorm(newp), main="QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
             cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5)
      abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5)
    }
  }
  if(scale=="normal") newp <- qnorm(newp)
  return(newp)
}


