#' Residuals for regression models with two-parts outcomes
#'
#' Calculates DPIT residuals for semi-continuous outcome regression such as 2PM.
#'
#' @usage resid_2pm(model0, model1, y, part0, part1, plot=TRUE, scale = "normal")
#'
#' @param model0 glm model for 0 outcomes
#' @param model1 model for continuous outcoms
#' @param y outcome variables
#' @param part0 cdf for y=0
#' @param part1 cdf for y>0
#' @param plot A logical value indicating whether or not to return QQ-plot
#' @param scale You can choose the scale of the residuals among `normal` and `uniform` scales. The default scale is `normal`.
#'
#' @returns The double probability integral transform residuals(DPIT residuals).
#'
#'
#' @importFrom stats ecdf
#' @importFrom MASS gamma.dispersion
#' @export
#'

resid_2pm <- function(model0, model1, y, part0, part1, plot=TRUE, scale = "normal"){

  if(!missing(model0) && !missing(model1) && !missing(y)){
    if(model1$family[[1]] != "Gamma") stop("The continuous part should follow Gamma() family.")
    else{
      n <- length(y)
      cdfgamma <- pgamma(y[y==0],scale = model1$fitted.values*gamma.dispersion(model1),
                         shape=1/gamma.dispersion(model1))
      p1f <- model0$fitted.values
      cdf1 <- rep(0,n)
      cdf1[y==0] <- model0$fitted.values[y==0]
      cdf1[y>0] <- model0$fitted.values[which(y>0)] + (1-model0$fitted.values[which(y>0)])*cdfgamma

      newp <- cdf1*ecdf(p1f)(cdf1)
    }
  }

  if(!missing(part0) && !missing(part1) && !missing(y)){
    n <- length(y)
    cdf1 <- rep(0,n)
    cdf1[y==0] <- part0[which(y==0)] # part1 should be shorter than part0 from y information
    cdf1[y>0] <- part0[which(y>0)] + (1-part0[which(y>0)])*part1
    newp <- cdf1*ecdf(part0)(cdf1)
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


