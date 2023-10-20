#' Residuals for regression models with two-parts outcomes
#'
#' Calculates DPIT proposed residuals for two parts model for semi-continuous outcomes.
#' `resid_2pm` can be used either with `model0` and `model1` or with `part0` and `part1`.
#'
#' @usage resid_2pm(model0, model1, y, part0, part1, plot=TRUE, scale = "normal")
#'
#'
#' @param model0 model object for 0 outcomes (e.g., logistic regression)
#' @param model1 model object for continuous outcomes (e.g., gamma regression)
#' @param y outcome variables
#' @param part0 alternate argument via cdf of model0
#' @param part1 alternate argument via cdf of model1
#' @param plot A logical value indicating whether or not to return QQ-plot
#' @param scale You can choose the scale of the residuals among `normal` and `uniform` scales. The default scale is `normal`.
#'
#'
#'
#'
#' @returns residuals. If plot=TRUE, also produces a QQ plot.
#' @seealso [resid_semiconti()]
#'
#' @importFrom stats ecdf
#' @importFrom MASS gamma.dispersion
#' @export
#'
#' @examples
#' n = 500
#' beta10= 1; beta11=-2; beta12=-1
#' beta13 <- -1; beta14 <- -1; beta15 <- -2
#' x11<-rnorm(n); x12<-rbinom(n,size=1,prob=0.4)
#'
#' p1<-1/(1+exp(-(beta10+x11*beta11+x12*beta12)))
#' lambda1<-exp(beta13+beta14*x11+beta15*x12)
#' y2 <- rgamma(n,scale=lambda1/2,shape=2) # gamma case
#' y <- rep(0,n)
#' u <-runif(n,0,1)
#' ind1 <- which(u>=p1)
#' y[ind1] <- y2[ind1]
#'
#' # Putting Model
#' mgamma <- glm(y[ind1]~x11[ind1]+x12[ind1],family=Gamma(link = "log"))
#' m10 <- glm(y==0~x12+x11,family=binomial(link = "logit"))
#' resid_2pm(model0 = m10, model1 = mgamma, y= y)
#'
#' # Alternative arguments: cdf
#' cdfgamma <- pgamma(y[ind1],scale = mgamma$fitted.values*gamma.dispersion(mgamma),
#'                   shape=1/gamma.dispersion(mgamma))
#' p1f <- m10$fitted.values
#' resid_2pm(y=y, part0= p1f, cdfgamma)


resid_2pm <- function(model0, model1, y, part0, part1, plot=TRUE, scale = "normal"){
  if(missing(y)) stop("argument y is missing")
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


