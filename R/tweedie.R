#' Residuals for semicontinuous outcome regression
#'
#' `resid.semiconti` is used to calculate residuals for semi-continuous outcomes regression such as tweedie model.
#' A model object of semicontinuous regression from `tweedie` package is recommended.
#'
#' @usage resid_semiconti(model, plot=TRUE)
#'
#' @param model model object(using tweedie family)
#' @param plot A logical value indicating whether or not to return QQ-plot
#'
#' @return residuals
#'
#' @importFrom stats ecdf
#' @import tweedie
#' @importFrom tweedie dtweedie
#' @importFrom tweedie ptweedie
#' @export
#' @examples
#' library(tweedie)
#' library(statmod)
#' n <- 1e3
#' x11<-rnorm(n)
#' x12 <- rnorm(n)
#' b1 <- 1; b2<- 1;b0 <- 5
#' lambda1<-exp(b0+b1*x11+b2*x12)
#' z <- runif(n,0,1)
#' y1 <- qtweedie(z, mu=lambda1, xi=1.6,phi=10)
#' out=tweedie.profile(y1~x11+x12,p.vec=seq(1.1,1.9,length=9),
#'                     method="interpolation",do.ci=TRUE,do.smooth=TRUE,do.plot=TRUE)
#' m1 <- glm(y1~x11+x12,family=tweedie(var.power=out$p.max,link.power=0))
#' tweedie.reisd <- resid_semiconti(m1) # with qqplot
#' tweedie.reisd <- resid_semiconti(m1, plot=FALSE) #without qqpplot


resid_semiconti <- function(model, plot=TRUE){
  model.family <- model$family[[1]]
  if(model.family == "Tweedie"){
    y1 <- model$y
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
      qqplot(qnorm(ppoints(n)),qnorm(empcdf), main=title, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
             cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5)
      abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5)
    }
  }
  return(newp)
}


