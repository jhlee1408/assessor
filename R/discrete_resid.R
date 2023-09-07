#' Residuals for discrete outcome regression
#'
#' `resid_disc` is used to calculate the residuals for a discrete outcome GLM
#' such as logistics(binary), ordinal, poisson, or negative binomial regression.
#'
#' @usage resid_disc(model, plot=TRUE)
#' @param model glm model object (eg. `glm`, `glm.nb`, `polr`)
#' @param plot A logical value indicating whether or not to return QQ-plot
#'
#' @returns residuals
#' @import stats
#' @import graphics
#' @export
#'
#'
#' @examples
#' library(glmnet)
#' library(MASS)
#' n <- 1e3
#'
#' ## negative binomial example
#' b <- c(2, 1, -2)
#' x1 <- rnorm(n); x2 <- rbinom(n,1,0.7)
#'
#' y <-  rnbinom(n, size= 2, mu=exp(b[1]+b[2]*x1+b[3]*x2))
#' m.nb <- glm.nb(y~x1+x2)
#' nb.resid <- resid_disc(m.nb) ## with qqplot
#' nb.resid <- resid_disc(m.nb, plot=FALSE) ## without qqplot
#'
#' ## poisson example
#' y <-  rpois(n, exp(b[1]+b[2]*x1+b[3]*x2))
#' m.pois <- glm(y~x1+x2, family=poisson)
#' pois.resid <- resid_disc(m.pois) ## with qqplot
#' pois.resid <- resid_disc(m.pois, plot=FALSE) ## without qqplot
#'
#' ## binary example
#' b <- c(-5,2,1,3)
#' x1<-rnorm(n,1,1); x2 <- rbinom(n,1,0.7)
#' q1<-1/(1+exp(b[1]+b[2]*x1+b[3]*x2+b[4]*x1*x2))
#' y1 <- rbinom(n,size=1,prob = 1-q1)
#' m.binary <- glm(y1~x1*x2,family =binomial(link = "logit") )
#' bin.resid <- resid_disc(m.binary) ## with qqplot
#' bin.resid <- resid_disc(m.binary, plot=FALSE) ## without qqplot
#
#'
#' ## ordinal example
#' x1<-rnorm(n,mean=2)
#' b1 <- 3
#' p0 <- plogis(1,location=b*x1); p1 <- plogis(4,location=b*x1)-p0; p2 <- 1-p0-p1
#' genemult <- function(p){
#'  rmultinom(1,size=1,prob=c(p[1],p[2],p[3]))
#' }
#' test <- apply(cbind(p0,p1,p2),1,genemult)
#' y1 <- rep(0,n)
#' y1[which(test[1,]==1)] <-0; y1[which(test[2,]==1)] <-1; y1[which(test[3,]==1)] <-2
#' m.ordi <- polr(as.factor(y1)~x1,method="logistic")
#' resid.ordi <- resid_disc(m.ordi) # with qqplot
#' resid.ordi <- resid_disc(m.ordi, plot=FALSE) # without qqplot


resid_disc <- function(model, plot=TRUE){

  # Model checking
  glm.test <- (paste(model$call)[1] %in% c("glm", "glm.nb"))
  polr.test <- (paste(model$call)[1] %in% c("polr"))

  if(glm.test) model.family <- paste(family(model)[[1]])
  if(polr.test) model.family <- "multi"

  if(paste(unlist(strsplit(model.family, split= ""))[1:17], collapse = '') == "Negative Binomial") model.family <- "nb"

  if(!glm.test && !polr.test){
    stop("model should have discrete random variable.")
  }
  if((glm.test || polr.test) && !(model.family %in% c("poisson", "binomial", "multi","nb"))){
      stop("Check if your model family is poisson, bionomial, multi or negative binomial. Otherwise, use resid.zeroinfl() for zeroinflated model or resid.semiconti() for tweedie model.")
  }

  #
  if(model.family == "nb" && glm.test){
    empcdf <- resid.nb(model)
  }
  if(model.family == "poisson" && glm.test){
    empcdf <- resid.pois(model)
  }

  if(model.family == "binomial" && glm.test){
    empcdf <- resid.bin(model)
  }
  if(model.family == "multi" && polr.test){
    empcdf <- resid.logi(model)
  }
  if(plot==T){
    empcdf2 <- empcdf[empcdf!=1]
    n <- length(empcdf2)
    qqplot(qnorm(ppoints(n)),qnorm(empcdf), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
           cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5)
    abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5)
  }
  return(empcdf)
}
