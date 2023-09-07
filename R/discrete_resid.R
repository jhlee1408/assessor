#' Residuals for discrete outcome regression
#'
#' `resid_disc` is used to calculate the newly proposed residuals for a discrete outcome GLM.
#' Specifically, the model assumption of logistics(binary), ordinal, poisson, or
#' negative binomial regression can be assessed through `resid_disc()`.
#'
#' @usage resid_disc(model, plot=TRUE)
#' @param model glm model object (eg. `glm`, `glm.nb`, `polr`)
#' @param plot A logical value indicating whether or not to return QQ-plot
#'
#' @returns The double probability integral transform residuals(DPIT residuals).
#'
#'
#' @import stats
#' @import graphics
#' @export
#'
#' @details
#' The double probability integral transform residuals (DPIT residuals) is defined as follows.
#' \deqn{\hat{G}_{M_i}(s) = \frac{1}{n-1}\sum_{j=1, j \neq i}^{n}\hat{F}_{M}\bigg(\hat{F}_{M}^{(-1)}(s|X_j)\bigg|X_j\bigg)}
#' \deqn{\hat{r}(Y_i|X_i) = \hat{G}_{M_i}\bigg(\hat{F}_{M}(Y_i|X_i)\bigg)}
#' where \eqn{\hat{F}_{M}} refers to the cdf of the fitted model.
#' DPIT residuals should closely follow a uniform distribution, otherwise it implies model deficiency.\cr
#'
#' Check refernce for more details.
#'
#'
#'
#'
#' @references Yang, Lu. "Double Probability Integral Transform Residuals for Regression Models with Discrete Outcomes." arXiv preprint arXiv:2308.15596 (2023).
#'
#' @examples
#' library(glmnet)
#' library(MASS)
#' n=500
#' set.seed(1234)
#' ## Negative Binomial example
#' # Covariates
#' x1<-rnorm(n); x2 <- rbinom(n,1,0.7)
#' ### Parameters
#' beta0 <- -2; beta1 <- 2; beta2<- 1
#' size1<- 2
#' lambda1<-exp(beta0+beta1*x1+beta2*x2)
#' # generate outcomes
#' y <- rnbinom(n, mu=lambda1, size=size1)
#'
#' # True model
#' model1 <- glm.nb(y~x1+x2)
#' resid_disc(model1,plot = TRUE)
#'
#' # Overdispersion
#' model2 <- glm(y~x1+x2,family = poisson(link = "log"))
#' resid_disc(model2,plot = TRUE)
#'
#' ## Binary example
#' n<- 500
#' set.seed(1234)
#' # Covariates
#' x1<-rnorm(n,1,1); x2 <- rbinom(n,1,0.7)
#' # Coefficients
#' beta0 <- -5; beta1 <- 2; beta2<- 1; beta3 <- 3
#' q1<-1/(1+exp(beta0+beta1*x1+beta2*x2+beta3*x1*x2))
#' y1 <- rbinom(n,size=1,prob = 1-q1)
#'
#' # True model
#' model01 <- glm(y1~x1*x2,family =binomial(link = "logit") )
#' resid_disc(model01,plot = TRUE)
#'
#' # Missing covariates
#' model02 <- glm(y1~x1,family =binomial(link = "logit") )
#' resid_disc(model02,plot = TRUE)

#' ## Poisson example
#' n <- 500
#' set.seed(1234)
#' # Covariates
#' x1<-rnorm(n); x2 <-  rbinom(n,1,0.7)
#' # Coefficients
#' beta0 <- -2; beta1 <- 2; beta2<- 1
#' lambda1<-exp(beta0+beta1*x1+beta2*x2)
#'
#' y <- rpois(n, lambda1)
#'
#' # True model
#' poismodel1 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' resid_disc(poismodel1,plot = TRUE)
#'
#' # Enlarge three outcomes
#' y <- rpois(n, lambda1)+c(rep(0,(n-3)),c(10,15,20))
#' poismodel2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' resid_disc(poismodel2,plot = TRUE)
#'
#' ## Ordinal example
#' n<- 500
#' set.seed(1234)
#' # Covariates
#' x1 <-rnorm(n,mean=2)
#' # Coefficient
#' beta1 <- 3
#'
#' # True model
#' p0 <- plogis(1,location=beta1*x1)
#' p1 <- plogis(4,location=beta1*x1)-p0
#' p2 <- 1-p0-p1
#' genemult <- function(p){
#'   rmultinom(1,size=1,prob=c(p[1],p[2],p[3]))
#'}
#' test <- apply(cbind(p0,p1,p2),1,genemult)
#' y1 <- rep(0,n)
#' y1[which(test[1,]==1)] <-0
#' y1[which(test[2,]==1)] <-1
#' y1[which(test[3,]==1)] <-2
#' multimodel <- polr(as.factor(y1)~x1,method="logistic")
#' resid_disc(multimodel,plot = TRUE)
#'
#' ## Non-Proportionality
#' n<- 500
#' set.seed(1234)
#' x1 <-rnorm(n,mean=2)
#' beta1 <- 3; beta2 <- 1
#' p0 <- plogis(1,location=beta1*x1)
#' p1 <- plogis(4,location=beta2*x1)-p0
#' p2 <- 1-p0-p1
#' genemult <- function(p){
#'  rmultinom(1,size=1,prob=c(p[1],p[2],p[3]))
#' }
#' test <- apply(cbind(p0,p1,p2),1,genemult)
#' y1 <- rep(0,n)
#' y1[which(test[1,]==1)] <-0
#' y1[which(test[2,]==1)] <-1
#' y1[which(test[3,]==1)] <-2
#' multimodel <- polr(as.factor(y1)~x1,method="logistic")
#' resid_disc(multimodel,plot = TRUE)


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
    qqplot(qnorm(ppoints(n)),qnorm(empcdf),main="QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
           cex.lab=1, cex.axis=1, cex.main=1.5,lwd=1.5)
    abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=1.5)
  }
  return(empcdf)
}
