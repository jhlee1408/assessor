#' Ordered Curve
#'
#' `ord.curve` is used to assess the mean structure of generalized linear models. The result
#' plot can be used to compare the cumulative sum of the response variable and its hypothesized value.
#' Deviation from the diagonal suggests possibility that the mean structure of the model is incorrect.
#'
#'
#' @details
#'  If the mean structure is well specified in the model,
#' \eqn{(L_1(t), L_2(t))} should be close to each other.
#' If the curve is distant from the diagonal, it suggests incorrectness in the mean structure.
#' Moreover, if the curve is above the diagonal, the summation of the response is larger than
#' the fitted mean, which implies that the mean is underestimated, and vice versa. \cr
#'
#' The role of `thr`(threshold variable) is to determine the rule for summing each mean and outcome observation.
#' The candidate for `thr` could be any function of predictors such as a single predictor(eg. `x1`),
#' a linear combination of predictor(eg.`x1+x2`), or fitted values(eg. `fitted(model)`). \cr
#'
#'
#' For more details, see the reference paper.
#'
#' @references Yang, Lu. "Double Probability Integral Transform Residuals for Regression Models with Discrete Outcomes." arXiv preprint arXiv:2308.15596 (2023).
#'
#' @usage ord_curve(model, thr)
#'
#' @param model glm model object (eg. `glm`, `glm.nb`, `polr`)
#' @param thr threshold variable (predictor, fitted values, or variable to be included as a covariate)
#'
#' @importFrom graphics abline
#'
#' @examples
#' ## Binomial example of Ordered curve
#' n <- 500
#' set.seed(1234)
#' x1 <-rnorm(n,1,1); x2 <- rbinom(n,1,0.7)
#' beta0 <- -5; beta1 <- 2; beta2<- 1; beta3 <- 3
#' q1 <-1/(1+exp(beta0+beta1*x1+beta2*x2+beta3*x1*x2))
#' y1 <- rbinom(n,size=1,prob = 1-q1)
#'
#' model0 <- glm(y1~x1*x2,family =binomial(link = "logit") )
#' ord_curve(model0,thr=model0$fitted.values) # set the threshold as fitted values
#' model1 <- glm(y1~x1,family =binomial(link = "logit") )
#' ord_curve(model1,thr=x2) # set the threshold as a covariate
#'
#' ## Poisson example of Ordered curve
#' n <- 500
#' set.seed(1234)
#' x1 <- rnorm(n); x2 <- rnorm(n)
#' beta0 <- 0; beta1 <- 2; beta2 <- 1
#' lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
#'
#' y <- rpois(n, lambda1)
#' poismodel1 <- glm(y ~ x1+x2, family = poisson(link = "log"))
#' ord_curve(poismodel1,thr=poismodel1$fitted.values)
#'
#' poismodel2 <- glm(y ~ x1, family = poisson(link = "log"))
#' ord_curve(poismodel2,thr=poismodel2$fitted.values)
#' ord_curve(poismodel2,thr=x2)
#'
#' @export

ord_curve <- function(model, thr){
  y1 <- model$y
  q10 <- model$fitted.values

  plot(cumsum((q10[sort(thr,index.return=TRUE)$ix]))/sum(q10), # X
       cumsum(y1[sort(thr,index.return=TRUE)$ix])/sum(y1),
       main = paste("Z:", deparse(substitute(thr))),
       xlab =expression(L[2](t)), ylab = expression(L[1](t)),
       cex.lab=1, cex.axis=1, cex.main=1.5,lwd=2,type="l",ylim=c(0,1),xlim=c(0,1))

  abline(0,1,col="red",lty=5,cex.lab=2, cex.axis=2, cex.main=2,lwd=2)
}
