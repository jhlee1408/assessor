#' Ordered Curve
#'
#' `ord.curve` is used to plot the ordered-curve plot.
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
#'
#' @usage ord_curve(model, thr)
#'
#' @param model glm model object (eg. `glm`, `glm.nb`, `polr`)
#' @param thr threshold variable (predictor, fitted values, or variable to be included as a covariate)
#'
#' @importFrom graphics abline
#'
#' @examples
#' library(glmnet)
#' n <- 1e2
#' b <- c(2, -1, -1)
#' x1 <- rnorm(n); x2 <- rnorm(n)
#' y <-  rpois(n, exp(b[1]+b[2]*x1+b[3]*x2))
#'
#' m1 <- glm(y~x1+x2, family=poisson)
#'
#' ord_curve(m1, x1) # Ordered curve with a threshold, x1
#' ord_curve(m1, x2) # Ordered curve with a threshold, x2
#' ord_curve(m1, fitted(m1)) # Ordered curve with a threshold, fitted values
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
