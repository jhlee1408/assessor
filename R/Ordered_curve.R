#' Ordered curve for assessing mean structures
#'
#' Creates a plot to assess the mean structure of regression models. The
#' plot compares the cumulative sum of the response variable and its hypothesized value.
#' Deviation from the diagonal suggests the possibility that the mean structure of the model is incorrect.
#'
#'
#' @details
#' The ordered curve plots \deqn{\hat{L}_1(t)=\frac{\sum_{i=1}^n\left[Y_i1(Z_i\leq t)\right]}{\sum_{i=1}^nY_i}} against
#'  \deqn{\hat{L}_2(t)=\frac{\sum_{i=1}^n\left[\hat{\lambda}_i1(Z_i\leq t)\right]}{\sum_{i=1}^n\hat{\lambda}_i},}
#' where \eqn{\hat{\lambda}_i} is the fitted mean, and \eqn{Z_i} is the threshold variable. \cr
#'
#'  If the mean structure is correctly specified in the model,
#' \eqn{\hat L_1(t)} and \eqn{\hat L_2(t)} should be close to each other.
#'
#' If the curve is distant from the diagonal, it suggests incorrectness in the mean structure.
#' Moreover, if the curve is above the diagonal, the summation of the response is larger than
#' the fitted mean, which implies that the mean is underestimated, and vice versa. \cr
#'
#' The role of `thr` (threshold variable \eqn{Z}) is to determine the rule  for accumulating \eqn{\hat{\lambda}_i} and \eqn{Y_i}, \eqn{i=1,\ldots,n}
#' for the ordered curve.
#' The candidate for `thr` could be any function of predictors such as a single predictor (e.g., `x1`),
#' a linear combination of predictor (e.g., `x1+x2`), or fitted values (e.g., `fitted(model)`).
#' It can also be a variable being considered to be included in the mean function.
#' If a variable  leads to a large discrepancy between the ordered curve and the diagonal,
#'  including this variable in the mean function should be considered.
#'
#' For more details, see the reference paper.
#'
#' @references Yang, Lu. "Double Probability Integral Transform Residuals for Regression Models with Discrete Outcomes." arXiv preprint arXiv:2308.15596 (2023).
#'
#' @usage ord_curve(model, thr)
#'
#' @param model Regression model object (e.g., `glm`, `glm.nb`, `polr`, `lm`)
#' @param thr Threshold variable (e.g., predictor, fitted values, or variable to be included as a covariate)
#'
#' @importFrom graphics abline
#'
#' @examples
#' ## Binary example of ordered curve
#' n <- 500
#' set.seed(1234)
#' x1 <- rnorm(n, 1, 1)
#' x2 <- rbinom(n, 1, 0.7)
#' beta0 <- -5
#' beta1 <- 2
#' beta2 <- 1
#' beta3 <- 3
#' q1 <- 1 / (1 + exp(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x1 * x2))
#' y1 <- rbinom(n, size = 1, prob = 1 - q1)
#'
#' ## True Model
#' model0 <- glm(y1 ~ x1 * x2, family = binomial(link = "logit"))
#' ord_curve(model0, thr = model0$fitted.values) # set the threshold as fitted values
#'
#' ## Missing a covariate
#' model1 <- glm(y1 ~ x1, family = binomial(link = "logit"))
#' ord_curve(model1, thr = x2) # set the threshold as a covariate
#'
#' ## Poisson example of ordered curve
#' n <- 500
#' set.seed(1234)
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' beta0 <- 0
#' beta1 <- 2
#' beta2 <- 1
#' lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
#'
#' y <- rpois(n, lambda1)
#'
#' ## True Model
#' poismodel1 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' ord_curve(poismodel1, thr = poismodel1$fitted.values)
#'
#' ## Missing a covariate
#' poismodel2 <- glm(y ~ x1, family = poisson(link = "log"))
#' ord_curve(poismodel2, thr = poismodel2$fitted.values)
#' ord_curve(poismodel2, thr = x2)
#'
#' @export

ord_curve <- function(model, thr) {

  glm.test <- (paste(model$call)[1] %in% c("glm", "glm.nb"))
  polr.test <- (paste(model$call)[1] %in% c("polr"))
  lm.test <- (paste(model$call)[1] %in% c("lm"))

  if(glm.test | polr.test) y1 <- model$y
  if(lm.test) y1 <- model$model[,1]


  if (length(thr) != length(y1)) stop("Length of thr and model fitted value has to match")

  q10 <- model$fitted.values

  plot(cumsum((q10[sort(thr, index.return = TRUE)$ix])) / sum(q10), # X
    cumsum(y1[sort(thr, index.return = TRUE)$ix]) / sum(y1),
    main = paste("Z:", deparse(substitute(thr))),
    xlab = expression(L[2](t)), ylab = expression(L[1](t)),
    cex.lab = 1, cex.axis = 1, cex.main = 1.5, lwd = 2, type = "l", ylim = c(0, 1), xlim = c(0, 1)
  )

  abline(0, 1, col = "red", lty = 5, cex.lab = 2, cex.axis = 2, cex.main = 2, lwd = 2)
}
