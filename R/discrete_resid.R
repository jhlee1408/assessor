#' Residuals for regression models with discrete outcomes
#'
#' Calculates the DPIT residuals for regression models with discrete outcomes.
#' Specifically, the model assumption of GLMs with binary, ordinal, Poisson,
#' and negative binomial outcomes
#' can be assessed using `resid_disc()`.
#'
#' @usage resid_disc(model, plot=TRUE, scale="normal")
#' @param model Model object (e.g., `glm`, `glm.nb`, `polr`)
#' @param plot A logical value indicating whether or not to return QQ-plot
#' @param scale You can choose the scale of the residuals among `normal` and `uniform` scales.
#' The sample quantiles of the residuals are plotted against
#' the theoretical quantiles of a standard normal distribution under the normal scale,
#' and against the theoretical quantiles of a uniform (0,1) distribution under the uniform scale.
#'  The defalut scale is `normal`.
#'
#' @returns DPIT residuals. If `plot=TRUE`, also produces a QQ plot.
#'
#'
#' @import stats
#' @import graphics
#' @export
#'
#' @details
#' The DPIT residual for the \eqn{i}th observation is defined as follows:
#' \deqn{\hat{r}(Y_i|X_i) = \hat{G}\bigg(\hat{F}(Y_i|\mathbf{X}_i)\bigg)}
#' where
#' \deqn{\hat{G}(s) = \frac{1}{n-1}\sum_{j=1, j \neq i}^{n}\hat{F}\bigg(\hat{F}^{(-1)}(\mathbf{X}_j)\bigg|\mathbf{X}_j\bigg)}
#' and \eqn{\hat{F}} refers to the fitted cumulative distribution function.
#' When `scale="uniform"`, DPIT residuals should closely follow a uniform distribution, otherwise it implies model deficiency.
#' When `scale="normal"`, it applies the normal quantile transformation to the DPIT residuals
#' \deqn{\Phi^{-1}\left[\hat{r}(Y_i|\mathbf{X}_i)\right],i=1,\ldots,n.} The null pattern is the standard normal distribution in this case.
#' \cr
#'
#' Check reference for more details.
#'
#'
#'
#'
#' @references Yang, Lu. "Double Probability Integral Transform Residuals for Regression Models with Discrete Outcomes." arXiv preprint arXiv:2308.15596 (2023).
#'
#' @examples
#' library(MASS)
#' n <- 500
#' set.seed(1234)
#' ## Negative Binomial example
#' # Covariates
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.7)
#' ### Parameters
#' beta0 <- -2
#' beta1 <- 2
#' beta2 <- 1
#' size1 <- 2
#' lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
#' # generate outcomes
#' y <- rnbinom(n, mu = lambda1, size = size1)
#'
#' # True model
#' model1 <- glm.nb(y ~ x1 + x2)
#' resid.nb1 <- resid_disc(model1, plot = TRUE, scale = "uniform")
#'
#' # Overdispersion
#' model2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' resid.nb2 <- resid_disc(model2, plot = TRUE, scale = "normal")
#'
#' ## Binary example
#' n <- 500
#' set.seed(1234)
#' # Covariates
#' x1 <- rnorm(n, 1, 1)
#' x2 <- rbinom(n, 1, 0.7)
#' # Coefficients
#' beta0 <- -5
#' beta1 <- 2
#' beta2 <- 1
#' beta3 <- 3
#' q1 <- 1 / (1 + exp(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x1 * x2))
#' y1 <- rbinom(n, size = 1, prob = 1 - q1)
#'
#' # True model
#' model01 <- glm(y1 ~ x1 * x2, family = binomial(link = "logit"))
#' resid.bin1 <- resid_disc(model01, plot = TRUE)
#'
#' # Missing covariates
#' model02 <- glm(y1 ~ x1, family = binomial(link = "logit"))
#' resid.bin2 <- resid_disc(model02, plot = TRUE)
#'
#' ## Poisson example
#' n <- 500
#' set.seed(1234)
#' # Covariates
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.7)
#' # Coefficients
#' beta0 <- -2
#' beta1 <- 2
#' beta2 <- 1
#' lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
#' y <- rpois(n, lambda1)
#'
#' # True model
#' poismodel1 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' resid.poi1 <- resid_disc(poismodel1, plot = TRUE)
#'
#' # Enlarge three outcomes
#' y <- rpois(n, lambda1) + c(rep(0, (n - 3)), c(10, 15, 20))
#' poismodel2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' resid.poi2 <- resid_disc(poismodel2, plot = TRUE)
#'
#' ## Ordinal example
#' n <- 500
#' set.seed(1234)
#' # Covariates
#' x1 <- rnorm(n, mean = 2)
#' # Coefficient
#' beta1 <- 3
#'
#' # True model
#' p0 <- plogis(1, location = beta1 * x1)
#' p1 <- plogis(4, location = beta1 * x1) - p0
#' p2 <- 1 - p0 - p1
#' genemult <- function(p) {
#'   rmultinom(1, size = 1, prob = c(p[1], p[2], p[3]))
#' }
#' test <- apply(cbind(p0, p1, p2), 1, genemult)
#' y1 <- rep(0, n)
#' y1[which(test[1, ] == 1)] <- 0
#' y1[which(test[2, ] == 1)] <- 1
#' y1[which(test[3, ] == 1)] <- 2
#' multimodel <- polr(as.factor(y1) ~ x1, method = "logistic")
#' resid.ord1 <- resid_disc(multimodel, plot = TRUE)
#'
#' ## Non-Proportionality
#' n <- 500
#' set.seed(1234)
#' x1 <- rnorm(n, mean = 2)
#' beta1 <- 3
#' beta2 <- 1
#' p0 <- plogis(1, location = beta1 * x1)
#' p1 <- plogis(4, location = beta2 * x1) - p0
#' p2 <- 1 - p0 - p1
#' genemult <- function(p) {
#'   rmultinom(1, size = 1, prob = c(p[1], p[2], p[3]))
#' }
#' test <- apply(cbind(p0, p1, p2), 1, genemult)
#' y1 <- rep(0, n)
#' y1[which(test[1, ] == 1)] <- 0
#' y1[which(test[2, ] == 1)] <- 1
#' y1[which(test[3, ] == 1)] <- 2
#' multimodel <- polr(as.factor(y1) ~ x1, method = "logistic")
#' resid.ord2 <- resid_disc(multimodel, plot = TRUE)
resid_disc <- function(model, plot = TRUE, scale = "normal") {
  if (!(scale %in% c("normal", "uniform"))) stop("scale has to be either normal or uniform")
  # Model checking
  glm.test <- (paste(model$call)[[1]] %in% c("glm", "glm.nb"))
  polr.test <- (paste(model$call)[[1]] %in% c("polr"))


  if (!glm.test && !polr.test) {
    stop("model has to be glm, glm.nb or polr")
  }

  if (glm.test) model.family <- paste(family(model)[[1]])
  if (polr.test) model.family <- "multi"
  if (paste(unlist(strsplit(model.family, split = ""))[1:17], collapse = "") == "Negative Binomial") model.family <- "nb"


  if (model.family == "nb" && glm.test) {
    empcdf <- resid.nb(model)
  }
  if (model.family == "poisson" && glm.test) {
    empcdf <- resid.pois(model)
  }

  if (model.family == "binomial" && glm.test) {
    empcdf <- resid.bin(model)
  }
  if (model.family == "multi" && polr.test) {
    empcdf <- resid.ordi(model)
  }

  if (plot == T) {
    if (scale == "normal") {
      empcdf <- qnorm(empcdf)
      n <- length(empcdf)
      qqplot(qnorm(ppoints(n)), (empcdf),
        main = "QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
        cex.lab = 1, cex.axis = 1, cex.main = 1.5, lwd = 1.5
      )
      abline(0, 1, col = "red", lty = 5, cex.lab = 2, cex.axis = 2, cex.main = 2, lwd = 1.5)
    }
    if (scale == "uniform") {
      n <- length(empcdf)
      qqplot(ppoints(n), empcdf,
        main = "QQ plot", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
        cex.lab = 1, cex.axis = 1, cex.main = 1.5, lwd = 1.5
      )
      abline(0, 1, col = "red", lty = 5, cex.lab = 2, cex.axis = 2, cex.main = 2, lwd = 1.5)
    }
  } else {
    if (scale == "normal") empcdf <- qnorm(empcdf)
    if (scale == "uniform") empcdf <- empcdf
  }
  return(list("DPIT" = empcdf,
              "CRPS" = crps_dpit(empcdf)))
}
