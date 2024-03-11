#' Quasi Emprical residuals functions
#'
#' Draw the QQ-plot for regression models with discrete outcomes using quasi-empirical residual distribution functions.
#' Specifically, the model assumption of GLMs with binary, ordinal, Poisson, negative binomial,
#' zero-inlated Poisson, and zero-inflated negative binomial outcomes can be applicable to `resid_quasi()`.
#'
#'
#' @usage resid_quasi(model)
#'
#' @param model Model object (e.g., `glm`, `glm.nb`, `polr`, `zeroinfl`)
#'
#' @import np
#'
#' @returns A QQ plot.
#' \itemize{
#'  \item x-axis: Theoretical quantiles
#'  \item y-axis: Sample quantiles
#' }
#'
#' @details
#' detail here
#'
#' @references Lu Yang (2021). Assessment of Regression Models with Discrete Outcomes Using Quasi-Empirical Residual Distribution Functions, Journal of Computational and Graphical Statistics, 30(4), 1019-1035.
#' @export
#'
#' @examples
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
#' resid.nb1 <- resid_quasi(model1)
#'
#' # Overdispersion
#' model2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' resid.nb2 <- resid_quasi(model2)
#'
#' ## Zero inflated Poisson example
#'
#'
resid_quasi <- function(model){
  glm.test <- (paste(model$call)[[1]] %in% c("glm", "glm.nb"))
  polr.test <- (paste(model$call)[[1]] %in% c("polr"))
  zero.test <- (paste(model$call)[1] %in% c("zeroinfl"))


  if (glm.test) {
    model.family <- paste(family(model)[[1]])
    if (paste(unlist(strsplit(model.family, split = ""))[1:17], collapse = "") == "Negative Binomial") model.family <- "nb"
  }
  if (polr.test) model.family <- "multi"
  if (zero.test) {
    model.family <- model$dist
  }

  if (model.family == "nb" && glm.test) {
    empcdf <- resid.nb_quasi(model)
  }
  if (model.family == "poisson" && glm.test) {
    empcdf <- resid.pois_quasi(model)
  }
  if (model.family == "binomial" && glm.test) {
    empcdf <- resid.bin_quasi(model)
  }
  if (model.family == "multi" && polr.test) {
    empcdf <- resid.ordi_quasi(model)
  }
  if (model.family == "poisson" && zero.test) {
    empcdf <- resid.zpois_quasi(model)
  }
  if (model.family == "negbin" && zero.test) {
    empcdf <- resid.znb_quasi(model)
  }
}
