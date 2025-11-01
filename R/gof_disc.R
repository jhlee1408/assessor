#' Goodness-of-fit test for discrete outcome regression models
#'
#' Goodness-of-fit diagnostics for discrete-outcome regression models.
#' Works with GLMs (Poisson, binomial/logistic, negative binomial),
#' ordinal outcome regression (`MASS::polr`), and
#' zero-inflated regressions (Poisson and negative binomial via `pscl::zeroinfl()`).
#'
#' @usage gof_disc(model, B=1e2, seed=NULL)
#' @param model A fitted model object (e.g., from `glm()`,`polr()`, `glm.nb()` or `zeroinfl()`).
#' @param B Number of bootstrap samples. Default is 1e2.
#' @param seed random seed for bootstrap.
#'
#' @details
#' Let \eqn{(Y_i,\mathbf{X}_i),\ i=1,\ldots,n} denote independent observations, and let
#' \eqn{\hat F_M(\cdot \mid \mathbf{X}_i)} be the fitted model-based CDF.
#' It was shown in \emph{Yang (2025)} that under the correctly specified model,
#' \deqn{\hat{H}(u) = \frac{1}{n}\sum_{i=1}^n \hat{h}(u, Y_i, \mathbf{X}_i)}
#' should be close to the identity function, where
#' \deqn{\hat{h}(u, y, \mathbf{x}) =
#' \frac{u - \hat{F}_M (y-1 \mid \mathbf{x})}
#' {\hat{F}_M (y \mid \mathbf{x}) - \hat{F}_M (y-1 \mid \mathbf{x})}
#' \,\mathbf{1}\{ \hat{F}_M (y-1 \mid \mathbf{x}) < u < \hat{F}_M (y \mid \mathbf{x}) \}
#' + \mathbf{1}\{ u \ge \hat{F}_M (y \mid \mathbf{x}) \}.}
#' The test statistic
#' \deqn{S_n = \int_0^1 \{ \hat{H}(u) - u \}^2 du} measures deviation from the
#' identity function, with p-values obtained by bootstrap. This method
#' complements residual-based diagnostics by providing an inferential
#' check of model adequacy.
#' @references Yang L, Genest C, Neslehova J (2025). “A goodness-of-fit test for regression models with discrete outcomes.” Canadian Journal of Statistics
#'
#' @returns test statistics and p-values
#' @importFrom stats family
#' @export
#' @examples
#' library(MASS)
#' library(pscl)
#' n <- 500
#' B <- 1000
#' beta1 <- 1;  beta2 <- 1
#' beta0 <- -2; beta00 <- -2; beta10 <- 2
#' size1 <- 2
#' set.seed(1)
#' x1 <- rnorm(n)
#' x2 <- rbinom(n,1,0.7)
#' lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
#' p0 <- 1 / (1 + exp(-(beta00 + beta10 * x1)))
#' y0 <- rbinom(n, size = 1, prob = 1 - p0)
#' y1 <- rnegbin(n, mu=lambda1, theta=size1)
#' y <- ifelse(y0 == 0, 0, y1)
#' model1 <- zeroinfl(y ~ x1 + x2 | x1, dist = "negbin", link = "logit")
#' gof_disc(model1)
gof_disc <- function(model, B=1e2, seed=NULL){
  # Model checking
  glm.test  <- inherits(model, "glm")  || inherits(model, "negbin")
  zero.test <- inherits(model, "zeroinfl")
  polr.test <- inherits(model, "polr")

  # Throw error if unsupported
  if (!(glm.test || zero.test || polr.test)) {
    stop("model input should be glm, glm.nb, zeroinfl, or polr.")
  }
  if (glm.test) {
    model.family <- family(model)[[1]]
    if (startsWith(model.family, "Negative Binomial")) {
      model.family <- "nb"
    }

  } else if (zero.test) {
    model.family <- model$dist

  } else if (polr.test) {
    model.family <- "polr"
  }
  ##
  if (model.family == "nb" && glm.test) {
    gof <- gof_nb(B=B, nbmodel=model, seed=seed)
  }
  if (model.family == "poisson" && glm.test) {
    gof <- gof_pois(B=B, poismodel=model, seed=seed)
  }
  if (model.family == "binomial" && glm.test) {
    gof <- gof_bin(B=B, bimodel=model, seed=seed)
  }
  if (model.family == "poisson" && zero.test) {
    gof <- gof_zpois(B=B, model1=model, seed=seed)
  }
  if (model.family == "negbin" && zero.test) {
    gof <- gof_znb(B=B, model1 =model, seed=seed)
  }
  if( model.family =="polr"&& polr.test){
    gof <- gof_ordi(B=B, multimodel = model, seed=seed)
  }
 return(gof)
}
