#' Goodness-of-fit test for discrete outcome regression models
#'
#' Goodness-of-fit test for discrete-outcome regression models.
#' Works with GLMs (Poisson, binomial/logistic, negative binomial),
#' ordinal outcome regression (`MASS::polr`), and
#' zero-inflated regressions (zero-inflated Poisson and negative binomial via `pscl::zeroinfl()`).
#'
#' @usage gof_disc(model, B=1e2, seed=NULL)
#' @param model A fitted model object (e.g., from `glm()`, `polr()`, `glm.nb()` or `zeroinfl()`).
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
#' \deqn{S_n = \int_0^1 \{ \hat{H}(u) - u \}^2 du} measures the deviation of \eqn{\hat{h}(u,y,\mathbf{x})} from the
#' identity function, with p-values obtained by bootstrap. This method
#' complements residual-based diagnostics by providing a formal check of model adequacy.
#' @references Yang L, Genest C, Neslehova J (2025). “A goodness-of-fit test for regression models with discrete outcomes.” Canadian Journal of Statistics
#'
#' @returns Test statistics and p-values
#' @importFrom stats family
#' @import tweedie
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
gof_disc <- function(model, B = 1e2, seed = NULL) UseMethod("gof_disc")

#' @rawNamespace S3method(gof_disc,default)
gof_disc.default <- function(model, B = 1e2, seed = NULL) {
  cls <- paste(class(model), collapse = ", ")
  stop(
    sprintf(
      "Unsupported model class for gof_disc(): %s. Supported classes include glm, zeroinfl, and polr.",
      cls
    ),
    call. = FALSE
  )
}

#' @rawNamespace S3method(gof_disc,glm)
gof_disc.glm <- function(model, B = 1e2, seed = NULL) {
  .gof_glm(model, B = B, seed = seed)
}

.gof_glm <- function(model, B = 1e2, seed = NULL) {
  fam_obj <- stats::family(model)
  fam_raw  <- fam_obj$family
  link_raw <- fam_obj$link
  fam_key <- fam_raw
  if (startsWith(fam_raw, "Negative Binomial")) {
    fam_key <- "negbin"
  }
  key <- list(
    model  = model,
    family = fam_raw,
    link   = link_raw,
    B      = B,
    seed   = seed
  )
  class(key) <- fam_key
  .gof_glm_key(key)
}

.gof_glm_key <- function(key) {
  UseMethod(".gof_glm_key")
}


#' @rawNamespace S3method(.gof_glm_key,default)
.gof_glm_key.default <- function(key) {
  stop(
    sprintf(
      "Unsupported GLM family/link for gof_disc(): %s / %s",
      key$family, key$link
    ),
    call. = FALSE
  )
}

#' @rawNamespace S3method(.gof_glm_key, poisson)
.gof_glm_key.poisson <- function(key) {
  if (!exists("gof_pois", mode = "function", inherits = TRUE)) {
    stop("gof_pois() not found.", call. = FALSE)
  }
  gof_pois(B = key$B, poismodel = key$model, seed = key$seed)
}

#' @rawNamespace S3method(.gof_glm_key,binomial)
.gof_glm_key.binomial <- function(key) {
  if (!exists("gof_bin", mode = "function", inherits = TRUE)) {
    stop("gof_bin() not found.", call. = FALSE)
  }
  gof_bin(B = key$B, bimodel = key$model, seed = key$seed)
}

#' @rawNamespace S3method(.gof_glm_key,negbin)
.gof_glm_key.negbin <- function(key) {
  if (!exists("gof_nb", mode = "function", inherits = TRUE)) {
    stop("gof_nb() not found.", call. = FALSE)
  }
  gof_nb(B = key$B, nbmodel = key$model, seed = key$seed)
}

#' @rawNamespace S3method(gof_disc,zeroinfl)
gof_disc.zeroinfl <- function(model, B = 1e2, seed = NULL) {
  .gof_zeroinfl(model, B = B, seed = seed)
}

.gof_zeroinfl <- function(model, B = 1e2, seed = NULL) {
  dist_raw <- model$dist

  key <- list(
    model = model,
    dist  = dist_raw,
    B     = B,
    seed  = seed
  )

  class(key) <- dist_raw
  .gof_zeroinfl_key(key)
}

.gof_zeroinfl_key <- function(key) {
  UseMethod(".gof_zeroinfl_key")
}

#' @rawNamespace S3method(.gof_zeroinfl_key,poisson)
.gof_zeroinfl_key.poisson <- function(key) {
  if (!exists("gof_zpois", mode = "function", inherits = TRUE)) {
    stop("gof_zpois() not found.", call. = FALSE)
  }
  gof_zpois(B = key$B, model1 = key$model, seed = key$seed)
}

#' @rawNamespace S3method(.gof_zeroinfl_key,negbin)
.gof_zeroinfl_key.negbin <- function(key) {
  if (!exists("gof_znb", mode = "function", inherits = TRUE)) {
    stop("gof_znb() not found.", call. = FALSE)
  }
  gof_znb(B = key$B, model1 = key$model, seed = key$seed)
}

#' @rawNamespace S3method(.gof_zeroinfl_key,default)
.gof_zeroinfl_key.default <- function(key) {
  stop(
    sprintf("Unsupported zeroinfl dist for gof_disc(): %s", key$dist),
    call. = FALSE
  )
}

#' @rawNamespace S3method(gof_disc,polr)
gof_disc.polr <- function(model, B = 1e2, seed = NULL) {
  if (!exists("gof_ordi", mode = "function", inherits = TRUE)) {
    stop("gof_ordi() not found.", call. = FALSE)
  }
  gof_ordi(B = B, multimodel = model, seed = seed)
}
