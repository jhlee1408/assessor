#' Quasi emprical residuals functions
#'
#' Draw the quasi-empirical residual distribution functions for regression models with discrete outcomes.
#' Specifically, the model assumption of GLMs with binary, ordinal, Poisson, negative binomial,
#' zero-inlated Poisson, and zero-inflated negative binomial outcomes can be assessed using `quasi_plot()`.
#' A plot far apart from the diagonal indicates lack of fit.
#'
#'
#' @usage quasi_plot(model, line_args=list(), ...)
#'
#' @param model Model object (e.g., `glm`, `glm.nb`, `polr`, `zeroinfl`)
#' @param line_args A named list of graphical parameters passed to
#'   \code{graphics::abline()} to modify the reference (red) 45° line
#'   in the QQ plot. If left empty, a default red dashed line is drawn.
#' @param ... Additional graphical arguments passed to
#'   \code{stats::qqplot()} for customizing the QQ plot (e.g., \code{lty},
#'   \code{col}, \code{lwd}, \code{xlab}, \code{ylab}).
#'
#' @import np
#'
#' @returns A plot of quasi-empirial residual distribution function \eqn{\hat{U}(s;\beta)} against \eqn{s}.
#'
#' @details
#' The quasi-empirical residual distribution function is defined as follows:
#' \deqn{\hat{U}(s; \beta) = \sum_{i=1}^{n} W_{n}(s;\mathbf{X}_{i},\beta) 1[F(Y_{i}| X_{i}) < H(s;X_{i})]}
#' where
#' \deqn{W_n(s; \mathbf{X}_i, \beta) = \frac{K[(H(s; \mathbf{X}_i)-s)/ \epsilon_n]}{\sum_{j=1}^{n} K[(H(s; \mathbf{X}_j)-s)/ \epsilon_n]},}
#' \eqn{\epsilon_n} is the bandwidth; #' \eqn{H(s, X_i) = \mathrm{argmin}_{F(k \mid X_i)} F(k \mid X_i - s)} and \eqn{K} is a bounded, symmetric, and Lipschitz continuous kernel.
#'
#' @references Lu Yang (2021). Assessment of Regression Models with Discrete Outcomes Using Quasi-Empirical Residual Distribution Functions, Journal of Computational and Graphical Statistics, 30(4), 1019-1035.
#' @export
#'
#' @examples
#' ## Negative Binomial example
#' library(MASS)
#' # Covariates
#' n <- 500
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
#' resid.nb1 <- quasi_plot(model1)
#'
#' # Overdispersion
#' model2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' resid.nb2 <- quasi_plot(model2)
#'
#' ## Zero inflated Poisson example
#' library(pscl)
#' n <- 500
#' set.seed(1234)
#' # Covariates
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.7)
#' # Coefficients
#' beta0 <- -2
#' beta1 <- 2
#' beta2 <- 1
#' beta00 <- -2
#' beta10 <- 2
#'
#' # Mean of Poisson part
#' lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
#' # Excess zero probability
#' p0 <- 1 / (1 + exp(-(beta00 + beta10 * x1)))
#' ## simulate outcomes
#' y0 <- rbinom(n, size = 1, prob = 1 - p0)
#' y1 <- rpois(n, lambda1)
#' y <- ifelse(y0 == 0, 0, y1)
#' ## True model
#' modelzero1 <- zeroinfl(y ~ x1 + x2 | x1, dist = "poisson", link = "logit")
#' resid.zero1 <- quasi_plot(modelzero1)
quasi_plot <- function(model, line_args = list(), ...) {
  UseMethod("quasi_plot")
}

#' @rawNamespace S3method(quasi_plot,glm)
quasi_plot.glm <- function(model, line_args = list(), ...) {
  .quasi_plot_glm(model, line_args = line_args, ...)
}

.quasi_plot_glm <- function(model, line_args = list(), ...) {
  fam_obj <- stats::family(model)
  fam_raw <- fam_obj$family
  link_raw <- fam_obj$link

  fam_key <- fam_raw
  if (startsWith(fam_raw, "Negative Binomial")) {
    fam_key <- "nb"
  }

  key <- list(
    model     = model,
    family    = fam_raw,
    link      = link_raw,
    line_args = line_args
  )

  class(key) <- fam_key
  .quasi_plot_glm_key(key, ...)
}

.quasi_plot_glm_key <- function(key, ...) {
  UseMethod(".quasi_plot_glm_key")
}

#' @rawNamespace S3method(.quasi_plot_glm_key,nb)
.quasi_plot_glm_key.nb <- function(key, ...) {
  if (!exists("resid.nb_quasi", mode = "function", inherits = TRUE)) {
    stop("resid.nb_quasi() not found.", call. = FALSE)
  }
  resid.nb_quasi(key$model, line_args = key$line_args, ...)
}

#' @rawNamespace S3method(.quasi_plot_glm_key,poisson)
.quasi_plot_glm_key.poisson <- function(key, ...) {
  if (!exists("resid.pois_quasi", mode = "function", inherits = TRUE)) {
    stop("resid.pois_quasi() not found.", call. = FALSE)
  }
  resid.pois_quasi(key$model, line_args = key$line_args, ...)
}

#' @rawNamespace S3method(.quasi_plot_glm_key,binomial)
.quasi_plot_glm_key.binomial <- function(key, ...) {
  if (!exists("resid.bin_quasi", mode = "function", inherits = TRUE)) {
    stop("resid.bin_quasi() not found.", call. = FALSE)
  }
  resid.bin_quasi(key$model, line_args = key$line_args, ...)
}

#' @rawNamespace S3method(.quasi_plot_glm_key,default)
.quasi_plot_glm_key.default <- function(key, ...) {
  stop(
    sprintf("Unsupported GLM family/link for quasi_plot(): %s / %s", key$family, key$link),
    call. = FALSE
  )
}


#' @rawNamespace S3method(quasi_plot,zeroinfl)
quasi_plot.zeroinfl <- function(model, line_args = list(), ...) {
  .quasi_plot_zeroinfl(model, line_args = line_args, ...)
}

.quasi_plot_zeroinfl <- function(model, line_args = list(), ...) {
  dist_raw <- model$dist

  key <- list(
    model     = model,
    dist      = dist_raw,
    line_args = line_args
  )

  class(key) <- dist_raw
  .quasi_plot_zeroinfl_key(key, ...)
}

.quasi_plot_zeroinfl_key <- function(key, ...) {
  UseMethod(".quasi_plot_zeroinfl_key")
}

#' @rawNamespace S3method(.quasi_plot_zeroinfl_key,poisson)
.quasi_plot_zeroinfl_key.poisson <- function(key, ...) {
  if (!exists("resid.zpois_quasi", mode = "function", inherits = TRUE)) {
    stop("resid.zpois_quasi() not found.", call. = FALSE)
  }
  resid.zpois_quasi(key$model, line_args = key$line_args, ...)
}

#' @rawNamespace S3method(.quasi_plot_zeroinfl_key,negbin)
.quasi_plot_zeroinfl_key.negbin <- function(key, ...) {
  if (!exists("resid.znb_quasi", mode = "function", inherits = TRUE)) {
    stop("resid.znb_quasi() not found.", call. = FALSE)
  }
  resid.znb_quasi(key$model, line_args = key$line_args, ...)
}

#' @rawNamespace S3method(.quasi_plot_zeroinfl_key,default)
.quasi_plot_zeroinfl_key.default <- function(key, ...) {
  stop(sprintf("Unsupported zeroinfl dist for quasi_plot(): %s", key$dist), call. = FALSE)
}


#' @rawNamespace S3method(quasi_plot,polr)
quasi_plot.polr <- function(model, line_args = list(), ...) {
  if (!exists("resid.ordi_quasi", mode = "function", inherits = TRUE)) {
    stop("resid.ordi_quasi() not found.", call. = FALSE)
  }
  resid.ordi_quasi(model, line_args = line_args, ...)
}


#' @rawNamespace S3method(quasi_plot,default)
quasi_plot.default <- function(model, line_args = list(), ...) {
  cls <- paste(class(model), collapse = ", ")
  stop(
    sprintf("Unsupported model class for quasi_plot(): %s. Supported: glm, zeroinfl, polr.", cls),
    call. = FALSE
  )
}
