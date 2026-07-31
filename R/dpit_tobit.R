#' Residuals for a tobit model
#'
#' Computes DPIT residuals for tobit regression models using the observed
#' responses (`y`) and their corresponding fitted distributional parameters (`mu`, `sd`).
#'
#' @usage dpit_tobit(y, mu, sd)
#' @param y An observed outcome vector.
#' @param mu A vector of fitted mean values of latent variables.
#' @param sd A standard deviation of latent variables.
#' @returns A `dpit` object containing DPIT residuals.
#'
#' @details
#' For formulation details on semicontinuous outcomes, see \code{\link{dpit}}.
#'
#' @examples
#' ## Tobit regression model
#' library(VGAM)
#' n <- 500
#' beta13 <- 1
#' beta14 <- -3
#' beta15 <- 3
#'
#' set.seed(1234)
#' x11 <- runif(n)
#' x12 <- runif(n)
#' lambda1 <- beta13 + beta14 * x11 + beta15 * x12
#' sd0 <- 0.3
#' yun <- rnorm(n, mean = lambda1, sd = sd0)
#' y <- ifelse(yun >= 0, yun, 0)
#'
#' # Using VGAM package
#' # True model
#' fit1 <- vglm(formula = y ~ x11 + x12,
#'              tobit(Upper = Inf, Lower = 0, lmu = "identitylink"))
#' # Missing covariate
#' fit1miss <- vglm(formula = y ~ x11,
#'                  tobit(Upper = Inf, Lower = 0, lmu = "identitylink"))
#'
#' dpit.tobit1 <- dpit_tobit(y = y, mu = VGAM::fitted(fit1), sd = sd0)
#' resid.tobit1 <- residuals(dpit.tobit1)
#' plot(dpit.tobit1)
#' dpit.tobit2 <- dpit_tobit(y = y, mu = VGAM::fitted(fit1miss), sd = sd0)
#' resid.tobit2 <- residuals(dpit.tobit2)
#' plot(dpit.tobit2)
#'
#' # Using AER package
#' library(AER)
#' # True model
#' fit2 <- tobit(y ~ x11 + x12, left = 0, right = Inf, dist = "gaussian")
#' # Missing covariate
#' fit2miss <- tobit(y ~ x11, left = 0, right = Inf, dist = "gaussian")
#'
#' dpit.aer1 <- dpit_tobit(y = y, mu = fitted(fit2), sd = sd0)
#' resid.aer1 <- residuals(dpit.aer1)
#' plot(dpit.aer1)
#' dpit.aer2 <- dpit_tobit(y = y, mu = fitted(fit2miss), sd = sd0)
#' resid.aer2 <- residuals(dpit.aer2)
#' plot(dpit.aer2)
#' @export
dpit_tobit <- function(y, mu, sd) {
  y  <- as.numeric(y)
  mu <- as.numeric(mu)
  sd <- as.numeric(sd)
  p1f  <- pnorm(0, mean = mu, sd = sd)
  cdf1 <- pnorm(y, mean = mu, sd = sd)
  Fhat <- stats::ecdf(p1f)
  newp <- as.vector(cdf1 * Fhat(cdf1))

  .new_dpit(newp, method = "Tobit")
}

#' @rawNamespace S3method(dpit,vglm)
dpit.vglm <- function(model) {
  if (!"tobit" %in% model@family@vfamily) {
    stop(
      "dpit() supports vglm objects fitted with VGAM::tobit() only.",
      call. = FALSE
    )
  }

  coef_matrix <- stats::coef(model, matrix = TRUE)
  if (ncol(coef_matrix) < 2L || !"(Intercept)" %in% rownames(coef_matrix)) {
    stop("Cannot extract the Tobit scale parameter from the vglm object.", call. = FALSE)
  }

  scale_coef <- coef_matrix[, 2L]
  non_intercept <- setdiff(names(scale_coef), "(Intercept)")
  if (
    length(non_intercept) > 0L &&
      any(abs(scale_coef[non_intercept]) > sqrt(.Machine$double.eps))
  ) {
    stop(
      "dpit() currently supports VGAM Tobit models with constant scale only.",
      call. = FALSE
    )
  }

  y <- model@y
  fitted <- VGAM::fitted(model)
  link_fun <- get(model@misc$link[2], envir = asNamespace("VGAM"))
  sd <- do.call(link_fun,
                args = list(
                  theta = unname(scale_coef["(Intercept)"]),
                  inverse = TRUE
                ))

  out <- dpit_tobit(y = y, mu = fitted, sd = sd)
  .attach_model_call(out, model)
}

#' @rawNamespace S3method(dpit,tobit)
dpit.tobit <- function(model) {
  y <- as.numeric(model$y)
  fitted <- VGAM::fitted(model)
  sd <- summary(model)$scale

  out <- dpit_tobit(y = y, mu = fitted, sd = sd)
  .attach_model_call(out, model)
}
