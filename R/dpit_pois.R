#' Residuals for regression models with poisson outcomes
#'
#' Computes DPIT residuals for Poisson outcomes regression using the observed counts (`y`) and their
#' corresponding fitted mean values (`mu`).
#'
#' @usage dpit_pois(y, mu)
#' @param y An observed outcome vector.
#' @param mu A vector of fitted mean values.
#' @returns A `dpit` object containing DPIT residuals.
#'
#' @details
#' For formulation details on discrete outcomes, see \code{\link{dpit}}.
#'
#' @examples
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
#' poismodel <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' y1 <- poismodel$y
#' p1f <- fitted(poismodel)
#' dpit.poi <- dpit_pois(y=y1, mu=p1f)
#' resid.poi <- residuals(dpit.poi)
#' plot(dpit.poi)
#'
#' @export
dpit_pois <- function(y, mu) {
  n <- length(y)
  lambda1f <- mu
  res <- ppois(y, lambda = lambda1f)
  empcdf <- rep(NA,n)

  for(i in 1:n){
    qres <- qpois(res[i], lambda=lambda1f)-1
    pres <- ppois(qres,lambda = lambda1f)
    pres[i] <- 0
    empcdf[i] <-sum(pres)/(n-1)
  }
  .new_dpit(empcdf, method = "Poisson")
}
