#' DPIT residuals for regression models with various non-continuous outcomes
#'
#' Calculates DPIT residuals for regression models with non-continuous outcomes.
#' In particular, model assumptions for GLMs with discrete outcomes (e.g., binary, Poisson, and negative binomial), ordinal
#' regression models, zero-inflated regression models, and semicontinuous outcome
#' models can be assessed using \code{dpit()}.
#'
#' @usage dpit(model)
#' @param model A model object.
#'
#' @details
#' This function deploys the appropriate computation based on the class of
#' \code{model}. The supported model objects and outcome types are listed below.
#'
#' In addition to the class-based interface, the package also provides
#' distribution-specific DPIT residual calculators. If a fitted model comes from a
#' different class but has a supported outcome distribution, users can call
#' the corresponding distribution-based function directly.
#' For instance, for a regression model with Poisson outcomes,
#' one can use dpit to calculate the residuals if
#' the model is fit using `glm` function, or to use `dpit_pois` upon supplying fitted mean values.
#'
#' \itemize{
#' \item \strong{Discrete outcomes}
#'   \itemize{
#'   \item \code{\link[stats]{glm}} with \code{family = binomial()} (see \code{\link{dpit_bin}}).
#'   \item \code{\link[stats]{glm}} with \code{family = poisson()}
#'         (see \code{\link{dpit_pois}}).
#'   \item \code{\link[MASS]{glm.nb}} for negative binomial regression
#'         (see \code{\link{dpit_nb}}).
#'   \item \code{\link[MASS]{polr}} for ordinal outcomes
#'         (see \code{\link{dpit_ordi}}).
#'   }
#'
#' \item \strong{Zero-inflated discrete outcomes}
#'   \itemize{
#'   \item \code{\link[pscl]{zeroinfl}} with \code{dist = "poisson"}
#'         (see \code{\link{dpit_zpois}}).
#'   \item \code{\link[pscl]{zeroinfl}} with \code{dist = "negbin"}
#'         (see \code{\link{dpit_znb}}).
#'   }
#'
#' \item \strong{Semicontinuous outcomes}
#'   \itemize{
#'   \item Tobit regression via \code{\link[AER]{tobit}} from `AER` or \code{\link[VGAM]{vglm}} from `VGAM`
#'         (see \code{\link{dpit_tobit}}).
#'   \item Tweedie regression via \code{\link[stats]{glm}} with a Tweedie family
#'         (see \code{\link{dpit_tweedie}}).
#'   }
#' }
#'
#' \strong{Formulation for Discrete and Zero-Inflated Outcomes:}
#' \cr
#' The DPIT residual for the \eqn{i}th observation is defined as follows:
#' \deqn{\hat{r}(Y_i|X_i) = \hat{G}\bigg(\hat{F}_M(Y_i|\mathbf{X}_i)\bigg)}
#' where
#' \deqn{\hat{G}(s) = \frac{1}{n-1}\sum_{j=1, j \neq i}^{n}\hat{F}_M\bigg(\hat{F}_M^{(-1)}(\mathbf{X}_j)\bigg|\mathbf{X}_j\bigg)}
#' and \eqn{\hat{F}_M} refers to the fitted cumulative distribution function.
#' The `scale` argument is supplied to `residuals()`, `summary()`, or `plot()`, methods for further displaying basic object information, extracting residuals, summarizing residuals, and producing a QQ-plot, respectively.
#' When `scale="uniform"`, DPIT residuals should closely follow a uniform distribution, otherwise it implies model deficiency.
#' When `scale="normal"`, it applies the normal quantile transformation to the DPIT residuals
#' \deqn{\Phi^{-1}\left[\hat{r}(Y_i|\mathbf{X}_i)\right],i=1,\ldots,n.} The null pattern is the standard normal distribution in this case.
#' \cr
#'
#' \strong{Formulation for Semicontinuous Outcomes:}
#' \cr
#' The DPIT residuals for regression models with semicontinuous outcomes are \deqn{\hat{r}_i=\frac{\hat{F}_M(Y_i|\mathbf{X}_i)}{n}\sum_{j=1}^n1\left(\hat{p}_0(\mathbf{X}_j)\leq \hat{F}_M(Y_i|\mathbf{X}_i)\right), i=1,\ldots,n,}
#' where \eqn{\hat{p}_0(\mathbf{X}_i)} is the fitted probability of zero, and \eqn{\hat{F}_M(\cdot|\mathbf{X}_i)} is the  fitted cumulative distribution function for the \eqn{i}th observation. Furthermore, \deqn{\hat{F}_M(y|\mathbf{x})=\hat{p}_0(\mathbf{x})+\left(1-\hat{p}_0(\mathbf{x})\right)\hat{G}_M(y|\mathbf{x})}
#' where \eqn{\hat{G}_M} is the fitted cumulative distribution for the positive data.
#'
#' @returns A `dpit` object containing DPIT residuals.
#'
#'
#' @import stats
#' @import graphics
#' @export
#'
#'
#' @references
#' Yang, L. (2024). "Double probability integral transform residuals for regression models with discrete outcomes." \emph{Journal of Computational and Graphical Statistics}, 33(3), 787--803. \cr
#' Yang, L. (2024). "Diagnostics for regression models with semicontinuous outcomes." \emph{Biometrics}, 80(1), ujae007.
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
#' dpit.nb1 <- dpit(model1)
#' dpit.nb1
#' resid.nb1 <- residuals(dpit.nb1, scale = "uniform")
#' summary(dpit.nb1, scale = "uniform")
#' plot(dpit.nb1, scale = "uniform")
#'
#' # Overdispersion
#' model2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' dpit.nb2 <- dpit(model2)
#' resid.nb2 <- residuals(dpit.nb2, scale = "normal")
#' plot(dpit.nb2, scale = "normal")
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
#' dpit.bin1 <- dpit(model01)
#' resid.bin1 <- residuals(dpit.bin1)
#' plot(dpit.bin1)
#'
#' # Missing covariates
#' model02 <- glm(y1 ~ x1, family = binomial(link = "logit"))
#' dpit.bin2 <- dpit(model02)
#' resid.bin2 <- residuals(dpit.bin2)
#' plot(dpit.bin2)
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
#' dpit.poi1 <- dpit(poismodel1)
#' resid.poi1 <- residuals(dpit.poi1)
#' plot(dpit.poi1)
#'
#' # Enlarge three outcomes
#' y <- rpois(n, lambda1) + c(rep(0, (n - 3)), c(10, 15, 20))
#' poismodel2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
#' dpit.poi2 <- dpit(poismodel2)
#' resid.poi2 <- residuals(dpit.poi2)
#' plot(dpit.poi2)
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
#' dpit.ord1 <- dpit(multimodel)
#' resid.ord1 <- residuals(dpit.ord1)
#' plot(dpit.ord1)
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
#' dpit.ord2 <- dpit(multimodel)
#' resid.ord2 <- residuals(dpit.ord2)
#' plot(dpit.ord2)
dpit <- function(model) {
  UseMethod("dpit")
}
