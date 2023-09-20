#' QQ-plots of DPIT residuals
#'
#' Makes a QQ-plot of the DPIT residuals calculated from [resid_disc()], [resid_semiconti()] or [resid_zeroinfl()].
#' The plot should be close to the diagonal if the model is correctly specified.
#' Note that this function does not return residuals. To get both residuals and QQ-plot,
#' use [resid_disc()], [resid_semiconti()] and [resid_zeroinfl()].
#'
#' @usage qqresid(model, scale="normal")
#'
#' @param model Fitted model object (e.g., `glm()`, `glm.nb()`, `zeroinfl()`, and `polr()`)
#' @param scale You can choose the scale of the residuals among `normal` and `uniform` scales.
#' The sample quantiles of the residuals are plotted against
#' the theoretical quantiles of a standard normal distribution under the normal scale,
#' and against the theoretical quantiles of a uniform (0,1) distribution under the uniform scale.
#'  The defalut scale is `normal`.
#'
#'
#' @importFrom stats qqplot
#' @importFrom stats ppoints
#' @importFrom stats qnorm
#' @importFrom graphics abline
#' @export
#' @seealso [resid_disc()], [resid_semiconti()], [resid_zeroinfl()]
#'
#' @examples
#' n <- 100
#' b <- c(2, 1, -2)
#' x1 <- rnorm(n); x2 <- rbinom(n,1,0.7)
#' y <-  rpois(n, exp(b[1]+b[2]*x1+b[3]*x2))
#'
#' m1 <- glm(y~x1+x2, family=poisson)
#' qqresid(m1, scale="normal")
#' qqresid(m1, scale="uniform")


qqresid <- function(model, scale="normal"){
    glm.test <- (paste(model$call)[1] %in% c("glm", "glm.nb"))
    zero.test <- (paste(model$call)[1] %in% c("zeroinfl"))
    tweedie.test <- ifelse(model$family[[1]]=='Tweedie', T,F)
    polr.test <- (paste(model$call)[1] %in% c("polr"))


    if(glm.test) model.family <- paste(family(model)[[1]])
    if(zero.test) model.family <- model$dist

    title <- paste(paste(model$call)[1], model.family, "Model QQ-plot" , sep=" ")

    if(glm.test && !tweedie.test) empcdf <- rƒesid_disc(model, plot=T, scale)
    if(zero.test) empcdf <- resid_zeroinfl(model, plot=T, scale)
    if(glm.test && tweedie.test) empcdf <- resid_semiconti(model, plot=T, scale)
}

