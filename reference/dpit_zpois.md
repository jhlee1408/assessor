# Residuals for regression models with zero-inflated Poisson outcomes

Computes DPIT residuals for regression models with zero-inflated Poisson
outcomes using the observed counts(`y`) and their fitted distributional
parameters(`mu`, `pzero`).

## Usage

``` r
dpit_zpois(y, mu, pzero, plot=TRUE, scale="normal", line_args=list(), ...)
```

## Arguments

- y:

  An observed outcome vector.

- mu:

  A vector of fitted mean values for the count (non-zero) component.

- pzero:

  A vector of fitted probabilities for the zero-inflation component.

- plot:

  A logical value indicating whether or not to return QQ-plot

- scale:

  You can choose the scale of the residuals among `normal` and
  `uniform`. The sample quantiles of the residuals are plotted against
  the theoretical quantiles of a standard normal distribution under the
  normal scale, and against the theoretical quantiles of a uniform (0,1)
  distribution under the uniform scale. The default scale is `normal`.

- line_args:

  A named list of graphical parameters passed to
  [`graphics::abline()`](https://rdrr.io/r/graphics/abline.html) to
  modify the reference (red) 45° line in the QQ plot. If left empty, a
  default red dashed line is drawn.

- ...:

  Additional graphical arguments passed to
  [`stats::qqplot()`](https://rdrr.io/r/stats/qqnorm.html) for
  customizing the QQ plot (e.g., `pch`, `col`, `cex`, `xlab`, `ylab`).

## Value

DPIT residuals.

## Details

For formulation details on discrete outcomes, see
[`dpit`](https://jhlee1408.github.io/assessor/reference/dpit.md).

## Examples

``` r
## Zero-Inflated Poisson
library(pscl)
n <- 500
set.seed(1234)
# Covariates
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.7)
# Coefficients
beta0 <- -2
beta1 <- 2
beta2 <- 1
beta00 <- -2
beta10 <- 2

# Mean of Poisson part
lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
# Excess zero probability
p0 <- 1 / (1 + exp(-(beta00 + beta10 * x1)))
## simulate outcomes
y0 <- rbinom(n, size = 1, prob = 1 - p0)
y1 <- rpois(n, lambda1)
y <- ifelse(y0 == 0, 0, y1)
## True model
modelzero1 <- zeroinfl(y ~ x1 + x2 | x1, dist = "poisson", link = "logit")
y1 <- modelzero1$y
mu1    <- stats::predict(modelzero1, type = "count")
pzero1 <- stats::predict(modelzero1, type = "zero")
resid.zero1 <- dpit_zpois(y= y1, pzero=pzero1, mu=mu1)


## Zero inflation
modelzero2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
y2 <- modelzero2$y
mu2    <- fitted(modelzero2)
resid.zero2 <- dpit_pois(y= y2, mu=mu2)


```
