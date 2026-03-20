# Residuals for regression models with poisson outcomes

Computes DPIT residuals for Poisson outcomes regression using the
observed counts (`y`) and their corresponding fitted mean values (`mu`).

## Usage

``` r
dpit_pois(y, mu, plot=TRUE, scale="normal", line_args=list(), ...)
```

## Arguments

- y:

  An observed outcome vector.

- mu:

  A vector of fitted mean values.

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
## Poisson example
n <- 500
set.seed(1234)
# Covariates
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.7)
# Coefficients
beta0 <- -2
beta1 <- 2
beta2 <- 1
lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
y <- rpois(n, lambda1)

# True model
poismodel <- glm(y ~ x1 + x2, family = poisson(link = "log"))
y1 <- poismodel$y
p1f <- fitted(poismodel)
resid.poi <- dpit_pois(y=y1, mu=p1f)

```
