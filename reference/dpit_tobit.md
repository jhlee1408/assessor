# Residuals for a tobit model

Computes DPIT residuals for tobit regression models using the observed
responses (`y`) and their corresponding fitted distributional parameters
(`mu`, `sd`).

## Usage

``` r
dpit_tobit(y, mu, sd, plot=TRUE, scale="normal", line_args=list(), ...)
```

## Arguments

- y:

  An observed outcome vector.

- mu:

  A vector of fitted mean values of latent variables.

- sd:

  A standard deviation of latent variables.

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

For formulation details on semicontinuous outcomes, see
[`dpit`](https://jhlee1408.github.io/assessor/reference/dpit.md).

## Examples

``` r
## Tobit regression model
library(VGAM)
#> Loading required package: stats4
#> Loading required package: splines
n <- 500
beta13 <- 1
beta14 <- -3
beta15 <- 3

set.seed(1234)
x11 <- runif(n)
x12 <- runif(n)
lambda1 <- beta13 + beta14 * x11 + beta15 * x12
sd0 <- 0.3
yun <- rnorm(n, mean = lambda1, sd = sd0)
y <- ifelse(yun >= 0, yun, 0)

# Using VGAM package
# True model
fit1 <- vglm(formula = y ~ x11 + x12,
             tobit(Upper = Inf, Lower = 0, lmu = "identitylink"))
# Missing covariate
fit1miss <- vglm(formula = y ~ x11,
                 tobit(Upper = Inf, Lower = 0, lmu = "identitylink"))

resid.tobit1 <- dpit_tobit(y = y, mu = VGAM::fitted(fit1), sd = sd0)
resid.tobit2 <- dpit_tobit(y = y, mu = VGAM::fitted(fit1miss), sd = sd0)

# Using AER package
library(AER)
#> Loading required package: car
#> Loading required package: carData
#> Loading required package: lmtest
#> Loading required package: zoo
#> 
#> Attaching package: ‘zoo’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.Date, as.Date.numeric
#> 
#> Attaching package: ‘lmtest’
#> The following object is masked from ‘package:VGAM’:
#> 
#>     lrtest
#> Loading required package: sandwich
#> Loading required package: survival
#> 
#> Attaching package: ‘AER’
#> The following object is masked from ‘package:VGAM’:
#> 
#>     tobit
# True model
fit2 <- tobit(y ~ x11 + x12, left = 0, right = Inf, dist = "gaussian")
# Missing covariate
fit2miss <- tobit(y ~ x11, left = 0, right = Inf, dist = "gaussian")

resid.aer1 <- dpit_tobit(y = y, mu = fitted(fit2), sd = sd0)
resid.aer2 <- dpit_tobit(y = y, mu = fitted(fit2miss), sd = sd0)
```
