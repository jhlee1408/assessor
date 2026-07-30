# Residuals for regression models with tweedie outcomes

Computes DPIT residuals for Tweedie-distributed outcomes using the
observed responses (`y`), their fitted mean values (`mu`), the variance
power parameter (\\\xi\\), and the dispersion parameter (\\\phi\\).

## Usage

``` r
dpit_tweedie(y, mu, xi, phi)
```

## Arguments

- y:

  Observed outcome vector.

- mu:

  Vector of fitted mean values of each outcomes.

- xi:

  Value of \\\xi\\ such that the variance is \\Var\[Y\] = \phi\mu^\xi\\

- phi:

  Dispersion parameter \\\phi\\.

## Value

A `dpit` object containing DPIT residuals.

## Details

For formulation details on semicontinuous outcomes, see
[`dpit`](https://jhlee1408.github.io/assessor/reference/dpit.md).

## Examples

``` r
## Tweedie model
library(tweedie)
library(statmod)
n <- 300
x11 <- rnorm(n)
x12 <- rnorm(n)
beta0 <- 5
beta1 <- 1
beta2 <- 1
lambda1 <- exp(beta0 + beta1 * x11 + beta2 * x12)
y1 <- rtweedie(n, mu = lambda1, xi = 1.6, phi = 10)
# Choose parameter p
# True model
model1 <-
  glm(y1 ~ x11 + x12,
    family = tweedie(var.power = 1.6, link.power = 0)
  )
y1 <- model1$y
p.max <- get("p", envir = environment(model1$family$variance))
lambda1f <- model1$fitted.values
phi1f <- summary(model1)$dis
dpit.tweedie <- dpit_tweedie(y= y1, mu=lambda1f, xi=p.max, phi=phi1f)
resid.tweedie <- residuals(dpit.tweedie)
plot(dpit.tweedie)
```
