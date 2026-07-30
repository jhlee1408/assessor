# Residuals for regression models with poisson outcomes

Computes DPIT residuals for Poisson outcomes regression using the
observed counts (`y`) and their corresponding fitted mean values (`mu`).

## Usage

``` r
dpit_pois(y, mu)
```

## Arguments

- y:

  An observed outcome vector.

- mu:

  A vector of fitted mean values.

## Value

A `dpit` object containing DPIT residuals.

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
dpit.poi <- dpit_pois(y=y1, mu=p1f)
resid.poi <- residuals(dpit.poi)
plot(dpit.poi)

```
