# Residuals for regression models with binary outcomes

Computes DPIT residuals for regression models with binary outcomes using
the observed responses (`y`) and their fitted distributional
parameters(`prob`).

## Usage

``` r
dpit_bin(y, prob)
```

## Arguments

- y:

  An observed outcome vector.

- prob:

  A vector of fitted probabilities of one.

## Value

A `dpit` object containing DPIT residuals.

## Details

For formulation details on discrete outcomes, see
[`dpit_pois`](https://jhlee1408.github.io/assessor/reference/dpit_pois.md).

## Examples

``` r
## Binary example
n <- 500
set.seed(1234)
# Covariates
x1 <- rnorm(n, 1, 1)
x2 <- rbinom(n, 1, 0.7)
# Coefficients
beta0 <- -5
beta1 <- 2
beta2 <- 1
beta3 <- 3
q1 <- 1 / (1 + exp(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x1 * x2))
y1 <- rbinom(n, size = 1, prob = 1 - q1)

# True model
model01 <- glm(y1 ~ x1 * x2, family = binomial(link = "logit"))
fitted1 <- fitted(model01)
y1 <- model01$y
dpit.bin1 <- dpit_bin(y=y1, prob=fitted1)
resid.bin1 <- residuals(dpit.bin1)
plot(dpit.bin1)


# Missing covariates
model02 <- glm(y1 ~ x1, family = binomial(link = "logit"))
y2 <- model02$y
fitted2 <- fitted(model02)
dpit.bin2 <- dpit_bin(y=y2, prob=fitted2)
resid.bin2 <- residuals(dpit.bin2)
plot(dpit.bin2)
```
