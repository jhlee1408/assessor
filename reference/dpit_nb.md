# Residuals for regression models with negative binomial outcomes

Computes DPIT residuals for regression models with negative binomial
outcomes using the observed counts (`y`) and their fitted distributional
parameters (`mu`, `size`).

## Usage

``` r
dpit_nb(y, mu, size)
```

## Arguments

- y:

  An observed outcome vector.

- mu:

  A vector of fitted mean values.

- size:

  A dispersion parameter of the negative binomial distribution.

## Value

A `dpit` object containing DPIT residuals.

## Details

For formulation details on discrete outcomes, see
[`dpit`](https://jhlee1408.github.io/assessor/reference/dpit.md).

## Examples

``` r
## Negative Binomial example
library(MASS)
n <- 500
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.7)
### Parameters
beta0 <- -2
beta1 <- 2
beta2 <- 1
size1 <- 2
lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
# generate outcomes
y <- rnbinom(n, mu = lambda1, size = size1)

# True model
model1 <- glm.nb(y ~ x1 + x2)
y1 <- model1$y
fitted1 <- fitted(model1)
size1 <- model1$theta
dpit.nb1 <- dpit_nb(y=y1, mu=fitted1, size=size1)
resid.nb1 <- residuals(dpit.nb1)
plot(dpit.nb1)


# Overdispersion
model2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
y2 <- model2$y
fitted2 <- fitted(model2)
dpit.nb2 <- dpit_pois(y=y2, mu=fitted2)
resid.nb2 <- residuals(dpit.nb2)
plot(dpit.nb2)
```
