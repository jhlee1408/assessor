# Quasi empirical residuals functions

Draw the quasi-empirical residual distribution functions for regression
models with discrete outcomes. Specifically, the model assumption of
GLMs with binary, ordinal, Poisson, negative binomial, zero-inflated
Poisson, and zero-inflated negative binomial outcomes can be assessed
using `quasi_plot()`. A plot far apart from the diagonal indicates lack
of fit.

## Usage

``` r
quasi_plot(model, line_args=list(), ...)
```

## Arguments

- model:

  Model object (e.g., `glm`, `glm.nb`, `polr`, `zeroinfl`)

- line_args:

  A named list of graphical parameters passed to
  [`graphics::abline()`](https://rdrr.io/r/graphics/abline.html) to
  modify the reference 45° line in the quasi-empirical plot. If left
  empty, a default red dashed line is drawn.

- ...:

  Additional graphical arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  customizing the quasi-empirical plot (e.g., `lty`, `col`, `lwd`,
  `xlab`, `ylab`).

## Value

Invisibly returns `NULL`. The function is called for its side effect of
plotting the quasi-empirical residual distribution function
\\\hat{U}(s;\beta)\\ against \\s\\.

## Details

The quasi-empirical residual distribution function is defined as
follows: \$\$\hat{U}(s; \beta) = \sum\_{i=1}^{n}
W\_{n}(s;\mathbf{X}\_{i},\beta) 1\[F(Y\_{i}\| X\_{i}) \<
H(s;X\_{i})\]\$\$ where \$\$W_n(s; \mathbf{X}\_i, \beta) =
\frac{K\[(H(s; \mathbf{X}\_i)-s)/ \epsilon_n\]}{\sum\_{j=1}^{n} K\[(H(s;
\mathbf{X}\_j)-s)/ \epsilon_n\]},\$\$ \\\epsilon_n\\ is the bandwidth
selected suing cross validation; \\H(s, X_i) = \mathrm{argmin}\_{F(k
\mid X_i)} \|F(k \mid X_i) - s\|\\, \\F\\ is the CDF and \\K\\ is a
bounded, symmetric, and Lipschitz continuous kernel.

## References

Yang, L. (2021). "Assessment of regression models with discrete outcomes
using quasi-empirical residual distribution functions." *Journal of
Computational and Graphical Statistics*, 30(4), 1019–1035.

## Examples

``` r
## Negative Binomial example
library(MASS)
# Covariates
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
resid.nb1 <- quasi_plot(model1)


# Overdispersion
model2 <- glm(y ~ x1 + x2, family = poisson(link = "log"))
resid.nb2 <- quasi_plot(model2)


## Zero inflated Poisson example
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
resid.zero1 <- quasi_plot(modelzero1)
```
