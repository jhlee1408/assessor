# Goodness-of-fit test for discrete outcome regression models

Goodness-of-fit test for discrete-outcome regression models. Works with
GLMs (Poisson, binomial, negative binomial), ordinal outcome regression
([`MASS::polr`](https://rdrr.io/pkg/MASS/man/polr.html)), and
zero-inflated regressions (zero-inflated Poisson and negative binomial
fit via
[`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)).

## Usage

``` r
gof_disc(model, B=1e2, seed=NULL)
```

## Arguments

- model:

  A fitted model object (e.g., from
  [`glm()`](https://rdrr.io/r/stats/glm.html),
  [`polr()`](https://rdrr.io/pkg/MASS/man/polr.html),
  [`glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) or
  [`zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html)).

- B:

  A positive integer giving the number of bootstrap samples. Default is
  1e2.

- seed:

  Random seed for bootstrap.

## Value

An object of class `"htest"` containing the test statistic (`S`), the
number of bootstrap samples (`B`), the p-value, the method description,
and the model call.

## Details

Let \\(Y_i,\mathbf{X}\_i),\\ i=1,\ldots,n\\ denote independent
observations, and let \\\hat F_M(\cdot \mid \mathbf{X}\_i)\\ be the
fitted model-based CDF. It was shown in *Yang et al. (2026)* that under
the correctly specified model, \$\$\hat{H}(u) = \frac{1}{n}\sum\_{i=1}^n
\hat{h}(u, Y_i, \mathbf{X}\_i)\$\$ should be close to the identity
function, where \$\$\hat{h}(u, y, \mathbf{x}) = \frac{u - \hat{F}\_M
(y-1 \mid \mathbf{x})} {\hat{F}\_M (y \mid \mathbf{x}) - \hat{F}\_M (y-1
\mid \mathbf{x})} \\\mathbf{1}\\ \hat{F}\_M (y-1 \mid \mathbf{x}) \< u
\< \hat{F}\_M (y \mid \mathbf{x}) \\ + \mathbf{1}\\ u \ge \hat{F}\_M (y
\mid \mathbf{x}) \\.\$\$ The test statistic \$\$S_n = \int_0^1 \\
\hat{H}(u) - u \\^2 du\$\$ measures the deviation of
\\\hat{h}(u,y,\mathbf{x})\\ from the identity function, with p-values
obtained by bootstrap. This method complements residual-based
diagnostics by providing a formal check of model adequacy.

## References

Yang, L., Genest, C., & Nešlehová, J. G. (2026). "A goodness-of-fit test
for regression models with discrete outcomes." *Canadian Journal of
Statistics*, 54(2), e70046.

## Examples

``` r
library(MASS)
library(pscl)
n <- 100
beta1 <- 1;  beta2 <- 1
beta0 <- -2; beta00 <- -2; beta10 <- 2
size1 <- 2
set.seed(1)
x1 <- rnorm(n)
x2 <- rbinom(n,1,0.7)
lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
p0 <- 1 / (1 + exp(-(beta00 + beta10 * x1)))
y0 <- rbinom(n, size = 1, prob = 1 - p0)
y1 <- rnegbin(n, mu=lambda1, theta=size1)
y <- ifelse(y0 == 0, 0, y1)
model1 <- zeroinfl(y ~ x1 + x2 | x1, dist = "negbin", link = "logit")
gof_disc(model1, B=50)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> 
#>  Goodness-of-fit test for regression models with discrete outcomes
#>  (zero-inflated negative binomial)
#> 
#> data:  zeroinfl(formula = y ~ x1 + x2 | x1, dist = "negbin", link = "logit")
#> S = 0.033878, B = 50, p-value = 0.88
#> 
```
