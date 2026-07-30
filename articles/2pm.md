# Two part

### Two part models

The input arguments for the
[`dpit_2pm()`](https://jhlee1408.github.io/assessor/reference/dpit_2pm.md)
function differ from those of other functions in `assessor` package.
Specifically, users can utilize this function with either models or
Probability Integral Transform (PIT) as input. The function returns a
`dpit` object; [`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) extract,
summarize, and plot its residuals on the selected scale.

For instance, in evaluating the distribution assumptions of a two-part
model that combines a logistic and a gamma regression, you should
provide the logistic regression model object as the argument for
`model0` and the gamma regression model object for `model1`. We
recommend utilizing the model input when assessing a gamma + logistic
two-part model. Alternatively, users can directly use the PIT as input
if their two-part model is not a gamma+logistic combination. In such
cases, users should first calculate the PIT and then input into `part0`
and `part1`, respectively.

This function accommodates two combinations: either `model0` in
conjunction with `model1` or `part0` in conjunction with `part1`. Note
that it is essential to specify the `y` (outcome) values in the function
arguments.

The underlying model is a two-part model. The probability of zero is
``` math
p_0(\mathbf{X})=\text{logit}^{-1}\left(\beta_0+X_{1}\beta_{1}+X_{2}\beta_{2} \right),
```
where $`X_1`$ is a standard normal variable, $`X_2`$ is binary with
probability of one as 0.4, and
$`(\beta_0,\beta_{1},\beta_{2})=(1,-2,-1)`$.  
A gamma distribution is employed to generate positive data. The mean
function of the positive part is described as
``` math
\lambda_S=\exp\left(\beta_{0S}+\beta_{1S}X_1+\beta_{2S}X_2\right).
```
We let $`(\beta_{0S},\beta_{1S},\beta_{2S})=(-1,-1,-2)`$. The dispersion
parameter is set to be 0.5.

- Models as input
- PIT as input

``` r

library(assessor)
library(MASS)
n <- 500
beta10 <- 1
beta11 <- -2
beta12 <- -1
beta13 <- -1
beta14 <- -1
beta15 <- -2
x11 <- rnorm(n)
x12 <- rbinom(n, size = 1, prob = 0.4)

p1 <- 1 / (1 + exp(-(beta10 + x11 * beta11 + x12 * beta12)))
lambda1 <- exp(beta13 + beta14 * x11 + beta15 * x12)
y2 <- rgamma(n, scale = lambda1 / 2, shape = 2)
y <- rep(0, n)
u <- runif(n, 0, 1)
ind1 <- which(u >= p1)
y[ind1] <- y2[ind1]
```

``` r

# models as input
mgamma <- glm(y[ind1] ~ x11[ind1] + x12[ind1], family = Gamma(link = "log")) # Gamma regression
m10 <- glm(y == 0 ~ x12 + x11, family = binomial(link = "logit")) # logistic regression

dpit.models <- dpit_2pm(model0 = m10, model1 = mgamma, y = y)
resid.models <- residuals(dpit.models, scale = "normal")
summary(dpit.models, scale = "normal")
#> Summary of DPIT residuals
#> 
#> Model calls:
#> model0:
#> glm(formula = y == 0 ~ x12 + x11, family = binomial(link = "logit"))
#> model1:
#> glm(formula = y[ind1] ~ x11[ind1] + x12[ind1], family = Gamma(link = "log"))
#> 
#> Residual scale: normal
#> Sample size: 500
#> 
#>          Min.       1st Qu.        Median          Mean       3rd Qu. 
#> -3.0345717620 -0.6707005462 -0.0110073257 -0.0004765426  0.6739872290 
#>          Max.     Std. Dev. 
#>  3.2263350600  1.0120292496
plot(dpit.models, scale = "normal")
```

![](2pm_files/figure-html/2pm%202-1.png)

``` r

library(assessor)
library(MASS)
n <- 500
beta10 <- 1
beta11 <- -2
beta12 <- -1
beta13 <- -1
beta14 <- -1
beta15 <- -2
x11 <- rnorm(n)
x12 <- rbinom(n, size = 1, prob = 0.4)

p1 <- 1 / (1 + exp(-(beta10 + x11 * beta11 + x12 * beta12)))
lambda1 <- exp(beta13 + beta14 * x11 + beta15 * x12)
y2 <- rgamma(n, scale = lambda1 / 2, shape = 2)
y <- rep(0, n)
u <- runif(n, 0, 1)
ind1 <- which(u >= p1)
y[ind1] <- y2[ind1]

# PIT as input
mgamma <- glm(y[ind1] ~ x11[ind1] + x12[ind1], family = Gamma(link = "log")) # gamma regression
m10 <- glm(y == 0 ~ x12 + x11, family = binomial(link = "logit")) # logistic regression

cdfgamma <- pgamma(y[ind1],
  scale = mgamma$fitted.values * gamma.dispersion(mgamma),
  shape = 1 / gamma.dispersion(mgamma)
)
p1f <- m10$fitted.values

dpit.pit <- dpit_2pm(part0 = p1f, part1 = cdfgamma, y = y)
resid.pit <- residuals(dpit.pit, scale = "uniform")
summary(dpit.pit, scale = "uniform")
#> Summary of DPIT residuals
#> 
#> Residual scale: uniform
#> Sample size: 500
#> 
#>         Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
#> 0.0005183729 0.2553949759 0.5076375532 0.5019392730 0.7571177366 0.9996803816 
#>    Std. Dev. 
#> 0.2889051341
plot(dpit.pit, scale = "uniform")
```

![](2pm_files/figure-html/pit%20as%20input-1.png)
