## two part test_that()\
#### Gamma + Logistic
library(MASS)
n <- 500
### Parameters
beta00 <-  -1
beta01 <- -2
beta02 <- -1
betas0 <- -1
betas1 <- -1
betas2 <- -2

### Covariates
set.seed(1234)
x11 <- rnorm(n)
x12 <- rbinom(n, size = 1, prob = 0.4)
p1 <- 1 / (1 + exp(-(beta00 + x11 * beta01 + x12 * beta02)))
lambda1 <- exp(betas0 + betas1 * x11 + betas2 * x12)

### Simulate positive data
y2 <- rgamma(n, scale = lambda1 / 2, shape = 2)

y <- rep(0, n)
u <- runif(n, 0, 1)
#### Simulate semicontinuous data
ind1 <- which(u >= p1)
y[ind1] <- y2[ind1]

#### Gamma regression for the positive part
mgamma <-
  glm(y[ind1] ~ x11[ind1] + x12[ind1], family = Gamma(link = "log"))
#### Logistic regression for the zero part
m10 <- glm(y == 0 ~ x11 + x12, family = binomial(link = "logit"))

cdfgamma <-
  pgamma(
    y[ind1],
    scale = mgamma$fitted.values * gamma.dispersion(mgamma),
    shape = 1 / gamma.dispersion(mgamma)
  )
p1f <- m10$fitted.values
cdf1 <- rep(0,n)
cdf1[y==0] <- m10$fitted.values[y==0]
cdf1[y>0] <- m10$fitted.values[which(y>0)] + (1-m10$fitted.values[which(y>0)])*cdfgamma
newp <- cdf1*ecdf(p1f)(cdf1)
newp <- qnorm(newp)

## test_that: 2pm = gamma + logitstic
test_that("2pm(gamma + logistic) residuals are the same with the reference residuals", {
  expect_identical(resid_2pm(part0 = m10$fitted.values,
                             part1 = cdfgamma,
                             y = y, plot=F), newp, tolerance =1e-5)
  expect_identical(resid_2pm(model0 = m10, model1 = mgamma,
                             y = y, plot=F), newp, tolerance =1e-5)
})
