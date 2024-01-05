## resid_zeroinfl() test_that code
## Zero-Inflated Poisson
library(pscl)
n <- 500
set.seed(1234)
### Covariates
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.7)

### Coefficients
beta0 <- -2
beta1 <- 2
beta2 <- 1


beta00 <- -2
beta10 <- 2

##### Mean of Poisson part
lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
#### Excess zero probability
p0 <- 1 / (1 + exp(-(beta00 + beta10 * x1)))

#### simulate outcomes
y0 <- rbinom(n, size = 1, prob = 1 - p0)
y1 <- rpois(n, lambda1)

y <- ifelse(y0 == 0, 0, y1)

modelzero1 <- zeroinfl(y ~ x1 + x2 | x1, dist = "poisson", link = "logit") # true
modelzero2 <- glm(y ~ x1 + x2, family = poisson(link = "log")) # simple poisson

meanpoisson <- exp(modelzero1$coefficients$count[1] + modelzero1$coefficients$count[2] * x1 + modelzero1$coefficients$count[3] * x2)
pzero <- 1 / (1 + exp(-(modelzero1$coefficients$zero[1] + modelzero1$coefficients$zero[2] * x1)))

res <-(pzero + (1 - pzero) * (ppois(y, meanpoisson)))
res <- matrix(rep(res,n),n,n,byrow=TRUE)
meanpoisson <- matrix(rep(meanpoisson,n),n,n,byrow = FALSE)
pzero <- matrix(rep(pzero,n),n,n,byrow = FALSE)
qres <- ifelse(res < (pzero + (1 - pzero) * (ppois(0, lambda = meanpoisson))), 0,
               qpois(pmax((res - pzero) / (1 - pzero), 0), lambda = meanpoisson))-1
pres <- ifelse(qres==-1,0,(pzero + (1 - pzero) * (ppois(qres, meanpoisson))))
diag(pres) <- 0
zero.resid1 <- qnorm(apply(pres,2,sum)/(n-1))
names(zero.resid1) <- 1:500



## test_that for zeroinfl
test_that("Zeroinflation-Poisson residuals are the same as the reference residuals", {
  expect_identical(resid_zeroinfl(modelzero1, plot=F),zero.resid1, tolerance= 1e-5)
})
