### Semicontinuous test_that()

## Tweedie
library(statmod)
library(tweedie)
n = 500
set.seed(1234)
x11 <- rnorm(n)
x12 <- rnorm(n)

beta0 <- 5
beta1 <- 1
beta2 <- 1
lambda1 <- exp(beta0 + beta1 * x11 + beta2 * x12)


set.seed(1234)
y1 <- rtweedie(n, mu = lambda1, xi = 1.6, phi = 10)

#### True model
model1 <-
  glm(y1 ~ x11 + x12,
      family = tweedie(var.power = 1.6, link.power = 0))
lambda1f <- model1$fitted.values
phi1f <- summary(model1)$dis
p1f <- dtweedie(rep(0,n),mu=lambda1f, xi=1.6,phi=phi1f)
cdf1 <- ptweedie(y1,mu=lambda1f, xi=1.6,phi=phi1f )
func <- ecdf(p1f)
newp <- cdf1*func(cdf1)
newp<- qnorm(newp)
names(newp) <- 1:n

#### Missing covariate
model2 <-
  glm(y1 ~ x11, family = tweedie(var.power = 1.6, link.power =0))

lambda1f <- model2$fitted.values
phi1f <- summary(model2)$dis
p1f <- dtweedie(rep(0,n),mu=lambda1f, xi=1.6,phi=phi1f)
cdf1 <- ptweedie(y1,mu=lambda1f, xi=1.6,phi=phi1f )
func <- ecdf(p1f)
newp2 <- cdf1*func(cdf1)
newp2 <- qnorm(newp2)
names(newp2) <- 1:n


# test_that: Tweedie
test_that("Tweedie residuals are the same with the reference residuals", {
  expect_identical(resid_semiconti(model1, plot=F), newp , tolerance =1e-5)
  expect_identical(resid_semiconti(model2, plot=F), newp2, tolerance =1e-5)
})


# ## Tobit: VGAM
# n <- 500
#
# #### Parameter
# beta0 <- 2
# beta1 <- 2
# beta2 <- 2
#
# #### Covariates
# set.seed(1234)
# x11 <- runif(n,-1, 1)
# x12 <- runif(n,-1, 1)
# lambda1 <- beta0 + beta1 * x11 + beta2 * x12
# sd0 <- 0.2
# yun <- rnorm(n, mean = lambda1, sd = sd0)
# y <- ifelse(yun >= 0, yun, 0)
#
# ##### True model
# library(VGAM)
# fit1 <-
#   vglm(y ~ x11 + x12, tobit(Upper = Inf, Lower = 0, lmu = "identitylink"))
# p1f <- pnorm(0,mean=fitted(fit1),sd=exp(coef(fit1)[2]))
# cdf1 <- pnorm(y,mean=fitted(fit1),sd=exp(coef(fit1)[2]))
# newp <- cdf1*ecdf(p1f)(cdf1)
#
#
# fit1miss <- vglm(y ~ x11, tobit(Upper = Inf, Lower = 0, lmu = "identitylink"))
# p1f <- pnorm(0,mean=fitted(fit1miss),sd=exp(coef(fit1miss)[2]))
# cdf1 <- pnorm(y,mean=fitted(fit1miss),sd=exp(coef(fit1miss)[2]))
# newp2 <- cdf1*ecdf(p1f)(cdf1)
#
#
#
# # test_that: Tobit: VGAM
# test_that("Tobit(VGAM) residuals are the same with the reference residuals", {
#   expect_identical(resid_semiconti(fit1, plot=F, scale="uniform"), newp , tolerance =1e-5)
#   expect_identical(resid_semiconti(fit1miss, plot=F,scale="uniform"), newp2, tolerance =1e-5)
# })
#
