#### resid_disc() test functions
## Negative Binomial
library(MASS)
n <- 500
set.seed(1234)

### Covariates
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
model1 <- glm.nb(y ~ x1 + x2)# True model
model2 <- glm(y ~ x1 + x2, family = poisson(link = "log")) # overdispersion

# Reference values
lambda1f <- model1$fitted.values
size1f <- summary(model1)$theta

res <- pnbinom(y,mu=lambda1f,size=size1f)
res <- matrix(rep(res,n),n,n,byrow=TRUE)
lambda1f <- matrix(rep(lambda1f,n),n,n,byrow = FALSE)

qres <- qnbinom(res,mu=lambda1f,size=size1f)-1
pres <- pnbinom(qres,mu=lambda1f,size=size1f)

diag(pres) <- 0
empcdf1 <- qnorm(apply(pres,2,sum)/(n-1))

## overdispersion
lambda2f <- model2$fitted.values

res <- ppois(y,lambda2f)
res <- matrix(rep(res,n),n,n,byrow=TRUE)
lambda2f <- matrix(rep(lambda2f,n),n,n,byrow = FALSE)

qres <- qpois(res,lambda2f)-1
pres <- ppois(qres,lambda2f)

diag(pres) <- 0
empcdf2 <- qnorm(apply(pres,2,sum)/(n-1))


# test_that: NB
test_that("NB residuals are the same with the reference residuals", {
  expect_identical(resid_disc(model1, plot=F), empcdf1, tolerance =1e-5)
  expect_identical(resid_disc(model2, plot=F), empcdf2, tolerance =1e-5)
})

## Poisson
n <- 500
set.seed(1234)
### Covariates
x1 <- rnorm(n);x2 <- rbinom(n, 1, 0.7)
### Coefficients
beta0 <- -2
beta1 <- 2
beta2 <- 1
lambda1 <- exp(beta0 + beta1 * x1 + beta2 * x2)
y <- rpois(n, lambda1)
poismodel1 <- glm(y ~ x1 + x2, family = poisson(link = "log")) # True
y3 <- rpois(n, lambda1) + c(rep(0, (n - 3)), c(10, 15, 20))
poismodel2 <- glm(y3 ~ x1 + x2, family = poisson(link = "log")) # enlarged

## Reference: Poisson
lambda1f3 <- poismodel1$fitted.values
res <- ppois(y,lambda=lambda1f3)
qres <- qpois(matrix(rep(res,n),n,n,byrow=TRUE),lambda=matrix(rep(lambda1f3,n),n,n,byrow = FALSE))-1
pres <- ppois(qres,lambda=matrix(rep(lambda1f3,n),n,n,byrow = FALSE))
diag(pres) <- 0
empcdf1 <- qnorm(apply(pres,2,sum)/(n-1))

lambda1f3 <- poismodel2$fitted.values
res <- ppois(y3,lambda=lambda1f3)
qres <- qpois(matrix(rep(res,n),n,n,byrow=TRUE),lambda=matrix(rep(lambda1f3,n),n,n,byrow = FALSE))-1
pres <- ppois(qres,lambda=matrix(rep(lambda1f3,n),n,n,byrow = FALSE))
diag(pres) <- 0
empcdf2 <- (qnorm(apply(pres,2,sum)/(n-1)))


# test_that: Poisson
test_that("Poisson residuals are the same with the reference residuals", {
  expect_identical(resid_disc(poismodel1, plot=F), empcdf1, tolerance =1e-5)
  expect_identical(resid_disc(poismodel2, plot=F), empcdf2, tolerance =1e-5)
})


## Binary
n <- 500
set.seed(1234)
### Covariates
x1 <- rnorm(n, 1, 1);x2 <- rbinom(n, 1, 0.7)
### Coefficients
beta0 <- -5
beta1 <- 2
beta2 <- 1
beta3 <- 3

q1 <- 1 / (1 + exp(beta0 + beta1 * x1 + beta2 * x2 + beta3 * x1 * x2))
y1 <- rbinom(n, size = 1, prob = 1 - q1)

model01 <- glm(y1 ~ x1 * x2, family = binomial(link = "logit")) # True
model02 <- glm(y1 ~ x1, family = binomial(link = "logit")) # Missing cov

## Reference: binary
q10 <- 1 - model01$fitted.values
res <- ifelse(y1==0,q10,1)
res <- matrix(rep(res,n),n,n,byrow=TRUE)
qres <- 1*(res>=q10)*(res<1)*0 + 1*(res<q10)*(res<1)*(-1)+ 1*(res==1)
q10 <- matrix(rep(q10,n),n,n,byrow = FALSE)
pres <-  ifelse(qres==0,q10,ifelse(qres==1,1,0))
diag(pres) <- 0

ses <- ifelse(y1==1,q10,0)
ses <- matrix(rep(ses,n),n,n,byrow=TRUE)
qses <- (1*(ses==0)*(-2)+1*(ses>=q10)*(ses<1)*0 + 1*(ses<q10)*(ses<1)*(-1)+1*(ses==1))+1
pses <-  ifelse(qses==0,q10,ifelse(qses==1,1,0))

diag(pses) <- 0
rempcdf01 <- qnorm(apply(pses,2,sum)/(n-1)) # 1
empcdf01 <- qnorm(apply(pres,2,sum)/(n-1)) # 0
cdf01 <- rep(NA,n)
cdf01[y1==1] <- rempcdf01[y1==1]
cdf01[y1==0] <- empcdf01[y1==0]

# missing cov
q10 <- 1 - model02$fitted.values
res <- ifelse(y1==0,q10,1)
res <- matrix(rep(res,n),n,n,byrow=TRUE)
qres <- 1*(res>=q10)*(res<1)*0 + 1*(res<q10)*(res<1)*(-1)+ 1*(res==1)
q10 <- matrix(rep(q10,n),n,n,byrow = FALSE)
pres <-  ifelse(qres==0,q10,ifelse(qres==1,1,0))
diag(pres) <- 0

ses <- ifelse(y1==1,q10,0)
ses <- matrix(rep(ses,n),n,n,byrow=TRUE)
qses <- (1*(ses==0)*(-2)+1*(ses>=q10)*(ses<1)*0 + 1*(ses<q10)*(ses<1)*(-1)+1*(ses==1))+1
pses <-  ifelse(qses==0,q10,ifelse(qses==1,1,0))

diag(pses) <- 0
rempcdf02 <- qnorm(apply(pses,2,sum)/(n-1)) # 1
empcdf02 <- qnorm(apply(pres,2,sum)/(n-1)) # 0
cdf02 <- rep(NA,n)
cdf02[y1==1] <- rempcdf02[y1==1]
cdf02[y1==0] <- empcdf02[y1==0]

# test_that: binary
test_that("Logistics residuals are the same with the reference residuals", {
  expect_identical(resid_disc(model01, plot=F), cdf01, tolerance =1e-5)
  expect_identical(resid_disc(model02, plot=F), cdf02, tolerance =1e-5)
})

## Ordinal
library(MASS)
n <- 500

set.seed(1234)
### Covariates
x1 <- rnorm(n, mean = 2)

### Coefficient
beta1 <- 3
p0 <- plogis(1, location = beta1 * x1)
p1 <- plogis(4, location = beta1 * x1) - p0
p2 <- 1 - p0 - p1

genemult <- function(p) {
  rmultinom(1, size = 1, prob = c(p[1], p[2], p[3]))
}


test <- apply(cbind(p0, p1, p2), 1, genemult)
y1 <- rep(0, n)
y1[which(test[1, ] == 1)] <- 0
y1[which(test[2, ] == 1)] <- 1
y1[which(test[3, ] == 1)] <- 2

multimodel <- polr(as.factor(y1) ~ x1, method = "logistic")

## non-prop
n <- 500
set.seed(1234)
x1 <- rnorm(n, mean = 2)

beta1 <- 3
beta2 <- 1

p0 <- plogis(1, location = beta1 * x1)
p1 <- plogis(4, location = beta2 * x1) - p0
p2 <- 1 - p0 - p1

test <- apply(cbind(p0, p1, p2), 1, genemult)
y1 <- rep(0, n)
y1[which(test[1, ] == 1)] <- 0
y1[which(test[2, ] == 1)] <- 1
y1[which(test[3, ] == 1)] <- 2

multimodel2 <- polr(as.factor(y1) ~ x1, method = "logistic")




