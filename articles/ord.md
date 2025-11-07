# Ordered curve

### Assessing the mean structure

[`ord_curve()`](https://jhlee1408.github.io/assessor/reference/ord_curve.md)creates
a plot to assess the mean structure of regression models. The plot
compares the cumulative sum of the response variable and its
hypothesized value. Deviation from the diagonal suggests the possibility
that the mean structure of the model is incorrect.

Ordered curve method is not restricted to a discrete outcome regression
model. This function also supports a continuous outcome regression such
as `lm` object as well as `glm`, `glm.nb`, or `polr`.

In the example below, the underlying model is a logistic regression with
the probability of 1 as
${logit}^{- 1}\left( \beta_{0} + \beta_{1}X_{1} + \beta_{2}X_{2} + \beta_{3}X_{1}X_{2} \right)$,
where
$\left( \beta_{0},\beta_{1},\beta_{2},\beta_{3} \right) = ( - 5,2,1,3)$,
$X_{1} \sim N(1,1)$, and $X_{2}$ is a dummy variable with a probability
of one equal to 0.7. For the misspecified model, the binary covariate
and the interaction term are omitted. \#### Example

``` r
library(assessor)
## Binary example of ordered curve
n <- 500
set.seed(1234)
x1 <-rnorm(n,1,1); x2 <- rbinom(n,1,0.7)
beta0 <- -5; beta1 <- 2; beta2<- 1; beta3 <- 3
q1 <-1/(1+exp(beta0+beta1*x1+beta2*x2+beta3*x1*x2))
y1 <- rbinom(n,size=1,prob = 1-q1)

par(mfrow=c(1,2))
model0 <- glm(y1~x1*x2,family =binomial(link = "logit") )
ord_curve(model0,thr=model0$fitted.values) 
model1 <- glm(y1~x1,family =binomial(link = "logit") )
ord_curve(model1,thr=x2) 
```

![](ord_files/figure-html/ord%20curve%20poisson-1.png) The figures above
illustrate ordered curves plots corresponding to `model0` and `model1`.
In the left panel, the curve closely aligns with the diagonal line,
indicating that the mean structure of `model0` is correctly specified.
On the contrary, `model1` exhibits a deviation from the diagonal line
due to the omission of the variable $x_{2}$. This misspecification,
coupled with the choice of the threshold value as $x_{2}$, results in
the observed discrepancy in the ordered curve in the right panel.

The substantial disparity between the observed curve in `model1` and the
diagonal line strongly suggests the necessity of including the variable
$x_{2}$ in the model. This inclusion is crucial for accurately capturing
the underlying mean structure.
