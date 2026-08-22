# Methods for DPIT residual objects

Methods for printing, summarizing, plotting, and extracting values from
objects returned by
[`dpit()`](https://jhlee1408.github.io/assessor/reference/dpit.md),
[`dpit_2pm()`](https://jhlee1408.github.io/assessor/reference/dpit_2pm.md),
and the distribution-specific DPIT calculators. The print method
displays the fitted model call or calls, when available, and the sample
size. The summary method reports quantiles, mean, and standard deviation
for the selected scale. The residuals method extracts the DPIT
residuals, and the plot method constructs a QQ plot of the DPIT
residuals against its reference distribution.

## Usage

``` r
# S3 method for class 'dpit'
print(x, ...)

# S3 method for class 'dpit'
residuals(object, scale = c("normal", "uniform"), ...)

# S3 method for class 'dpit'
summary(object, scale = c("normal", "uniform"), ...)

# S3 method for class 'summary.dpit'
print(x, ...)

# S3 method for class 'dpit'
plot(x, scale = c("normal", "uniform"), line_args = list(), ...)
```

## Arguments

- x:

  A `dpit` object to be printed or plotted by `print.dpit()` or
  `plot.dpit()`, or a `summary.dpit` object to be printed by
  `print.summary.dpit()`.

- ...:

  Additional arguments passed to or from methods. For `plot.dpit()`,
  graphical arguments are passed to
  [`stats::qqplot()`](https://rdrr.io/r/stats/qqnorm.html) for
  customizing the QQ plot (e.g., `pch`, `col`, `cex`, `xlab`, `ylab`).

- object:

  A `dpit` object whose residuals are extracted by `residuals.dpit()` or
  summarized by `summary.dpit()`.

- scale:

  You can choose the scale of the residuals among `normal` and
  `uniform`. The sample quantiles of the residuals are plotted against
  the theoretical quantiles of a standard normal distribution under the
  normal scale, and against the theoretical quantiles of a uniform (0,1)
  distribution under the uniform scale. The default scale is `normal`.

- line_args:

  A named list of graphical parameters passed to
  [`graphics::abline()`](https://rdrr.io/r/graphics/abline.html) to
  modify the reference (red) 45° line in the QQ plot. If left empty, a
  default red dashed line is drawn.

## Value

`residuals.dpit()` returns a numeric vector of DPIT residuals.
`summary.dpit()` returns an object of class `summary.dpit`. The print
and plot methods return their input invisibly.
