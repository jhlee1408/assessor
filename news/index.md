# Changelog

## assessor 1.3.1

CRAN release: 2026-04-20

- [`dpit()`](https://jhlee1408.github.io/assessor/reference/dpit.md) now
  rejects non-Tobit
  [`VGAM::vglm`](https://rdrr.io/pkg/VGAM/man/vglm.html) objects and
  non-constant-scale VGAM Tobit fits with informative errors.
- [`dpit_2pm()`](https://jhlee1408.github.io/assessor/reference/dpit_2pm.md)
  now requires exactly one model or supplied-probability input for each
  component and uses one common calculation path.
- Corrected the
  [`ord_curve()`](https://jhlee1408.github.io/assessor/reference/ord_curve.md)
  axis documentation and threshold label, and clarified the plotting
  documentation for
  [`ord_curve()`](https://jhlee1408.github.io/assessor/reference/ord_curve.md)
  and
  [`quasi_plot()`](https://jhlee1408.github.io/assessor/reference/quasi_plot.md).
- [`dpit_ordi()`](https://jhlee1408.github.io/assessor/reference/dpit_ordi.md)
  now respects the ordering supplied in `level` and validates the fitted
  category-probability matrix.
- Corrected the zero-inflated negative-binomial bandwidth calculation to
  use the negative-binomial CDF.
- Goodness-of-fit bootstrap refits now preserve offsets for binomial,
  Poisson, and negative-binomial regression models.
- [`gof_disc()`](https://jhlee1408.github.io/assessor/reference/gof_disc.md)
  now requires `B` to be a positive integer.
- Fixed an error where `scale = "uniform"` did not produce a uniformly
  scaled QQ plot.
- Fixed the negative-binomial goodness-of-fit bootstrap so that each
  replicate is refitted using the simulated response.
- [`dpit()`](https://jhlee1408.github.io/assessor/reference/dpit.md),
  [`dpit_2pm()`](https://jhlee1408.github.io/assessor/reference/dpit_2pm.md),
  and the distribution-specific calculators for binary, Poisson,
  negative binomial, ordinal, zero-inflated, Tobit, and Tweedie outcomes
  now return `dpit` objects. Use
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) to extract
  values, [`summary()`](https://rdrr.io/r/base/summary.html) to
  summarize them, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) to draw QQ
  plots. The `scale` argument is now supplied to
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), or
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
  graphical arguments are supplied to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).
- The [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods now label
  the number of residuals as the sample size.

## assessor 1.3.0

CRAN release: 2026-03-22

### Breaking changes

- Renamed all functions with prefix `resid_*` to `dpit_*`.
- Unified `resid_disc()`, `resid_zeroinfl()`, and `resid_semiconti()`
  into a single S3 generic function
  [`dpit()`](https://jhlee1408.github.io/assessor/reference/dpit.md).
- Renamed `resid_quasi()` to
  [`quasi_plot()`](https://jhlee1408.github.io/assessor/reference/quasi_plot.md).

### New features

- Added
  [`dpit_2pm()`](https://jhlee1408.github.io/assessor/reference/dpit_2pm.md)
  (formerly `resid_2pm()`).
- Introduced direct DPIT computation functions:
  - [`dpit_bin()`](https://jhlee1408.github.io/assessor/reference/dpit_bin.md)
  - [`dpit_pois()`](https://jhlee1408.github.io/assessor/reference/dpit_pois.md)
  - [`dpit_nb()`](https://jhlee1408.github.io/assessor/reference/dpit_nb.md)
  - [`dpit_ordi()`](https://jhlee1408.github.io/assessor/reference/dpit_ordi.md)
  - [`dpit_znb()`](https://jhlee1408.github.io/assessor/reference/dpit_znb.md)
  - [`dpit_zpois()`](https://jhlee1408.github.io/assessor/reference/dpit_zpois.md)
  - [`dpit_tobit()`](https://jhlee1408.github.io/assessor/reference/dpit_tobit.md)
  - [`dpit_tweedie()`](https://jhlee1408.github.io/assessor/reference/dpit_tweedie.md)
- Added
  [`gof_disc()`](https://jhlee1408.github.io/assessor/reference/gof_disc.md),
  an S3-based goodness-of-fit testing function for discrete outcome
  regression models.
- Updated
  [`quasi_plot()`](https://jhlee1408.github.io/assessor/reference/quasi_plot.md)
  to an S3-based interface.

## assessor 1.2.0

- Add `bballHR` data set which is used as an example for assessing a
  discrete outcome regression.

## assessor 1.1.0

CRAN release: 2024-04-03

- Add `resid_quasi()` function for QQ-plot by using kernel methods.
- Add `MEPS` data set which is used as an example for assessing a
  semicontinuous outcome regression.

## assessor 1.0.0

CRAN release: 2024-01-28

- Initial CRAN release
