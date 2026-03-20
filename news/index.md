# Changelog

## assessor 1.3.0

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
