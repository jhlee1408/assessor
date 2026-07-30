# assessor 1.3.1
* Fixed an error where `scale = "uniform"` did not produce a uniformly scaled QQ plot.
* Fixed the negative-binomial goodness-of-fit bootstrap so that each replicate is refitted using the simulated response.
* `dpit()`, `dpit_2pm()`, and the distribution-specific calculators for binary,
  Poisson, negative binomial, ordinal, zero-inflated, Tobit, and Tweedie
  outcomes now return `dpit` objects. Use `residuals()` to extract values,
  `summary()` to summarize them, and `plot()` to draw QQ plots. The `scale`
  argument is now supplied to `residuals()`, `summary()`, or `plot()`, and
  graphical arguments are supplied to `plot()`.
* The `print()` and `summary()` methods now label the number of residuals as
  the sample size.

# assessor 1.3.0
## Breaking changes
* Renamed all functions with prefix `resid_*` to `dpit_*`.
* Unified `resid_disc()`, `resid_zeroinfl()`, and `resid_semiconti()` into a single S3 generic function `dpit()`.
* Renamed `resid_quasi()` to `quasi_plot()`.

## New features
* Added `dpit_2pm()` (formerly `resid_2pm()`).
* Introduced direct DPIT computation functions:
  - `dpit_bin()`
  - `dpit_pois()`
  - `dpit_nb()`
  - `dpit_ordi()`
  - `dpit_znb()`
  - `dpit_zpois()`
  - `dpit_tobit()`
  - `dpit_tweedie()`
* Added `gof_disc()`, an S3-based goodness-of-fit testing function for discrete outcome regression models.
* Updated `quasi_plot()` to an S3-based interface.

# assessor 1.2.0
* Add `bballHR` data set which is used as an example for assessing a discrete outcome regression. 

# assessor 1.1.0

* Add `resid_quasi()` function for QQ-plot by using kernel methods. 
* Add `MEPS` data set which is used as an example for assessing a semicontinuous outcome regression. 

# assessor 1.0.0

* Initial CRAN release
