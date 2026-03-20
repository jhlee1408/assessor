
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
