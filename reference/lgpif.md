# LGPIF Data

Data from the Wisconsin Local Government Property Insurance Fund
(LGPIF). The LGPIF was established to provide property insurance for
local government entities that include counties, cities, towns,
villages, school districts, fire departments, and other miscellaneous
entities, and is administered by the Wisconsin Office of the Insurance
Commissioner. Properties covered under this fund include government
buildings, vehicles, and equipment.

## Usage

``` r
LGPIF
```

## Format

A data frame with 5677 rows and 41 variables:

- `PolicyNum`:

  Policy number

- `Year`:

  Policy year

- `ClaimBC`:

  Total building and contents (BC) claims in the year

- `ClaimIM`:

  Total inland marine (IM) claims in the year (contractor’s equipment)

- `ClaimPN`:

  Total comprehensive claims from new motor vehicles in the year

- `ClaimPO`:

  Total comprehensive claims from old motor vehicles in the year

- `ClaimCN`:

  Total collision claims from new vehicles in the year

- `ClaimCO`:

  Total collision claims from old vehicles in the year

- `TypeCity`:

  Indicator for city entity

- `TypeCounty`:

  Indicator for county entity

- `TypeMisc`:

  Indicator for miscellaneous entity

- `TypeSchool`:

  Indicator for school entity

- `TypeTown`:

  Indicator for town entity

- `TypeVillage`:

  Indicator for village entity

- `IsRC`:

  Indicator for replacement cost (motor vehicles)

- `CoverageBC`:

  Log coverage amount for building and contents (in millions of dollars)

- `lnDeductBC`:

  Log deductible amount for building and contents

- `NoClaimCreditBC`:

  Indicator for no BC claims in prior year

- `yAvgBC`:

  Average BC claim amount

- `FreqBC`:

  Frequency of BC claims

- `CoverageIM`:

  Log coverage amount for inland marine (in millions of dollars)

- `lnDeductIM`:

  Log deductible amount for inland marine

- `NoClaimCreditIM`:

  Indicator for no IM claims in prior year

- `yAvgIM`:

  Average IM claim amount

- `FreqIM`:

  Frequency of IM claims

- `CoveragePN`:

  Log coverage amount for comprehensive new vehicles (in millions of
  dollars)

- `NoClaimCreditPN`:

  Indicator for no PN claims in prior year

- `yAvgPN`:

  Average PN claim amount

- `FreqPN`:

  Frequency of PN claims

- `CoveragePO`:

  Log coverage amount for comprehensive old vehicles (in millions of
  dollars)

- `NoClaimCreditPO`:

  Indicator for no PO claims in prior year

- `yAvgPO`:

  Average PO claim amount

- `FreqPO`:

  Frequency of PO claims

- `CoverageCN`:

  Log coverage amount for collision of new vehicles (in millions of
  dollars)

- `NoClaimCreditCN`:

  Indicator for no CN claims in prior year

- `yAvgCN`:

  Average CN claim amount

- `FreqCN`:

  Frequency of CN claims

- `CoverageCO`:

  Log coverage amount for collision of old vehicles (in millions of
  dollars)

- `NoClaimCreditCO`:

  Indicator for no CO claims in prior year

- `yAvgCO`:

  Average CO claim amount

- `FreqCO`:

  Frequency of CO claims

## Source

https://sites.google.com/a/wisc.edu/jed-frees/

## References

Frees, E. W., Lee, G., & Yang, L. (2016). Multivariate
frequency-severity regression models in insurance. Risks, 4(1), 4.
