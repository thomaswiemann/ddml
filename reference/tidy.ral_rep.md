# Tidy a RAL Rep Object

Tidy a RAL Rep Object

## Usage

``` r
# S3 method for class 'ral_rep'
tidy(
  x,
  fit_idx = 1,
  aggregation = c("median", "mean", "spectral"),
  type = "HC1",
  conf.int = FALSE,
  conf.level = 0.95,
  uniform = FALSE,
  bootstraps = 999L,
  ...
)
```

## Arguments

- x:

  An object inheriting from class `ral_rep`.

- fit_idx:

  Integer index of the fit. Defaults to 1. Set to `NULL` for all fits.

- aggregation:

  Character string. Aggregation rule.

- type:

  Character. HC type. Default `"HC1"`.

- conf.int:

  Logical. Include CIs? Default `FALSE`.

- conf.level:

  Confidence level. Default 0.95.

- uniform:

  Logical. Uniform CIs? Default `FALSE`.

- bootstraps:

  Integer. Bootstrap draws. Default 999.

- ...:

  Currently unused.

## Value

A `data.frame` with columns `term`, `estimate`, `std.error`,
`statistic`, `p.value`, `fit_label`, and `aggregation`.
