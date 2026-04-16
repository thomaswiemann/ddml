# Tidy a RAL Object

Extracts coefficient estimates, standard errors, test statistics, and
p-values in a tidy data frame.

## Usage

``` r
# S3 method for class 'ral'
tidy(
  x,
  fit_idx = 1,
  conf.int = FALSE,
  conf.level = 0.95,
  type = "HC1",
  uniform = FALSE,
  bootstraps = 999L,
  ...
)
```

## Arguments

- x:

  An object inheriting from class `ral`.

- fit_idx:

  Integer index of the fit to report. Defaults to 1. Set to `NULL` for
  all fits.

- conf.int:

  Logical. Include confidence intervals? Default `FALSE`.

- conf.level:

  Confidence level. Default 0.95.

- type:

  Character. HC type. Default `"HC1"`.

- uniform:

  Logical. Uniform confidence bands? Default `FALSE`.

- bootstraps:

  Integer. Bootstrap draws. Default 999.

- ...:

  Currently unused.

## Value

A `data.frame` with columns `term`, `estimate`, `std.error`,
`statistic`, `p.value`, and `fit_label`.

## References

Chernozhukov V, Chetverikov D, Kato K (2013). "Gaussian approximations
and multiplier bootstrap for maxima of sums of high-dimensional random
vectors." Annals of Statistics, 41(6), 2786-2819.
