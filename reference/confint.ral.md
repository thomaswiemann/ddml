# Confidence Intervals for RAL Estimators

Computes confidence intervals for one or more parameters.

## Usage

``` r
# S3 method for class 'ral'
confint(
  object,
  parm = NULL,
  level = 0.95,
  fit_idx = 1,
  type = "HC1",
  uniform = FALSE,
  bootstraps = 999L,
  ...
)
```

## Arguments

- object:

  An object inheriting from class `ral`.

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  Confidence level. Default 0.95.

- fit_idx:

  Integer index of the fit. Defaults to 1.

- type:

  Character. HC type. Default `"HC1"`.

- uniform:

  Logical. If `TRUE`, computes uniform confidence bands using the
  multiplier bootstrap. Default `FALSE`.

- bootstraps:

  Integer number of bootstrap draws. Only used when `uniform = TRUE`.
  Default 999.

- ...:

  Currently unused.

## Value

A matrix with columns for lower and upper bounds. When `uniform = TRUE`,
the attribute `"crit_val"` contains the critical value.

## References

Chernozhukov V, Chetverikov D, Kato K (2013). "Gaussian approximations
and multiplier bootstrap for maxima of sums of high-dimensional random
vectors." Annals of Statistics, 41(6), 2786-2819.

## See also

[`vcov.ral`](https://www.thomaswiemann.com/ddml/reference/vcov.ral.md)
