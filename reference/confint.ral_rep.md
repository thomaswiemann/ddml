# Confidence Intervals for RAL Rep Objects

Confidence Intervals for RAL Rep Objects

## Usage

``` r
# S3 method for class 'ral_rep'
confint(
  object,
  parm = NULL,
  level = 0.95,
  fit_idx = 1,
  aggregation = c("median", "mean", "spectral"),
  type = "HC1",
  uniform = FALSE,
  bootstraps = 999L,
  ...
)
```

## Arguments

- object:

  An object inheriting from class `ral_rep`.

- parm:

  Parameter specification (names or indices).

- level:

  Confidence level. Default 0.95.

- fit_idx:

  Integer index of the fit. Defaults to 1.

- aggregation:

  Character string. Aggregation rule.

- type:

  Character. HC type. Default `"HC1"`.

- uniform:

  Logical. Uniform bands via multiplier bootstrap? Default `FALSE`.

- bootstraps:

  Integer. Bootstrap draws. Default 999.

- ...:

  Currently unused.

## Value

A matrix with columns for lower and upper bounds. When `uniform = TRUE`,
the attribute `"crit_val"` contains the aggregated critical value.

## References

Chernozhukov V, Chetverikov D, Kato K (2013). "Gaussian approximations
and multiplier bootstrap for maxima of sums of high-dimensional random
vectors." Annals of Statistics, 41(6), 2786-2819.
