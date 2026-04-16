# Tidy a ddml_rep Object

DML-specific tidy method. Adds `ensemble_type` and `aggregation`
columns. Delegates to `tidy.ral_rep` for the base table computation.

## Usage

``` r
# S3 method for class 'ddml_rep'
tidy(
  x,
  ensemble_idx = 1,
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

  A `ddml_rep` object.

- ensemble_idx:

  Integer index of the ensemble type to report. Defaults to 1. Set to
  `NULL` for all ensemble types.

- aggregation:

  Character string. Aggregation method.

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
`statistic`, `p.value`, `ensemble_type`, and `aggregation`.

## See also

[`summary.ddml_rep`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_rep.md)
for the aggregation equations.

## Examples

``` r
# \donttest{
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace")]
reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
                      learners = list(what = ols),
                      sample_folds = 2,
                      resamples = 3, silent = TRUE)
tidy(reps)
#>          term      estimate   std.error     statistic      p.value
#> 1          D1 -1.538460e-01 0.014767084 -10.418174102 2.048480e-25
#> 2 (Intercept) -4.980567e-05 0.006910411  -0.007207337 9.942494e-01
#>   ensemble_type aggregation
#> 1          nnls      median
#> 2          nnls      median
tidy(reps, conf.int = TRUE)
#>          term      estimate   std.error     statistic      p.value
#> 1          D1 -1.538460e-01 0.014767084 -10.418174102 2.048480e-25
#> 2 (Intercept) -4.980567e-05 0.006910411  -0.007207337 9.942494e-01
#>   ensemble_type aggregation    conf.low   conf.high
#> 1          nnls      median -0.18278900 -0.12490310
#> 2          nnls      median -0.01359396  0.01349435
# }
```
