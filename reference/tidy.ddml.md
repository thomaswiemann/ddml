# Tidy a DDML Object

DML-specific tidy method. Adds `ensemble_type` column labeling based on
the estimator's learner/ensemble configuration. Delegates to `tidy.ral`
for the base table computation.

## Usage

``` r
# S3 method for class 'ddml'
tidy(
  x,
  ensemble_idx = 1,
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

  A `ddml` object.

- ensemble_idx:

  Integer index of the ensemble type to report. Defaults to 1. Set to
  `NULL` for all.

- conf.int:

  Logical. Include confidence intervals? Default `FALSE`.

- conf.level:

  Confidence level. Default 0.95.

- type:

  Character. HC type. Default `"HC1"`.

- uniform:

  Logical. Uniform CIs? Default `FALSE`.

- bootstraps:

  Integer. Bootstrap draws. Default 999.

- ...:

  Currently unused.

## Value

A `data.frame` with columns `term`, `estimate`, `std.error`,
`statistic`, `p.value`, and `ensemble_type`.

## Examples

``` r
# \donttest{
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace")]
plm_fit = ddml_plm(y, D, X,
                learners = list(what = ols),
                sample_folds = 2, silent = TRUE)
tidy(plm_fit)
#>          term      estimate   std.error    statistic      p.value
#> 1          D1 -0.1545738155 0.014698010 -10.51664935 7.240272e-26
#> 2 (Intercept) -0.0001075169 0.006908269  -0.01556351 9.875826e-01
#>         ensemble_type
#> 1 single base learner
#> 2 single base learner
tidy(plm_fit, conf.int = TRUE)
#>          term      estimate   std.error    statistic      p.value
#> 1          D1 -0.1545738155 0.014698010 -10.51664935 7.240272e-26
#> 2 (Intercept) -0.0001075169 0.006908269  -0.01556351 9.875826e-01
#>         ensemble_type    conf.low   conf.high
#> 1 single base learner -0.18338139 -0.12576625
#> 2 single base learner -0.01364748  0.01343244
# }
```
