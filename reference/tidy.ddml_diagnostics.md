# Tidy Stacking Diagnostics

Returns a flat data.frame of per-learner stacking diagnostics for all
nuisance equations. Suitable for table creation with `kable()`, `gt()`,
or `modelsummary::datasummary()`.

## Usage

``` r
# S3 method for class 'ddml_diagnostics'
tidy(x, ...)
```

## Arguments

- x:

  An object of class `ddml_diagnostics`.

- ...:

  Currently unused.

## Value

A `data.frame` with columns `equation`, `learner`, `mspe`, `r2`,
`weight`, and optionally `cvc_pval`.

## Examples

``` r
# \donttest{
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace")]
learners = list(list(what = ols),
               list(what = mdl_glmnet))
plm_fit = ddml_plm(y, D, X,
                    learners = learners,
                    sample_folds = 2, silent = TRUE)
tidy(diagnostics(plm_fit, cvc = TRUE))
#>   equation   learner      mspe         r2 weight_nnls cvc_pval
#> 1      y_X learner_1 0.2439130 0.02151811   0.5786286    0.976
#> 2      y_X learner_2 0.2439839 0.02123374   0.4196157    0.028
#> 3      y_X      nnls 0.2439830 0.02123744          NA       NA
#> 4     D1_X learner_1 0.2182311 0.07742177   0.7225156    0.480
#> 5     D1_X learner_2 0.2182298 0.07742748   0.2736754    0.528
#> 6     D1_X      nnls 0.2182235 0.07745395          NA       NA
# }
```
