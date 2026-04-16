# Glance at a ddml_rep Object

DML-specific glance method. Includes DML fields.

## Usage

``` r
# S3 method for class 'ddml_rep'
glance(x, ...)
```

## Arguments

- x:

  A `ddml_rep` object.

- ...:

  Currently unused.

## Value

A one-row `data.frame`.

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
glance(reps)
#>   nobs sample_folds shortstack ensemble_type model_type         estimator_name
#> 1 5000            2      FALSE          nnls   ddml_plm Partially Linear Model
#>   nresamples
#> 1          3
# }
```
