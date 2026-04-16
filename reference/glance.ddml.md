# Glance at a DDML Object

DML-specific glance method. Includes DML fields like `sample_folds`,
`shortstack`, and `model_type`.

## Usage

``` r
# S3 method for class 'ddml'
glance(x, ...)
```

## Arguments

- x:

  A `ddml` object.

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
plm_fit = ddml_plm(y, D, X,
                learners = list(what = ols),
                sample_folds = 2, silent = TRUE)
glance(plm_fit)
#>   nobs sample_folds shortstack ensemble_type model_type         estimator_name
#> 1 5000            2      FALSE          nnls   ddml_plm Partially Linear Model
# }
```
