# Split a ddml_rep Object by Ensemble Type

Returns a named list of single-ensemble `ddml_rep` objects.

## Usage

``` r
# S3 method for class 'ddml_rep'
as.list(x, ...)
```

## Arguments

- x:

  A `ddml_rep` object.

- ...:

  Currently unused.

## Value

A named list of `ddml_rep` objects.

## Examples

``` r
# \donttest{
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace")]
reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
                      learners = list(what = ols),
                      resamples = 3,
                      sample_folds = 2,
                      silent = TRUE)
as.list(reps)
#> $nnls
#> DDML replicated fits: Partially Linear Model 
#>   Resamples: 3   Obs: 5000   Folds: 2 
#> 
#> Use summary() for aggregated inference.
#> Use x[[i]] to access individual fits.
#> 
# }
```
