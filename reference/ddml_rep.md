# Construct a Multi-Resample DDML Object

Validates a list of `ddml` fits and stamps class `"ddml_rep"` for
multi-resample aggregation.

## Usage

``` r
ddml_rep(fits)

# S3 method for class 'ddml_rep'
print(x, ...)
```

## Arguments

- fits:

  A list of at least 2 objects inheriting from class `"ddml"`. All fits
  must share the same primary class, coefficient names, ensemble type,
  and number of observations.

- x:

  A `ddml_rep` object.

- ...:

  Currently unused.

## Value

An object of class `c("ddml_rep", "ral_rep")` with fields:

- fits:

  List of `ddml` objects.

- nresamples:

  Number of resamples.

- model_type:

  Primary class of the fits.

- coef_names:

  Coefficient names.

- ensemble_type:

  Ensemble types.

- nobs:

  Number of observations.

- sample_folds:

  Number of cross-fitting folds.

- shortstack:

  Logical, whether short-stacking was used.

## See also

[`ddml_replicate()`](https://www.thomaswiemann.com/ddml/reference/ddml_replicate.md)

Other ddml replication:
[`ddml_replicate()`](https://www.thomaswiemann.com/ddml/reference/ddml_replicate.md)

## Examples

``` r
# \donttest{
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace")]
fits = lapply(1:3, function(r) {
  ddml_plm(y, D, X,
           learners = list(what = ols),
           sample_folds = 2, silent = TRUE)
})
reps = ddml_rep(fits)
summary(reps)
#> DDML estimation: Partially Linear Model 
#> Obs: 5000   Folds: 2   Resamples: 3   Aggregation: median
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> D1          -1.54e-01   1.48e-02  -10.43   <2e-16 ***
#> (Intercept) -9.98e-05   6.91e-03   -0.01     0.99    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
