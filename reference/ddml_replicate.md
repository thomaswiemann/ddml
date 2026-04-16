# Replicate a DDML Estimator Across Multiple Resamples

Convenience wrapper that calls a `ddml_*` estimator function multiple
times with independent sample splits and returns a `ddml_rep` object for
aggregated inference.

## Usage

``` r
ddml_replicate(fn, ..., resamples = 5, silent = FALSE)
```

## Arguments

- fn:

  A `ddml_*` estimator function (e.g., `ddml_plm`).

- ...:

  Arguments passed to `fn`.

- resamples:

  Integer number of independent resamples. Must be \>= 2. Default 5.

- silent:

  Logical. If `TRUE`, suppresses all output at both the resample level
  and within each estimator call. Default `FALSE`.

## Value

An object of class `"ddml_rep"`.

## See also

[`ddml_rep()`](https://www.thomaswiemann.com/ddml/reference/ddml_rep.md)

Other ddml replication:
[`ddml_rep()`](https://www.thomaswiemann.com/ddml/reference/ddml_rep.md)

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
summary(reps)
#> DDML estimation: Partially Linear Model 
#> Obs: 5000   Folds: 2   Resamples: 3   Aggregation: median
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.155060   0.014738  -10.52   <2e-16 ***
#> (Intercept) -0.000106   0.006947   -0.02     0.99    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
