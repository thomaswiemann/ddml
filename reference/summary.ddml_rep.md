# Summary for ddml_rep Objects

DML-specific summary override. Adds ensemble type labels, folds,
shortstack status to the base `ral_rep` summary.

## Usage

``` r
# S3 method for class 'ddml_rep'
summary(
  object,
  aggregation = c("median", "mean", "spectral"),
  type = "HC1",
  ...
)

# S3 method for class 'summary.ddml_rep'
print(x, digits = 3, ...)
```

## Arguments

- object:

  A `ddml_rep` object.

- aggregation:

  Character string: `"median"` (default), `"mean"`, or `"spectral"`.

- type:

  Character. HC type. Default `"HC1"`.

- ...:

  Currently unused.

- x:

  An object of class `summary.ddml_rep`.

- digits:

  Number of significant digits. Default 3.

## Value

An object of class `"summary.ddml_rep"`.

## Details

See
[`summary.ral_rep`](https://www.thomaswiemann.com/ddml/reference/summary.ral_rep.md)
for the aggregation formulas.

## References

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). "Double/debiased machine learning for treatment and
structural parameters." The Econometrics Journal, 21(1), C1-C68.

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
#> D1          -1.55e-01   1.47e-02  -10.50   <2e-16 ***
#> (Intercept) -8.61e-05   6.91e-03   -0.01     0.99    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(reps, aggregation = "mean")
#> DDML estimation: Partially Linear Model 
#> Obs: 5000   Folds: 2   Resamples: 3   Aggregation: mean
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> D1          -1.54e-01   1.48e-02   -10.4   <2e-16 ***
#> (Intercept)  1.62e-05   6.91e-03     0.0        1    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
