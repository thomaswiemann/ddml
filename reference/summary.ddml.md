# Summary for DDML Estimators

Computes a coefficient table with estimates, standard errors, z-values,
and p-values for all ensemble types. Standard errors are based on a
heteroskedasticity-robust sandwich variance; see
[`vcov.ral`](https://www.thomaswiemann.com/ddml/reference/vcov.ral.md)
for the HC0/HC1/HC3 formulas.

## Usage

``` r
# S3 method for class 'ddml'
summary(object, type = "HC1", ...)

# S3 method for class 'summary.ddml'
print(x, digits = 3, ...)
```

## Arguments

- object:

  An object of class `ddml`.

- type:

  Character. HC type (`"HC0"`, `"HC1"`, or `"HC3"`). Default `"HC1"`.

- ...:

  Currently unused.

- x:

  An object of class `summary.ddml`.

- digits:

  Number of significant digits. Default 3.

## Value

An object of class `summary.ddml` with:

- `coefficients`:

  A 3-dimensional array (\\p \times 4 \times\\ nensb) of estimates,
  standard errors, z-values, and p-values.

- `type`:

  The HC type used.

- `nobs`:

  Number of observations.

- `sample_folds`:

  Number of cross-fitting folds.

- `ensemble_type`:

  Ensemble type labels.

## See also

[`vcov.ral`](https://www.thomaswiemann.com/ddml/reference/vcov.ral.md)

Other ddml inference:
[`lincom()`](https://www.thomaswiemann.com/ddml/reference/lincom.md)

## Examples

``` r
# \donttest{
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace")]
plm_fit = ddml_plm(y, D, X,
                learners = list(what = ols),
                sample_folds = 2, silent = TRUE)
summary(plm_fit)
#> DDML estimation: Partially Linear Model 
#> Obs: 5000   Folds: 2
#> 
#>             Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.15376    0.01471  -10.45   <2e-16 ***
#> (Intercept) -0.00011    0.00690   -0.02     0.99    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(plm_fit, type = "HC3")
#> DDML estimation: Partially Linear Model 
#> Obs: 5000   Folds: 2  SE: HC3
#> 
#>             Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.15376    0.01472  -10.45   <2e-16 ***
#> (Intercept) -0.00011    0.00690   -0.02     0.99    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
