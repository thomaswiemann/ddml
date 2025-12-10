# Print Methods for Treatment Effect Estimators.

Print methods for treatment effect estimators.

## Usage

``` r
# S3 method for class 'summary.ddml_fpliv'
print(x, digits = 3, ...)

# S3 method for class 'summary.ddml_pliv'
print(x, digits = 3, ...)

# S3 method for class 'summary.ddml_plm'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `summary.ddml_plm`, `summary.ddml_pliv`, and
  `summary.ddml_fpliv`, as returned by
  [`summary.ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_plm.md),
  [`summary.ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_plm.md),
  and
  [`summary.ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_plm.md),
  respectively.

- digits:

  Number of significant digits used for priniting.

- ...:

  Currently unused.

## Value

NULL.

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]

# Estimate the partially linear model using a single base learner, ridge.
plm_fit <- ddml_plm(y, D, X,
                    learners = list(what = mdl_glmnet,
                                    args = list(alpha = 0)),
                    sample_folds = 2,
                    silent = TRUE)
summary(plm_fit)
#> PLM estimation results: 
#>  
#> , , single base learner
#> 
#>              Estimate Std. Error t value Pr(>|t|)
#> (Intercept) -0.000385    0.00687  -0.056 9.55e-01
#> D_r         -0.147853    0.01468 -10.071 7.43e-24
#> 
```
