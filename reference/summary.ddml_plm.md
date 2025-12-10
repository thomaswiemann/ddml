# Inference Methods for Partially Linear Estimators.

Inference methods for partially linear estimators. Simple wrapper for
[`sandwich::vcovHC()`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html)
and
[`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html).
Default standard errors are heteroskedasiticty-robust. If the `ddml`
estimator was computed using a `cluster_variable`, the standard errors
are also cluster-robust by default.

## Usage

``` r
# S3 method for class 'ddml_fpliv'
summary(object, ...)

# S3 method for class 'ddml_pliv'
summary(object, ...)

# S3 method for class 'ddml_plm'
summary(object, ...)
```

## Arguments

- object:

  An object of class `ddml_plm`, `ddml_pliv`, or `ddml_fpliv` as fitted
  by
  [`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md),
  [`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md),
  and
  [`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
  respectively.

- ...:

  Additional arguments passed to `vcovHC` and `vcovCL`. See
  [`sandwich::vcovHC()`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html)
  and
  [`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)
  for a complete list of arguments.

## Value

An array with inference results for each `ensemble_type`.

## References

Zeileis A (2004). "Econometric Computing with HC and HAC Covariance
Matrix Estimators.” Journal of Statistical Software, 11(10), 1-17.

Zeileis A (2006). “Object-Oriented Computation of Sandwich Estimators.”
Journal of Statistical Software, 16(9), 1-16.

Zeileis A, Köll S, Graham N (2020). “Various Versatile Variances: An
Object-Oriented Implementation of Clustered Covariances in R.” Journal
of Statistical Software, 95(1), 1-36.

## See also

[`sandwich::vcovHC()`](https://sandwich.R-Forge.R-project.org/reference/vcovHC.html),
[`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)

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
#> (Intercept) -0.000812    0.00688  -0.118 9.06e-01
#> D_r         -0.147777    0.01471 -10.044 9.76e-24
#> 
```
