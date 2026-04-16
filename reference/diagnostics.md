# Stacking Diagnostics for DDML Estimators

Computes per-learner diagnostics including MSPE, R-squared, ensemble
weights, and optionally cross-validation comparison (CVC) p-values for
each nuisance equation.

## Usage

``` r
diagnostics(object, cvc = FALSE, bootnum = 500, ...)
```

## Arguments

- object:

  An object of class `ddml`.

- cvc:

  Logical. Compute CVC p-values via multiplier bootstrap? Default
  `FALSE`. CVC tests whether each learner is significantly outperformed
  by the others.

- bootnum:

  Number of bootstrap replications for CVC. Default 500. Ignored when
  `cvc = FALSE`.

- ...:

  Currently unused.

## Value

An object of class `ddml_diagnostics` containing per-equation learner
diagnostics. Use [`print()`](https://rdrr.io/r/base/print.html) for
formatted output or
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) for a flat
data.frame.

## References

Lei J (2020). "Cross-Validation With Confidence." Journal of the
American Statistical Association, 115(532), 1978-1997.

## See also

Other utilities:
[`crosspred()`](https://www.thomaswiemann.com/ddml/reference/crosspred.md),
[`crossval()`](https://www.thomaswiemann.com/ddml/reference/crossval.md),
[`ddml()`](https://www.thomaswiemann.com/ddml/reference/ddml.md),
[`ensemble()`](https://www.thomaswiemann.com/ddml/reference/ensemble.md),
[`ensemble_weights()`](https://www.thomaswiemann.com/ddml/reference/ensemble_weights.md),
[`shortstacking()`](https://www.thomaswiemann.com/ddml/reference/shortstacking.md)

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
diagnostics(plm_fit, cvc = TRUE)
#> Stacking diagnostics: Partially Linear Model 
#> Obs: 5000 
#> 
#>   y_X:
#>    learner   mspe     r2 weight_nnls cvc_pval
#>  learner_1 0.2440 0.0211      0.4987    0.494
#>  learner_2 0.2440 0.0211      0.5001    0.458
#>       nnls 0.2439 0.0214          NA       NA
#> 
#>   D1_X:
#>    learner  mspe     r2 weight_nnls cvc_pval
#>  learner_1 0.218 0.0785      0.9968    0.478
#>  learner_2 0.218 0.0785      0.0000    0.488
#>       nnls 0.218 0.0785          NA       NA
#> 
tidy(diagnostics(plm_fit, cvc = TRUE))
#>   equation   learner      mspe         r2 weight_nnls cvc_pval
#> 1      y_X learner_1 0.2440140 0.02111301   0.4987039    0.454
#> 2      y_X learner_2 0.2440149 0.02110934   0.5000849    0.506
#> 3      y_X      nnls 0.2439484 0.02137638          NA       NA
#> 4     D1_X learner_1 0.2179779 0.07849230   0.9967816    0.480
#> 5     D1_X learner_2 0.2179783 0.07849050   0.0000000    0.480
#> 6     D1_X      nnls 0.2179782 0.07849105          NA       NA
# }
```
