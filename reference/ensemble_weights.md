# Compute Stacking Weights for Base Learners

Computes the stacking weights for an ensemble of base learners using
cross-validated out-of-sample predictions.

## Usage

``` r
ensemble_weights(
  y,
  X,
  type = "average",
  learners = NULL,
  cv_folds = 5,
  cv_subsamples = NULL,
  cv_results = NULL,
  custom_weights = NULL,
  silent = FALSE
)
```

## Arguments

- y:

  The outcome variable.

- X:

  The feature matrix.

- type:

  A character string or vector indicating the type(s) of ensemble
  weights to compute. Default is `"average"`.

- learners:

  Optional list of base learners. Required when `cv_results` is not
  supplied (learners are needed to run cross-validation). When
  `cv_results` is supplied, `learners` may be omitted; the number of
  learners is inferred from the cross-validation residuals. See
  [`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
  for the full specification.

- cv_folds:

  Number of cross-validation folds.

- cv_subsamples:

  Optional list of subsamples for cross-validation.

- cv_results:

  Optional pre-computed cross-validation results.

- custom_weights:

  Optional custom weights matrix.

- silent:

  A boolean indicating whether to suppress progress messages.

## Value

A list containing:

- `weights`:

  A matrix of computed ensemble weights.

- `cv_results`:

  Cross-validation results used for computing weights.

## See also

Other utilities:
[`crosspred()`](https://www.thomaswiemann.com/ddml/reference/crosspred.md),
[`crossval()`](https://www.thomaswiemann.com/ddml/reference/crossval.md),
[`ddml()`](https://www.thomaswiemann.com/ddml/reference/ddml.md),
[`diagnostics()`](https://www.thomaswiemann.com/ddml/reference/diagnostics.md),
[`ensemble()`](https://www.thomaswiemann.com/ddml/reference/ensemble.md),
[`shortstacking()`](https://www.thomaswiemann.com/ddml/reference/shortstacking.md)

## Examples

``` r
# \donttest{
y = AE98[, "worked"]
X = AE98[, c("age","agefst","black","hisp","othrace")]

# Compute stacking weights via NNLS
ew = ensemble_weights(y, X,
                      type = "nnls",
                      learners = list(list(what = ols),
                                     list(what = mdl_glmnet)),
                      cv_folds = 5,
                      silent = TRUE)
ew$weights
#>          nnls
#> [1,] 0.000000
#> [2,] 0.999721
# }
```
