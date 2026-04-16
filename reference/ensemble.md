# Stacking Estimator Using Combinations of Base Learners

Computes an ensemble of learners based on the specified aggregation type
and computes cross-validated out-of-sample predictions to inform the
weights.

## Usage

``` r
ensemble(
  y,
  X,
  type = "average",
  learners,
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

  A character string indicating the type of ensemble to compute. Default
  is `"average"`.

- learners:

  A list of base learners. See
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

An object of class `ensemble` containing:

- `mdl_fits`:

  List of fitted base learners.

- `weights`:

  Computed ensemble weights.

- `learners`:

  The base learners used.

- `cv_results`:

  Cross-validation results if computed.

- `mean_y`:

  Mean of the outcome variable.

- `constant_y`:

  Boolean indicating if y is constant.

## See also

Other utilities:
[`crosspred()`](https://www.thomaswiemann.com/ddml/reference/crosspred.md),
[`crossval()`](https://www.thomaswiemann.com/ddml/reference/crossval.md),
[`ddml()`](https://www.thomaswiemann.com/ddml/reference/ddml.md),
[`diagnostics()`](https://www.thomaswiemann.com/ddml/reference/diagnostics.md),
[`ensemble_weights()`](https://www.thomaswiemann.com/ddml/reference/ensemble_weights.md),
[`shortstacking()`](https://www.thomaswiemann.com/ddml/reference/shortstacking.md)

## Examples

``` r
# \donttest{
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
X = AE98[, c("age","agefst","black","hisp","othrace")]

# Fit an ensemble of ols, lasso, and ridge
ens_fit = ensemble(y, X,
                   type = "nnls",
                   learners = list(list(what = ols),
                                  list(what = mdl_glmnet),
                                  list(what = mdl_glmnet,
                                       args = list(alpha = 0))),
                   cv_folds = 5,
                   silent = TRUE)
ens_fit$weights
#>           nnls
#> [1,] 0.0000000
#> [2,] 0.9994217
#> [3,] 0.0000000
predict(ens_fit, newdata = X)[1:5]
#> [1] 0.5202216 0.5770087 0.5512930 0.6168059 0.5401202
# }
```
