# Predictions using Short-Stacking.

Predictions using short-stacking.

## Usage

``` r
shortstacking(
  y,
  X,
  Z = NULL,
  learners,
  sample_folds = 2,
  ensemble_type = "average",
  custom_ensemble_weights = NULL,
  compute_insample_predictions = FALSE,
  subsamples = NULL,
  silent = FALSE,
  progress = NULL,
  auxiliary_X = NULL,
  shortstack_y = y
)
```

## Arguments

- y:

  The outcome variable.

- X:

  A (sparse) matrix of predictive variables.

- Z:

  Optional additional (sparse) matrix of predictive variables.

- learners:

  May take one of two forms, depending on whether a single learner or
  stacking with multiple learners is used for estimation of the
  predictor. If a single learner is used, `learners` is a list with two
  named elements:

  - `what` The base learner function. The function must be such that it
    predicts a named input `y` using a named input `X`.

  - `args` Optional arguments to be passed to `what`.

  If stacking with multiple learners is used, `learners` is a list of
  lists, each containing four named elements:

  - `fun` The base learner function. The function must be such that it
    predicts a named input `y` using a named input `X`.

  - `args` Optional arguments to be passed to `fun`.

  - `assign_X` An optional vector of column indices corresponding to
    predictive variables in `X` that are passed to the base learner.

  - `assign_Z` An optional vector of column indices corresponding to
    predictive in `Z` that are passed to the base learner.

  Omission of the `args` element results in default arguments being used
  in `fun`. Omission of `assign_X` (and/or `assign_Z`) results in
  inclusion of all variables in `X` (and/or `Z`).

- sample_folds:

  Number of cross-fitting folds.

- ensemble_type:

  Ensemble method to combine base learners into final estimate of the
  conditional expectation functions. Possible values are:

  - `"nnls"` Non-negative least squares.

  - `"nnls1"` Non-negative least squares with the constraint that all
    weights sum to one.

  - `"singlebest"` Select base learner with minimum MSPE.

  - `"ols"` Ordinary least squares.

  - `"average"` Simple average over base learners.

  Multiple ensemble types may be passed as a vector of strings.

- custom_ensemble_weights:

  A numerical matrix with user-specified ensemble weights. Each column
  corresponds to a custom ensemble specification, each row corresponds
  to a base learner in `learners` (in chronological order). Optional
  column names are used to name the estimation results corresponding the
  custom ensemble specification.

- compute_insample_predictions:

  Indicator equal to 1 if in-sample predictions should also be computed.

- subsamples:

  List of vectors with sample indices for cross-fitting.

- silent:

  Boolean to silence estimation updates.

- progress:

  String to print before learner and cv fold progress.

- auxiliary_X:

  An optional list of matrices of length `sample_folds`, each containing
  additional observations to calculate predictions for.

- shortstack_y:

  Optional vector of the outcome variable to form short-stacking
  predictions for. Base learners are always trained on `y`.

## Value

`shortstack` returns a list containing the following components:

- `oos_fitted`:

  A matrix of out-of-sample predictions, each column corresponding to an
  ensemble type (in chronological order).

- `weights`:

  An array, providing the weight assigned to each base learner (in
  chronological order) by the ensemble procedures.

- `is_fitted`:

  When `compute_insample_predictions = T`. a list of matrices with
  in-sample predictions by sample fold.

- `auxiliary_fitted`:

  When `auxiliary_X` is not `NULL`, a list of matrices with additional
  predictions.

- `oos_fitted_bylearner`:

  A matrix of out-of-sample predictions, each column corresponding to a
  base learner (in chronological order).

- `is_fitted_bylearner`:

  When `compute_insample_predictions = T`, a list of matrices with
  in-sample predictions by sample fold.

- `auxiliary_fitted_bylearner`:

  When `auxiliary_X` is not `NULL`, a list of matrices with additional
  predictions for each learner.

Note that unlike `crosspred`, `shortstack` always computes out-of-sample
predictions for each base learner (at no additional computational cost).

## References

Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). "Model Averaging
and Double Machine Learning." Journal of Applied Econometrics, 40(3):
249-269.

Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2),
241-259.

## See also

Other utilities:
[`crosspred()`](https://www.thomaswiemann.com/ddml/reference/crosspred.md),
[`crossval()`](https://www.thomaswiemann.com/ddml/reference/crossval.md)

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
X = AE98[, c("morekids", "age","agefst","black","hisp","othrace","educ")]

# Compute predictions using shortstacking with base learners ols and lasso.
#     Two stacking approaches are simultaneously computed: Equally
#     weighted (ensemble_type = "average") and MSPE-minimizing with weights
#     in the unit simplex (ensemble_type = "nnls1"). Predictions for each
#     learner are also calculated.
shortstack_res <- shortstacking(y, X,
                                learners = list(list(fun = ols),
                                                list(fun = mdl_glmnet)),
                                ensemble_type = c("average",
                                                  "nnls1",
                                                  "singlebest"),
                                sample_folds = 2,
                                silent = TRUE)
dim(shortstack_res$oos_fitted) # = length(y) by length(ensemble_type)
#> [1] 5000    3
dim(shortstack_res$oos_fitted_bylearner) # = length(y) by length(learners)
#> [1] 5000    2
```
