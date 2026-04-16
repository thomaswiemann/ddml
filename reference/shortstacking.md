# Predictions using Short-Stacking

Predictions using short-stacking.

## Usage

``` r
shortstacking(
  y,
  X,
  learners,
  sample_folds = 2,
  ensemble_type = "average",
  custom_ensemble_weights = NULL,
  cluster_variable = seq_along(y),
  subsamples = NULL,
  silent = FALSE,
  auxiliary_X = NULL,
  parallel = NULL
)
```

## Arguments

- y:

  The outcome variable.

- X:

  A (sparse) matrix of predictive variables.

- learners:

  `learners` is a list of lists, each containing three named elements:

  - `what` The base learner function. The function must be such that it
    predicts a named input `y` using a named input `X`.

  - `args` Optional arguments to be passed to `what`.

  - `assign_X` An optional vector of column indices corresponding to
    variables in `X` that are passed to the base learner.

  Omission of the `args` element results in default arguments being used
  in `what`. Omission of `assign_X` results in inclusion of all
  predictive variables in `X`.

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

- cluster_variable:

  A vector of cluster indices.

- subsamples:

  List of vectors with sample indices for cross-fitting.

- silent:

  Boolean to silence estimation updates.

- auxiliary_X:

  An optional list of matrices of length `sample_folds`, each containing
  additional observations to calculate predictions for.

- parallel:

  An optional named list with parallel processing options. When `NULL`
  (the default), computation is sequential. Supported fields:

  `cores`

  :   Number of cores to use.

  `export`

  :   Character vector of object names to export to parallel workers
      (for custom learners that reference global objects).

  `packages`

  :   Character vector of additional package names to load on workers
      (for custom learners that use packages not imported by `ddml`).

## Value

`shortstack` returns a list containing the following components:

- `cf_fitted`:

  A matrix of out-of-sample predictions, each column corresponding to an
  ensemble type (in chronological order).

- `weights`:

  An array, providing the weight assigned to each base learner (in
  chronological order) by the ensemble procedures.

- `mspe`:

  A numeric vector of per-learner out-of-sample MSPEs, computed from
  cross-fitted residuals.

- `r2`:

  A numeric vector of per-learner out-of-sample R-squared values.

- `auxiliary_fitted`:

  When `auxiliary_X` is not `NULL`, a list of matrices with additional
  predictions.

- `cf_fitted_bylearner`:

  A matrix of out-of-sample predictions, each column corresponding to a
  base learner (in chronological order).

- `cf_resid_bylearner`:

  A matrix of per-learner out-of-sample residuals used for weight
  estimation.

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
[`crossval()`](https://www.thomaswiemann.com/ddml/reference/crossval.md),
[`ddml()`](https://www.thomaswiemann.com/ddml/reference/ddml.md),
[`diagnostics()`](https://www.thomaswiemann.com/ddml/reference/diagnostics.md),
[`ensemble()`](https://www.thomaswiemann.com/ddml/reference/ensemble.md),
[`ensemble_weights()`](https://www.thomaswiemann.com/ddml/reference/ensemble_weights.md)

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
                                learners = list(list(what = ols),
                                                list(what = mdl_glmnet)),
                                ensemble_type = c("average",
                                                  "nnls1",
                                                  "singlebest"),
                                sample_folds = 2,
                                silent = TRUE)
dim(shortstack_res$cf_fitted) # = length(y) by length(ensemble_type)
#> [1] 5000    3
dim(shortstack_res$cf_fitted_bylearner) # = length(y) by length(learners)
#> [1] 5000    2
```
