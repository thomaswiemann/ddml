# Cross-Fitted Predictions Using Stacking

Cross-fitted predictions using stacking.

## Usage

``` r
crosspred(
  y,
  X,
  learners,
  sample_folds = 10,
  ensemble_type = "average",
  cv_folds = 10,
  custom_ensemble_weights = NULL,
  cluster_variable = seq_along(y),
  subsamples = NULL,
  cv_subsamples = NULL,
  cv_subsamples_list = NULL,
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

- cv_folds:

  Number of folds used for cross-validation.

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

- cv_subsamples:

  List of lists, each corresponding to a subsample containing vectors
  with subsample indices for cross-validation.

- cv_subsamples_list:

  Deprecated; use `cv_subsamples` instead.

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

`crosspred` returns a list containing the following components:

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

- `cv_resid_byfold`:

  A list (length `sample_folds`) of inner cross-validation residual
  matrices used for ensemble weight estimation. `NULL` when a single
  learner is used.

- `auxiliary_fitted`:

  When `auxiliary_X` is not `NULL`, a list of matrices with additional
  predictions.

- `cf_fitted_bylearner`:

  A matrix of out-of-sample predictions, each column corresponding to a
  base learner (in chronological order).

- `cf_resid_bylearner`:

  A matrix of out-of-sample residuals (`y - cf_fitted_bylearner`), each
  column corresponding to a base learner.

- `auxiliary_fitted_bylearner`:

  When `auxiliary_X` is not `NULL`, a list of matrices with additional
  predictions for each learner.

## Details

`crosspred` implements the cross-fitting step of the Double/Debiased
Machine Learning procedure combined with stacking. It produces the
cross-fitted nuisance estimates \\\hat{\eta}(X_i)\\ used in the Neyman
orthogonal scores of all `ddml_*` estimators.

Let \\\\I_1, \ldots, I_S\\\\ be an \\S\\-fold partition of \\\\1,
\ldots, n\\\\, and denote the training set for fold \\s\\ by
\\\mathcal{T}\_s = \\1, \ldots, n\\ \setminus I_s\\. Given \\J\\ base
learners, the procedure operates on each cross-fitting fold \\s\\ in
three steps:

**Step 1 (Stacking weights).** Run \\K\\-fold cross-validation on
\\\mathcal{T}\_s\\ (via
[`crossval`](https://www.thomaswiemann.com/ddml/reference/crossval.md))
to estimate the MSPE of each base learner, and solve for fold-specific
stacking weights \\\hat{w}\_s = (\hat{w}\_{1,s}, \ldots,
\hat{w}\_{J,s})'\\.

**Step 2 (Fit).** Fit each base learner \\j\\ on the full training set
\\\mathcal{T}\_s\\, yielding \\\hat{f}\_{j,s}(\cdot)\\.

**Step 3 (Predict).** For each \\i \in I_s\\, compute the ensemble
cross-fitted prediction

\\\hat{\eta}(X_i) = \sum\_{j=1}^{J} \hat{w}\_{j,s}
\hat{f}\_{j,s}(X_i).\\

Since every observation belongs to exactly one fold, the result is a
complete \\n\\-vector of out-of-sample predictions. Crucially, both the
stacking weights \\\hat{w}\_s\\ and the base learner fits
\\\hat{f}\_{j,s}\\ depend only on \\\mathcal{T}\_s\\, which does not
contain observation \\i\\.

When a single learner is used (\\J = 1\\), no stacking or inner
cross-validation is performed: the learner is simply fitted on
\\\mathcal{T}\_s\\ and predictions are made for \\I_s\\.

## References

Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). "Model Averaging
and Double Machine Learning." Journal of Applied Econometrics, 40(3):
249-269.

Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2),
241-259.

## See also

Other utilities:
[`crossval()`](https://www.thomaswiemann.com/ddml/reference/crossval.md),
[`ddml()`](https://www.thomaswiemann.com/ddml/reference/ddml.md),
[`diagnostics()`](https://www.thomaswiemann.com/ddml/reference/diagnostics.md),
[`ensemble()`](https://www.thomaswiemann.com/ddml/reference/ensemble.md),
[`ensemble_weights()`](https://www.thomaswiemann.com/ddml/reference/ensemble_weights.md),
[`shortstacking()`](https://www.thomaswiemann.com/ddml/reference/shortstacking.md)

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
X = AE98[, c("morekids", "age","agefst","black","hisp","othrace","educ")]

# Compute cross-predictions using stacking with base learners ols and lasso.
#     Two stacking approaches are simultaneously computed: Equally
#     weighted (ensemble_type = "average") and MSPE-minimizing with weights
#     in the unit simplex (ensemble_type = "nnls1"). Predictions for each
#     learner are also calculated.
crosspred_res <- crosspred(y, X,
                           learners = list(list(what = ols),
                                           list(what = mdl_glmnet)),
                           ensemble_type = c("average",
                                             "nnls1",
                                             "singlebest"),
                           sample_folds = 2,
                           cv_folds = 2,
                           silent = TRUE)
dim(crosspred_res$cf_fitted) # = length(y) by length(ensemble_type)
#> [1] 5000    3
dim(crosspred_res$cf_fitted_bylearner) # = length(y) by length(learners)
#> [1] 5000    2
```
