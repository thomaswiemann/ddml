# Estimator of the Mean Squared Prediction Error Using Cross-Validation

Estimator of the mean squared prediction error of different learners
using cross-validation.

## Usage

``` r
crossval(
  y,
  X,
  learners,
  cv_folds = 10,
  cluster_variable = seq_along(y),
  cv_subsamples = NULL,
  silent = FALSE,
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

- cv_folds:

  Number of folds used for cross-validation.

- cluster_variable:

  A vector of cluster indices.

- cv_subsamples:

  List of vectors with sample indices for cross-validation.

- silent:

  Boolean to silence estimation updates.

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

`crossval` returns a list containing the following components:

- `mspe`:

  A vector of MSPE estimates, each corresponding to a base learner (in
  chronological order).

- `r2`:

  A vector of cross-validated \\R^2\\ values, each corresponding to a
  base learner (in chronological order).

- `cv_resid`:

  A matrix of out-of-sample residuals, each column corresponding to a
  base learner (in chronological order).

- `cv_subsamples`:

  Pass-through of `cv_subsamples`. See above.

## Details

`crossval` estimates the mean squared prediction error (MSPE) of \\J\\
base learners via \\K\\-fold cross-validation. It is the inner workhorse
of the stacking machinery used by
[`ensemble_weights`](https://www.thomaswiemann.com/ddml/reference/ensemble_weights.md)
to determine ensemble weights.

Given a generic conditional expectation function \\f_0(\cdot)\\ (e.g.,
\\E\[Y\vert X\]\\, \\E\[D\vert X\]\\), let \\\\I_1, \ldots, I_K\\\\ be a
\\K\\-fold partition of \\\\1, \ldots, n\\\\ and let
\\\hat{f}\_j^{(-k)}\\ denote learner \\j\\ trained on all observations
outside fold \\I_k\\. The out-of-sample residual for observation \\i \in
I_k\\ is

\\\hat{e}\_{i,j} = y_i - \hat{f}\_j^{(-k)}(X_i).\\

Since every observation belongs to exactly one fold, this yields a
complete \\n \times J\\ residual matrix. The cross-validated MSPE for
learner \\j\\ is

\\\widehat{\textrm{MSPE}}\_j = n^{-1} \sum\_{i=1}^{n}
\hat{e}\_{i,j}^2,\\

and the cross-validated \\R^2\\ is

\\\hat{R}^2_j = 1 - \widehat{\textrm{MSPE}}\_j \\/\\ \hat{\sigma}^2_y,\\

where \\\hat{\sigma}^2_y\\ is the sample variance of \\y\\.

## See also

Other utilities:
[`crosspred()`](https://www.thomaswiemann.com/ddml/reference/crosspred.md),
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

# Compare ols, lasso, and ridge using 4-fold cross-validation
cv_res <- crossval(y, X,
                   learners = list(list(what = ols),
                                   list(what = mdl_glmnet),
                                   list(what = mdl_glmnet,
                                        args = list(alpha = 0))),
                   cv_folds = 4,
                   silent = TRUE)
cv_res$mspe
#> [1] 0.2363514 0.2363506 0.2363641
```
