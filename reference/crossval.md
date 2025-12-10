# Estimator of the Mean Squared Prediction Error using Cross-Validation.

Estimator of the mean squared prediction error of different learners
using cross-validation.

## Usage

``` r
crossval(
  y,
  X,
  Z = NULL,
  learners,
  cv_folds = 5,
  cv_subsamples = NULL,
  silent = FALSE,
  progress = NULL
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

  `learners` is a list of lists, each containing four named elements:

  - `fun` The base learner function. The function must be such that it
    predicts a named input `y` using a named input `X`.

  - `args` Optional arguments to be passed to `fun`.

  - `assign_X` An optional vector of column indices corresponding to
    variables in `X` that are passed to the base learner.

  - `assign_Z` An optional vector of column indices corresponding to
    variables in `Z` that are passed to the base learner.

  Omission of the `args` element results in default arguments being used
  in `fun`. Omission of `assign_X` (and/or `assign_Z`) results in
  inclusion of all predictive variables in `X` (and/or `Z`).

- cv_folds:

  Number of folds used for cross-validation.

- cv_subsamples:

  List of vectors with sample indices for cross-validation.

- silent:

  Boolean to silence estimation updates.

- progress:

  String to print before learner and cv fold progress.

## Value

`crossval` returns a list containing the following components:

- `mspe`:

  A vector of MSPE estimates, each corresponding to a base learners (in
  chronological order).

- `oos_resid`:

  A matrix of out-of-sample prediction errors, each column corresponding
  to a base learners (in chronological order).

- `cv_subsamples`:

  Pass-through of `cv_subsamples`. See above.

## See also

Other utilities:
[`crosspred()`](https://www.thomaswiemann.com/ddml/reference/crosspred.md),
[`shortstacking()`](https://www.thomaswiemann.com/ddml/reference/shortstacking.md)

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
X = AE98[, c("morekids", "age","agefst","black","hisp","othrace","educ")]

# Compare ols, lasso, and ridge using 4-fold cross-validation
cv_res <- crossval(y, X,
                   learners = list(list(fun = ols),
                                   list(fun = mdl_glmnet),
                                   list(fun = mdl_glmnet,
                                        args = list(alpha = 0))),
                   cv_folds = 4,
                   silent = TRUE)
cv_res$mspe
#> [1] 0.2365375 0.2365378 0.2365485
```
