# Estimator for the Partially Linear IV Model.

Estimator for the partially linear IV model.

## Usage

``` r
ddml_pliv(
  y,
  D,
  Z,
  X,
  learners,
  learners_DX = learners,
  learners_ZX = learners,
  sample_folds = 10,
  ensemble_type = "nnls",
  shortstack = FALSE,
  cv_folds = 10,
  custom_ensemble_weights = NULL,
  custom_ensemble_weights_DX = custom_ensemble_weights,
  custom_ensemble_weights_ZX = custom_ensemble_weights,
  cluster_variable = seq_along(y),
  subsamples = NULL,
  cv_subsamples_list = NULL,
  silent = FALSE
)
```

## Arguments

- y:

  The outcome variable.

- D:

  A matrix of endogenous variables.

- Z:

  A matrix of instruments.

- X:

  A (sparse) matrix of control variables.

- learners:

  May take one of two forms, depending on whether a single learner or
  stacking with multiple learners is used for estimation of the
  conditional expectation functions. If a single learner is used,
  `learners` is a list with two named elements:

  - `what` The base learner function. The function must be such that it
    predicts a named input `y` using a named input `X`.

  - `args` Optional arguments to be passed to `what`.

  If stacking with multiple learners is used, `learners` is a list of
  lists, each containing four named elements:

  - `fun` The base learner function. The function must be such that it
    predicts a named input `y` using a named input `X`.

  - `args` Optional arguments to be passed to `fun`.

  - `assign_X` An optional vector of column indices corresponding to
    control variables in `X` that are passed to the base learner.

  - `assign_Z` An optional vector of column indices corresponding to
    instruments in `Z` that are passed to the base learner.

  Omission of the `args` element results in default arguments being used
  in `fun`. Omission of `assign_X` (and/or `assign_Z`) results in
  inclusion of all variables in `X` (and/or `Z`).

- learners_DX, learners_ZX:

  Optional arguments to allow for different base learners for estimation
  of \\E\[D\|X\]\\, \\E\[Z\|X\]\\. Setup is identical to `learners`.

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

- shortstack:

  Boolean to use short-stacking.

- cv_folds:

  Number of folds used for cross-validation in ensemble construction.

- custom_ensemble_weights:

  A numerical matrix with user-specified ensemble weights. Each column
  corresponds to a custom ensemble specification, each row corresponds
  to a base learner in `learners` (in chronological order). Optional
  column names are used to name the estimation results corresponding the
  custom ensemble specification.

- custom_ensemble_weights_DX, custom_ensemble_weights_ZX:

  Optional arguments to allow for different custom ensemble weights for
  `learners_DX`,`learners_ZX`. Setup is identical to
  `custom_ensemble_weights`. Note: `custom_ensemble_weights` and
  `custom_ensemble_weights_DX`,`custom_ensemble_weights_ZX` must have
  the same number of columns.

- cluster_variable:

  A vector of cluster indices.

- subsamples:

  List of vectors with sample indices for cross-fitting.

- cv_subsamples_list:

  List of lists, each corresponding to a subsample containing vectors
  with subsample indices for cross-validation.

- silent:

  Boolean to silence estimation updates.

## Value

`ddml_pliv` returns an object of S3 class `ddml_pliv`. An object of
class `ddml_pliv` is a list containing the following components:

- `coef`:

  A vector with the \\\theta_0\\ estimates.

- `weights`:

  A list of matrices, providing the weight assigned to each base learner
  (in chronological order) by the ensemble procedure.

- `mspe`:

  A list of matrices, providing the MSPE of each base learner (in
  chronological order) computed by the cross-validation step in the
  ensemble construction.

- `iv_fit`:

  Object of class `ivreg` from the IV regression of \\Y -
  \hat{E}\[Y\vert X\]\\ on \\D - \hat{E}\[D\vert X\]\\ using \\Z -
  \hat{E}\[Z\vert X\]\\ as the instrument. See also
  [`AER::ivreg()`](https://rdrr.io/pkg/AER/man/ivreg.html) for details.

- `learners`,`learners_DX`,`learners_ZX`, `cluster_variable`,
  `subsamples`, `cv_subsamples_list`,`ensemble_type`:

  Pass-through of selected user-provided arguments. See above.

## Details

`ddml_pliv` provides a double/debiased machine learning estimator for
the parameter of interest \\\theta_0\\ in the partially linear IV model
given by

\\Y = \theta_0D + g_0(X) + U,\\

where \\(Y, D, X, Z, U)\\ is a random vector such that \\E\[Cov(U,
Z\vert X)\] = 0\\ and \\E\[Cov(D, Z\vert X)\] \neq 0\\, and \\g_0\\ is
an unknown nuisance function.

## References

Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). "Model Averaging
and Double Machine Learning." Journal of Applied Econometrics, 40(3):
249-269.

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). "Double/debiased machine learning for treatment and
structural parameters." The Econometrics Journal, 21(1), C1-C68.

Kleiber C, Zeileis A (2008). Applied Econometrics with R.
Springer-Verlag, New York.

Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2),
241-259.

## See also

[`summary.ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_plm.md),
[`AER::ivreg()`](https://rdrr.io/pkg/AER/man/ivreg.html)

Other ddml:
[`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md),
[`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md)

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
Z = AE98[, "samesex"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]

# Estimate the partially linear IV model using a single base learner, ridge.
pliv_fit <- ddml_pliv(y, D, Z, X,
                      learners = list(what = mdl_glmnet,
                                      args = list(alpha = 0)),
                      sample_folds = 2,
                      silent = TRUE)
summary(pliv_fit)
#> PLIV estimation results: 
#>  
#> , , single base learner
#> 
#>              Estimate Std. Error   t value Pr(>|t|)
#> (Intercept) -3.44e-07     0.0069 -4.99e-05    1.000
#> D_r         -2.35e-01     0.1893 -1.24e+00    0.214
#> 
```
