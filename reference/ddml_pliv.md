# Estimator for the Partially Linear IV Coefficient

Estimator for the partially linear IV coefficient.

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
  silent = FALSE,
  parallel = NULL,
  fitted = NULL,
  splits = NULL,
  save_crossval = TRUE,
  ...
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
  lists, each containing three named elements:

  - `what` The base learner function. The function must be such that it
    predicts a named input `y` using a named input `X`.

  - `args` Optional arguments to be passed to `what`.

  - `assign_X` An optional vector of column indices corresponding to
    control variables in `X` that are passed to the base learner.

  Omission of the `args` element results in default arguments being used
  in `what`. Omission of `assign_X` results in inclusion of all
  variables in `X`.

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

- fitted:

  An optional named list of per-equation cross-fitted predictions,
  typically obtained from a previous fit via `fit$fitted`. When supplied
  (together with `splits`), base learners are not re-fitted; only
  ensemble weights are recomputed. This allows fast re-estimation with a
  different `ensemble_type`. See
  [`ddml_plm`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md)
  for an example.

- splits:

  An optional list of sample split objects, typically obtained from a
  previous fit via `fit$splits`. Must be supplied when `fitted` is
  provided. Can also be used standalone to provide pre-computed sample
  folds.

- save_crossval:

  Logical indicating whether to store the inner cross-validation
  residuals used for ensemble weight computation. Default `TRUE`. When
  `TRUE`, subsequent pass-through calls with data-driven ensembles
  (e.g., `"nnls"`) reproduce per-fold weights exactly. Set to `FALSE` to
  reduce object size at the cost of approximate weight recomputation.

- ...:

  Additional arguments passed to internal methods.

## Value

`ddml_pliv` returns an object of S3 class `ddml_pliv` and `ddml`. See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for the common output structure. Additional pass-through fields:
`learners`, `learners_DX`, `learners_ZX`.

## Details

**Parameter of Interest:** `ddml_pliv` provides a Double/Debiased
Machine Learning estimator for the partially linear instrumental
variable (IV) coefficient \\\theta_0\\, defined by the partially linear
IV model:

\$\$Y = \theta_0 D + g_0(X) + \varepsilon, \quad E\[Z\varepsilon\] = 0,
\quad E\[\varepsilon\|X\] = 0,\$\$

where \\W \equiv (Y, D, X, Z, \varepsilon)\\ is a random vector such
that \\E\[Cov(D, Z\|X)\] \neq 0\\, and \\g_0(X)\\ is an unknown nuisance
function.

**Neyman Orthogonal Score:** The Neyman orthogonal score is:

\$\$m(W; \theta, \eta) = \[(Y - \ell(X)) - \theta(D - r_D(X))\](Z -
r_Z(X))\$\$

where the nuisance parameters are \\\eta = (\ell, r_D, r_Z)\\ taking
true values \\\ell_0(X) = E\[Y\|X\]\\, \\r\_{D,0}(X) = E\[D\|X\]\\, and
\\r\_{Z,0}(X) = E\[Z\|X\]\\.

**Jacobian:**

\$\$J = -E\[(D - r_D(X))(Z - r_Z(X))^\top\]\$\$

See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for how the influence function and inference are derived from these
components.

## See also

Other ddml estimators:
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md),
[`ddml_apo()`](https://www.thomaswiemann.com/ddml/reference/ddml_apo.md),
[`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_attgt()`](https://www.thomaswiemann.com/ddml/reference/ddml_attgt.md),
[`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md),
[`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md),
[`ddml_policy()`](https://www.thomaswiemann.com/ddml/reference/ddml_policy.md)

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
#> DDML estimation: Partially Linear IV Model 
#> Obs: 5000   Folds: 2
#> 
#>             Estimate Std. Error z value Pr(>|z|)
#> D1          -0.21726    0.18642   -1.17     0.24
#> (Intercept) -0.00080    0.00691   -0.12     0.91
```
