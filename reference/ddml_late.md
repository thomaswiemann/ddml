# Estimator for the Local Average Treatment Effect

Estimator for the local average treatment effect.

## Usage

``` r
ddml_late(
  y,
  D,
  Z,
  X,
  learners,
  learners_DXZ = learners,
  learners_ZX = learners,
  sample_folds = 10,
  ensemble_type = "nnls",
  shortstack = FALSE,
  cv_folds = 10,
  custom_ensemble_weights = NULL,
  custom_ensemble_weights_DXZ = custom_ensemble_weights,
  custom_ensemble_weights_ZX = custom_ensemble_weights,
  cluster_variable = seq_along(y),
  stratify = TRUE,
  trim = 0.01,
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

  Binary instrumental variable.

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

- learners_DXZ, learners_ZX:

  Optional arguments to allow for different base learners for estimation
  of \\E\[D \vert X, Z\]\\, \\E\[Z \vert X\]\\. Setup is identical to
  `learners`.

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

- custom_ensemble_weights_DXZ, custom_ensemble_weights_ZX:

  Optional arguments to allow for different custom ensemble weights for
  `learners_DXZ`,`learners_ZX`. Setup is identical to
  `custom_ensemble_weights`. Note: `custom_ensemble_weights` and
  `custom_ensemble_weights_DXZ`,`custom_ensemble_weights_ZX` must have
  the same number of columns.

- cluster_variable:

  A vector of cluster indices.

- stratify:

  Boolean for stratified cross-fitting: if `TRUE`, subsamples are
  constructed to be balanced across treatment levels.

- trim:

  Number in (0, 1) for trimming the estimated propensity scores at
  `trim` and `1-trim`.

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

  An optional list of sample split objects. For `ddml_late`, recommended
  keys are `subsamples`, `subsamples_byZ`, `cv_subsamples`, and
  `cv_subsamples_byZ`.

- save_crossval:

  Logical indicating whether to store the inner cross-validation
  residuals used for ensemble weight computation. Default `TRUE`. When
  `TRUE`, subsequent pass-through calls with data-driven ensembles
  (e.g., `"nnls"`) reproduce per-fold weights exactly. Set to `FALSE` to
  reduce object size at the cost of approximate weight recomputation.

- ...:

  Additional arguments passed to internal methods.

## Value

`ddml_late` returns an object of S3 class `ddml_late` and `ddml`. See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for the common output structure. Additional pass-through fields:
`learners`, `learners_DXZ`, `learners_ZX`.

## Details

**Parameter of Interest:** `ddml_late` provides a Double/Debiased
Machine Learning estimator for the local average treatment effect. Under
the standard instrumental variable assumptions (conditional
independence, exclusion restriction, relevance, and monotonicity) with a
binary instrument \\Z\\ and a binary treatment \\D\\, the parameter is
identified by the following reduced form parameter:

\$\$\theta_0^{\textrm{LATE}} = \frac{E\[E\[Y\|Z=1, X\] - E\[Y\|Z=0,
X\]\]}{E\[E\[D\|Z=1, X\] - E\[D\|Z=0, X\]\]}\$\$

where \\W \equiv (Y, D, X, Z)\\ is the observed random vector.

**Nuisance Parameters:** The nuisance parameters are \\\eta = (\ell_0,
\ell_1, r_0, r_1, p)\\ taking true values \\\ell\_{z,0}(X) = E\[Y\|Z=z,
X\]\\, \\r\_{z,0}(X) = E\[D\|Z=z, X\]\\, and \\p_0(X) = \Pr(Z=1\|X)\\.

**Neyman Orthogonal Score / Moment Equation:** The Neyman orthogonal
score is:

\$\$m(W; \theta, \eta) = \frac{Z(Y - \ell_1(X))}{p(X)} -
\frac{(1-Z)(Y-\ell_0(X))}{1-p(X)} + \ell_1(X) - \ell_0(X) -
\theta\left(\frac{Z(D - r_1(X))}{p(X)} -
\frac{(1-Z)(D-r_0(X))}{1-p(X)} + r_1(X) - r_0(X)\right)\$\$

**Jacobian:**

\$\$J = -E\[r_1(X) - r_0(X)\]\$\$

See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for how the influence function and inference are derived from these
components.

## References

Imbens G, Angrist J (1994). "Identification and Estimation of Local
Average Treatment Effects." Econometrica, 62(2), 467-475.

## See also

Other ddml estimators:
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md),
[`ddml_apo()`](https://www.thomaswiemann.com/ddml/reference/ddml_apo.md),
[`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_attgt()`](https://www.thomaswiemann.com/ddml/reference/ddml_attgt.md),
[`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md),
[`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md),
[`ddml_policy()`](https://www.thomaswiemann.com/ddml/reference/ddml_policy.md)

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
Z = AE98[, "samesex"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]

# Estimate the local average treatment effect using a single base learner,
#     ridge.
late_fit <- ddml_late(y, D, Z, X,
                      learners = list(what = mdl_glmnet,
                                      args = list(alpha = 0)),
                      sample_folds = 2,
                      silent = TRUE)
summary(late_fit)
#> DDML estimation: Local Average Treatment Effect 
#> Obs: 5000   Folds: 2
#> 
#>      Estimate Std. Error z value Pr(>|z|)
#> LATE   -0.231      0.190   -1.22     0.22

# \donttest{
# Estimate the local average treatment effect using short-stacking with base
#     learners ols, lasso, and ridge. We can also use custom_ensemble_weights
#     to estimate the LATE using every individual base learner.
weights_everylearner <- diag(1, 3)
colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
late_fit <- ddml_late(y, D, Z, X,
                      learners = list(list(what = ols),
                                      list(what = mdl_glmnet),
                                      list(what = mdl_glmnet,
                                           args = list(alpha = 0))),
                      ensemble_type = 'nnls',
                      custom_ensemble_weights = weights_everylearner,
                      shortstack = TRUE,
                      sample_folds = 2,
                      silent = TRUE)
summary(late_fit)
#> DDML estimation: Local Average Treatment Effect 
#> Obs: 5000   Folds: 2  Stacking: short-stack
#> 
#> Ensemble type: nnls
#>      Estimate Std. Error z value Pr(>|z|)
#> LATE   -0.232      0.184   -1.26     0.21
#> 
#> Ensemble type: mdl:ols
#>      Estimate Std. Error z value Pr(>|z|)
#> LATE   -0.241      0.184   -1.31     0.19
#> 
#> Ensemble type: mdl:lasso
#>      Estimate Std. Error z value Pr(>|z|)
#> LATE   -0.233      0.184   -1.27     0.21
#> 
#> Ensemble type: mdl:ridge
#>      Estimate Std. Error z value Pr(>|z|)
#> LATE   -0.234      0.185   -1.27     0.21
# }
```
