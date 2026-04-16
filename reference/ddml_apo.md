# Estimator for the Average Potential Outcome

Estimator for the average potential outcome, allowing for custom weights
\\\omega(X)\\.

## Usage

``` r
ddml_apo(
  y,
  D,
  X,
  d = 1,
  weights = NULL,
  learners,
  learners_DX = learners,
  sample_folds = 10,
  ensemble_type = "nnls",
  shortstack = FALSE,
  cv_folds = 10,
  custom_ensemble_weights = NULL,
  custom_ensemble_weights_DX = custom_ensemble_weights,
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

  The endogenous variable of interest. Can be discrete or continuous.

- X:

  A (sparse) matrix of control variables.

- d:

  The treatment level of interest. The default is `d = 1`.

- weights:

  A numeric vector of length `nobs` specifying the weights
  \\\omega(X)\\. If `weights = NULL` (the default), a vector of 1s is
  used, which estimates the Average Potential Outcome (APO).

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

- learners_DX:

  Optional argument to allow for different estimators of \\E\[D\|X\]\\.
  Setup is identical to `learners`.

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

- custom_ensemble_weights_DX:

  Optional argument to allow for different custom ensemble weights for
  `learners_DX`. Setup is identical to `custom_ensemble_weights`. Note:
  `custom_ensemble_weights` and `custom_ensemble_weights_DX` must have
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

  An optional list of sample split objects. For `ddml_apo`, this must be
  a list with elements `subsamples` and `cv_subsamples` (and optionally
  `subsamples_byD` and `cv_subsamples_byD` for stratified splitting).
  Typically obtained from a previous fit via `fit$splits`.

- save_crossval:

  Logical indicating whether to store the inner cross-validation
  residuals used for ensemble weight computation. Default `TRUE`. When
  `TRUE`, subsequent pass-through calls with data-driven ensembles
  (e.g., `"nnls"`) reproduce per-fold weights exactly. Set to `FALSE` to
  reduce object size at the cost of approximate weight recomputation.

- ...:

  Additional arguments passed to internal methods.

## Value

`ddml_apo` returns an object of S3 class `ddml_apo` and `ddml`. See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for the common output structure. Additional pass-through fields:
`learners`, `learners_DX`.

## Details

**Parameter of Interest:** `ddml_apo` provides a Double/Debiased Machine
Learning estimator for the average potential outcome. Under conditional
unconfoundedness and overlap, the parameter is identified by the
following reduced form conditional expectation:

\$\$\theta_0^{\textrm{APO}} = E\[\omega(X) E\[Y\|D=d, X\]\],\$\$

where \\W \equiv (Y, D, X)\\ is the observed random vector and
\\\omega(X)\\ is a known weighting function. If \\\omega(X) = 1\\, this
parameter corresponds to the average potential outcome at treatment
level \\d\\.

**Nuisance Parameters:** The nuisance parameters are \\\eta = (\ell,
r)\\ taking true values \\\ell_0(X) = E\[Y\|D=d, X\]\\ and \\r_0(X) =
\Pr(D=d\|X)\\.

**Neyman Orthogonal Score / Moment Equation:** The Neyman orthogonal
score is:

\$\$m(W; \theta, \eta) = \left( \frac{\mathbf{1}\\D=d\\ (Y -
\ell(X))}{r(X)} + \ell(X) \right) \omega(X) - \theta\$\$

**Jacobian:**

\$\$J = -1\$\$

See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for how the influence function and inference are derived from these
components.

## See also

Other ddml estimators:
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md),
[`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_attgt()`](https://www.thomaswiemann.com/ddml/reference/ddml_attgt.md),
[`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md),
[`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md),
[`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md),
[`ddml_policy()`](https://www.thomaswiemann.com/ddml/reference/ddml_policy.md)

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]

# Estimate the APO for d = 1 using a single base learner, ridge.
apo_fit <- ddml_apo(y, D, X,
                    learners = list(what = mdl_glmnet),
                    sample_folds = 2,
                    silent = TRUE)
summary(apo_fit)
#> DDML estimation: Average Potential Outcome 
#> Obs: 5000   Folds: 2
#> 
#>     Estimate Std. Error z value Pr(>|z|)    
#> APO   0.4352     0.0125    34.8   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
