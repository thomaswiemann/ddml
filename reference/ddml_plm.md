# Estimator for the Partially Linear Regression Coefficient

Estimator for the partially linear regression coefficient.

## Usage

``` r
ddml_plm(
  y,
  D,
  X,
  learners,
  learners_DX = learners,
  sample_folds = 10,
  ensemble_type = "nnls",
  shortstack = FALSE,
  cv_folds = 10,
  custom_ensemble_weights = NULL,
  custom_ensemble_weights_DX = custom_ensemble_weights,
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
  different `ensemble_type`. See `ddml_plm` for an example.

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

`ddml_plm` returns an object of S3 class `ddml_plm` and `ddml`. See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for the common output structure. Additional pass-through fields:
`learners`, `learners_DX`.

## Details

**Parameter of Interest:** `ddml_plm` provides a Double/Debiased Machine
Learning estimator for the partially linear regression coefficient
\\\theta_0\\, defined by the partially linear regression model:

\$\$Y = \theta_0 D + g_0(X) + \varepsilon, \quad E\[D\varepsilon\] = 0,
\quad E\[\varepsilon\|X\] = 0,\$\$

where \\W\equiv(Y, D, X, \varepsilon)\\ is a random vector such that
\\E\[Var(D\|X)\] \neq 0\\, and \\g_0(X)\\ is an unknown nuisance
function.

**Neyman Orthogonal Score:** The Neyman orthogonal score is:

\$\$m(W; \theta, \eta) = \[(Y - \ell(X)) - \theta(D - r(X))\](D -
r(X))\$\$

where the nuisance parameters are \\\eta = (\ell, r)\\ taking true
values \\\ell_0(X) = E\[Y\|X\]\\ and \\r_0(X) = E\[D\|X\]\\.

**Jacobian:**

\$\$J = -E\[(D - r_0(X))(D - r_0(X))^\top\]\$\$

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
[`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md),
[`ddml_policy()`](https://www.thomaswiemann.com/ddml/reference/ddml_policy.md)

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]

# Estimate the partially linear model using a single base learner, ridge.
plm_fit <- ddml_plm(y, D, X,
                    learners = list(what = mdl_glmnet,
                                    args = list(alpha = 0)),
                    sample_folds = 2,
                    silent = TRUE)
summary(plm_fit)
#> DDML estimation: Partially Linear Model 
#> Obs: 5000   Folds: 2
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.146551   0.014646  -10.01   <2e-16 ***
#> (Intercept)  0.000292   0.006872    0.04     0.97    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Estimate the partially linear model using short-stacking with base learners
#     ols, lasso, and ridge. We can also use custom_ensemble_weights
#     to estimate the ATE using every individual base learner.
weights_everylearner <- diag(1, 3)
colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
plm_fit <- ddml_plm(y, D, X,
                    learners = list(list(what = ols),
                                    list(what = mdl_glmnet),
                                    list(what = mdl_glmnet,
                                         args = list(alpha = 0))),
                    ensemble_type = 'nnls',
                    custom_ensemble_weights = weights_everylearner,
                    shortstack = TRUE,
                    sample_folds = 2,
                    silent = TRUE)
summary(plm_fit)
#> DDML estimation: Partially Linear Model 
#> Obs: 5000   Folds: 2  Stacking: short-stack
#> 
#> Ensemble type: nnls
#>             Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.14851    0.01468  -10.12   <2e-16 ***
#> (Intercept)  0.00166    0.00688    0.24     0.81    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Ensemble type: mdl:ols
#>              Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.148462   0.014677  -10.12   <2e-16 ***
#> (Intercept) -0.000578   0.006882   -0.08     0.93    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Ensemble type: mdl:lasso
#>              Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.148749   0.014680   -10.1   <2e-16 ***
#> (Intercept) -0.000689   0.006884    -0.1     0.92    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Ensemble type: mdl:ridge
#>              Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.148367   0.014675  -10.11   <2e-16 ***
#> (Intercept) -0.000737   0.006884   -0.11     0.91    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# \donttest{
# Re-estimate with a different ensemble type using pass-through
#     (skips cross-fitting, only recomputes ensemble weights).
plm_fit2 <- ddml_plm(y, D, X,
                     learners = list(list(what = ols),
                                     list(what = mdl_glmnet),
                                     list(what = mdl_glmnet,
                                          args = list(alpha = 0))),
                     ensemble_type = 'average',
                     shortstack = TRUE,
                     sample_folds = 2,
                     silent = TRUE,
                     fitted = plm_fit$fitted,
                     splits = plm_fit$splits)
summary(plm_fit2)
#> DDML estimation: Partially Linear Model 
#> Obs: 5000   Folds: 2  Stacking: short-stack
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.148554   0.014677   -10.1   <2e-16 ***
#> (Intercept) -0.000668   0.006883    -0.1     0.92    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
