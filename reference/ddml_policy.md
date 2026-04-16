# Estimator for the Multi-Action Policy Value

Estimator for the expected value of a multi-action policy, with optional
per-level margins.

## Usage

``` r
ddml_policy(
  y,
  D,
  X,
  policy,
  margins = NULL,
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

  The observed discrete (potentially multi-valued) treatment variable.

- X:

  A (sparse) matrix of control variables.

- policy:

  A vector of length `nobs` giving the policy-assigned treatment level
  for each unit. Values must be a subset of those observed in `D`.

- margins:

  An optional numeric vector of length \\K\\ (the number of unique
  values in `policy`) giving per-level multipliers \\c_k\\. If `NULL`
  (the default), all margins are set to one, yielding the policy value
  \\E\[Y(\pi(X))\]\\.

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

  An optional list of sample split objects. For `ddml_policy`, this is
  typically a named list keyed by treatment level. Typically obtained
  from a previous fit via `fit$splits`.

- save_crossval:

  Logical indicating whether to store the inner cross-validation
  residuals used for ensemble weight computation. Default `TRUE`. When
  `TRUE`, subsequent pass-through calls with data-driven ensembles
  (e.g., `"nnls"`) reproduce per-fold weights exactly. Set to `FALSE` to
  reduce object size at the cost of approximate weight recomputation.

- ...:

  Additional arguments passed to internal methods.

## Value

`ddml_policy` returns an object of S3 class `ddml_policy` and `ddml`.
See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for the common output structure. Additional pass-through fields:
`learners`, `learners_DX`, `policy`, `margins`.

## Details

**Parameter of Interest:** `ddml_policy` provides a Double/Debiased
Machine Learning estimator for the expected value of a multi-action
policy \\\pi:\operatorname{supp}(X) \to \\d_1, \ldots, d_K\\\\ that
assigns each unit to one of \\K\\ treatment levels. The target parameter
is

\$\$\theta_0 = \sum\_{k=1}^{K} E\\\left\[\omega_k(X)\\ E\[Y \mid D =
d_k, X\]\right\],\$\$

where the known weight functions are \\\omega_k(X) = c_k
\\\mathbf{1}\\\pi(X) = d_k\\\\, and \\c_1, \ldots, c_K\\ are
user-supplied margins. When all margins equal one (`margins = NULL`),
the parameter reduces to the policy value \\E\[Y(\pi(X))\]\\.

Each term in the sum is a weighted average potential outcome (wAPO),
estimated internally via
[`ddml_apo`](https://www.thomaswiemann.com/ddml/reference/ddml_apo.md).

**Nuisance Parameters:** For each treatment level \\d_k\\, the nuisance
parameters are \\\eta_k = (g_k, m_k)\\ taking true values \\g\_{k,0}(X)
= E\[Y \mid D = d_k, X\]\\ and \\m\_{k,0}(X) = \Pr(D = d_k \mid X)\\.
Only \\K-1\\ propensity models are estimated; the last is derived as
\\m_K(X) = 1 - \sum\_{k=1}^{K-1} m_k(X)\\.

**Neyman Orthogonal Score / Moment Equation:** The Neyman orthogonal
score is:

\$\$m(W; \theta, \eta) = \sum\_{k=1}^{K} \omega_k(X) \left\[
\frac{\mathbf{1}\\D = d_k\\\\(Y - g_k(X))}{m_k(X)} + g_k(X) \right\] -
\theta\$\$

**Jacobian:**

\$\$J = -1\$\$

See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for how the influence function and inference are derived from these
components.

## References

Dudik M, Langford J, Li L (2011). "Doubly Robust Policy Evaluation and
Learning." Proceedings of the 28th International Conference on Machine
Learning, 1097-1104.

Zhou Z, Athey S, Wager S (2023). "Offline Multi-Action Policy Learning:
Generalization and Optimization." Operations Research, 71(2), 698-722.

## See also

Other ddml estimators:
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md),
[`ddml_apo()`](https://www.thomaswiemann.com/ddml/reference/ddml_apo.md),
[`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_attgt()`](https://www.thomaswiemann.com/ddml/reference/ddml_attgt.md),
[`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md),
[`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md),
[`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md)

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]

# Define a simple policy: assign D=1 if age > median, else D=0
policy <- ifelse(X[, "age"] > median(X[, "age"]), 1, 0)

# Estimate the policy value using a single base learner, ridge.
policy_fit <- ddml_policy(y, D, X,
                          policy = policy,
                          learners = list(what = mdl_glmnet),
                          sample_folds = 2,
                          silent = TRUE)
summary(policy_fit)
#> DDML estimation: Multi-Action Policy Value 
#> Obs: 5000   Folds: 2
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> Policy value   0.5257     0.0102    51.3   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
