# Estimator for Group-Time Average Treatment Effects

Estimator for group-time average treatment effects on the treated
(GT-ATT) in staggered Difference-in-Differences designs.

## Usage

``` r
ddml_attgt(
  y,
  X = NULL,
  t,
  G,
  learners,
  learners_qX = learners,
  sample_folds = 10,
  ensemble_type = "nnls",
  shortstack = FALSE,
  cv_folds = 10,
  custom_ensemble_weights = NULL,
  custom_ensemble_weights_qX = custom_ensemble_weights,
  cluster_variable = seq_len(nrow(as.matrix(y))),
  trim = 0.01,
  control_group = c("notyettreated", "nevertreated"),
  anticipation = 0,
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

  An \\n \times T\\ numeric matrix of outcomes. Row \\i\\ corresponds to
  unit \\i\\, column \\j\\ to time period `t[j]`.

- X:

  An \\n \times p\\ matrix of time-invariant covariates, or `NULL`.

- t:

  A numeric vector of length \\T\\ giving the time period labels (must
  match columns of `y`).

- G:

  A numeric vector of length \\n\\. Entry \\i\\ is the first treatment
  period for unit \\i\\. Use `0` or `Inf` for never-treated units.

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

- learners_qX:

  Optional argument to allow for different estimators of the cell-level
  propensity score \\q^{(g,t)}(X)\\. Setup is identical to `learners`.

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

- custom_ensemble_weights_qX:

  Optional argument to allow for different custom ensemble weights for
  `learners_qX`. Setup is identical to `custom_ensemble_weights`.

- cluster_variable:

  A vector of cluster indices.

- trim:

  Number in (0, 1) for trimming the estimated propensity scores at
  `trim` and `1-trim`.

- control_group:

  Character. `"notyettreated"` (default) uses never-treated and
  not-yet-treated units as controls. `"nevertreated"` uses only
  never-treated units.

- anticipation:

  Non-negative integer. Number of periods before treatment where
  anticipation effects may occur. Default 0.

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

`ddml_attgt` returns an object of S3 class `ddml_attgt` and `ddml`. See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for the common output structure. Additional pass-through fields:
`learners`, `learners_qX`, `cell_info`, `control_group`, `anticipation`.

## Details

**Parameter of Interest:** `ddml_attgt` provides a Double/Debiased
Machine Learning estimator for the group-time average treatment effects
on the treated (GT-ATT) in the staggered adoption model. For each group
\\g\\ and time period \\t\\, define the differenced outcome \\\Delta_g
Y\_{i,t} = Y\_{i,t} - Y\_{i,g^\*}\\ where \\g^\*\\ is the universal base
period. The GT-ATT is:

\$\$\theta_0^{(g,t)} = E\[\Delta_g Y\_{i,t} \| G_i = g\] -
E\[E\[\Delta_g Y\_{i,t} \| X_i, G_i \ne g, G_i \> t\] \| G_i = g\]\$\$

where \\W_i \equiv (Y\_{i,1}, \dots, Y\_{i,T}, G_i, X_i)\\ is the
observed random vector.

**Neyman Orthogonal Score:** The Neyman orthogonal score is:

\$\$m^{(g,t)}(W_i; \theta, \eta) = \frac{\mathbf{1}\\G_i = g\\ (\Delta_g
Y\_{i,t} - \ell^{(g,t)}(X_i))}{\pi^g} - \frac{q^{(g,t)}(X_i)
\mathbf{1}\\G_i \ne g\\ \mathbf{1}\\G_i \> t\\ (\Delta_g Y\_{i,t} -
\ell^{(g,t)}(X_i))}{\pi^g (1 - q^{(g,t)}(X_i))} - \frac{\mathbf{1}\\G_i
= g\\}{\pi^g} \theta\$\$

where the nuisance parameters are \\\eta = (\ell, q, \pi)\\ taking true
values \\\ell_0^{(g,t)}(X) = E\[\Delta_g Y\_{i,t} \mid G_i \ne g, G_i \>
t, X_i\]\\, \\q_0^{(g,t)}(X) = \Pr(G_i = g \mid X_i, \\G_i = g\\ \cup
\\G_i \> t\\)\\, and \\\pi_0^g = \Pr(G_i = g)\\.

**Jacobian:**

\$\$J^{(g,t)} = -1\$\$

See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for how the influence function and inference are derived from these
components.

## References

Callaway B, Sant'Anna P H C (2021). "Difference-in-Differences with
multiple time periods." Journal of Econometrics, 225(2), 200-230.

Chang N-C (2020). "Double/debiased machine learning for
difference-in-differences models." Econometrics Journal, 23(2), 177-191.

Ahrens A, Chernozhukov V, Hansen C B, Kozbur D, Schaffer M E, Wiemann T
(2026). "An Introduction to Double/Debiased Machine Learning." Journal
of Economic Literature, forthcoming.

## See also

Other ddml estimators:
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md),
[`ddml_apo()`](https://www.thomaswiemann.com/ddml/reference/ddml_apo.md),
[`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md),
[`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md),
[`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md),
[`ddml_policy()`](https://www.thomaswiemann.com/ddml/reference/ddml_policy.md)

## Examples

``` r
# \donttest{
set.seed(42)
n <- 200; T_ <- 4
X <- matrix(rnorm(n * 2), n, 2)
G <- sample(c(3, 4, Inf), n, replace = TRUE,
            prob = c(0.3, 0.3, 0.4))
y <- matrix(rnorm(n * T_), n, T_)
# Add treatment effect for treated units
for (i in seq_len(n)) {
  if (is.finite(G[i])) {
    for (j in seq_len(T_)) {
      if (j >= G[i]) y[i, j] <- y[i, j] + 1
    }
  }
}
fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                learners = list(what = ols),
                sample_folds = 2,
                silent = TRUE)
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
summary(fit)
#> DDML estimation: Group-Time Average Treatment Effects on the Treated 
#> Obs: 200   Folds: 2
#> 
#>          Estimate Std. Error z value Pr(>|z|)    
#> ATT(3,1)  -0.1127     0.2426   -0.46   0.6423    
#> ATT(3,3)   1.0521     0.2422    4.34  1.4e-05 ***
#> ATT(3,4)   1.0890     0.2678    4.07  4.8e-05 ***
#> ATT(4,1)  -0.0144     0.3107   -0.05   0.9629    
#> ATT(4,2)  -0.2359     0.2786   -0.85   0.3972    
#> ATT(4,4)   1.0364     0.3022    3.43   0.0006 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
