# Intro to Double/Debiased Machine Learning

All `ddml_*` estimators
([`ddml_plm`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md),
[`ddml_pliv`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md),
[`ddml_fpliv`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_ate`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_att`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_late`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md),
[`ddml_apo`](https://www.thomaswiemann.com/ddml/reference/ddml_apo.md),
[`ddml_policy`](https://www.thomaswiemann.com/ddml/reference/ddml_policy.md))
return objects that inherit from S3 class `"ddml"`.

Each object is a list containing the components described below.
Estimator-specific fields (e.g., pass-through learner arguments) are
documented on the individual estimator pages.

The [`ddml()`](https://www.thomaswiemann.com/ddml/reference/ddml.md)
constructor can also be used directly to build a `"ddml"` object from
user-supplied score components, enabling implementation of custom DML
estimators that inherit all S3 methods.

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

## Details

All `ddml_*` estimators target a low-dimensional parameter \\\theta_0\\
identified by a moment condition

\$\$E\[m(W; \theta_0, \eta_0)\] = 0,\$\$

where \\W\\ denotes observed random variables and \\\eta_0\\ is a
(potentially high-dimensional) nuisance parameter. Throughout, the score
\\m\\ is assumed to be *Neyman orthogonal*.

Estimation proceeds via cross-fitting: the sample is randomly
partitioned into \\K\\ folds \\\\I_k\\\_{k=1}^K\\. For each fold \\k\\,
nuisance parameters are estimated on the complementary folds
(\\\hat\eta\_{-k}\\) and the scores are evaluated on fold \\k\\. The DML
estimator \\\hat\theta\\ solves

\$\$\frac{1}{n} \sum\_{k=1}^{K} \sum\_{i \in I_k} m(W_i; \hat\theta,
\hat\eta\_{-k}) = 0.\$\$

Inference is based on the influence function. Define the Jacobian

\$\$J(\theta, \eta) = E\\\left\[ \frac{\partial m(W; \theta, \eta)}
{\partial \theta'}\right\]\$\$

and the influence function

\$\$\phi\_\theta(W_i; \theta, \eta, J) = -J^{-1}\\m(W_i; \theta,
\eta).\$\$

The variance of \\\hat\theta\\ is then estimated by

\$\$\hat{V} = \frac{1}{n^2} \sum_i \phi\_\theta(W_i; \hat\theta,
\hat\eta\_{-k(i)}, \hat{J})\\\phi\_\theta(W_i; \hat\theta,
\hat\eta\_{-k(i)}, \hat{J})'\$\$,

where \\\hat{J}\\ is the sample analog of the Jacobian:

\$\$\hat{J} = \frac{1}{n} \sum_i \frac{\partial m(W_i; \hat\theta,
\hat\eta\_{-k(i)})} {\partial \theta'}.\$\$

HC1 and HC3 variance estimators are described in
[`vcov.ral`](https://www.thomaswiemann.com/ddml/reference/vcov.ral.md).
The leverage (see
[`hatvalues.ral`](https://www.thomaswiemann.com/ddml/reference/hatvalues.ral.md))
for the DML estimator is

\$\$h\_\theta(W_i; \theta, \eta, J) = \mathrm{tr}\\\left( -J^{-1}
\frac{1}{n} \frac{\partial m(W_i; \theta, \eta)} {\partial
\theta'}\right),\$\$

and its sample analog is \\\hat{h}\_{\theta,i} = h\_\theta(W_i;
\hat\theta, \hat\eta\_{-k(i)}, \hat{J})\\, stored in `dinf_dtheta`.

Under regularity conditions and sufficient convergence of \\\hat\eta\\,
the DML estimator is asymptotically normal:

\$\$\sqrt{n}\\\hat{V}^{-1/2}(\hat\theta - \theta_0) \overset{d}{\to}
N(0, I).\$\$

Further details and regularity conditions are given in Chernozhukov et
al. (2018). The specific forms of the score \\m\\ and Jacobian \\J\\ for
each estimator are documented on their respective help pages (e.g.,
[`ddml_plm`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md),
[`ddml_ate`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md)).

## Common output components

- `coefficients`:

  A matrix of estimated target parameters: rows correspond to components
  of \\\theta\\, columns to ensemble types.

- `ensemble_weights`:

  A named list. Each element is a weight matrix (or 3D array when
  `shortstack = TRUE`) showing the weight assigned to each base learner
  by the ensemble procedure for the corresponding nuisance equation.

- `mspe`:

  A named list of numeric vectors containing per-learner out-of-sample
  MSPEs, computed from cross-fitted residuals.

- `r2`:

  A named list of numeric vectors containing per-learner out-of-sample
  R-squared values.

- `inf_func`:

  A 3D array of evaluated influence functions (`n x p x nensb`).

- `dinf_dtheta`:

  An optional 4D array of dimension `(n x p x p x nensb)` containing the
  derivatives of the influence functions with respect to \\\theta\\.
  Used internally by
  [`hatvalues.ral`](https://www.thomaswiemann.com/ddml/reference/hatvalues.ral.md)
  for HC3 inference.

- `scores`:

  A 3D array of evaluated Neyman orthogonal scores (`n x p x nensb`).

- `J`:

  A 3D array of evaluated Jacobians (`p x p x nensb`).

- `fitted`:

  A named list of per-equation cross-fitted prediction objects. Can be
  passed back via the `fitted` argument together with `splits` to skip
  cross-fitting on re-estimation.

- `splits`:

  The data splitting structure (subsamples, CV subsamples, and any
  stratification indices).

- `ensemble_type`:

  Character vector of ensemble types used.

- `cluster_variable`:

  The cluster variable vector used for sample splitting and inference.

- `nobs`:

  Number of observations.

- `sample_folds`:

  Number of cross-fitting folds.

- `shortstack`:

  Logical indicating whether short-stacking was used.

- `call`:

  The matched call.

- `coef_names`:

  Character vector of coefficient names.

- `estimator_name`:

  Character string identifying the estimator (e.g.,
  `"Partially Linear Model"`).

## S3 methods

The following generic methods are available for all `ddml` objects:
[`summary.ddml`](https://www.thomaswiemann.com/ddml/reference/summary.ddml.md),
[`coef.ral`](https://www.thomaswiemann.com/ddml/reference/coef.ral.md),
[`vcov.ral`](https://www.thomaswiemann.com/ddml/reference/vcov.ral.md),
[`confint.ral`](https://www.thomaswiemann.com/ddml/reference/confint.ral.md),
[`hatvalues.ral`](https://www.thomaswiemann.com/ddml/reference/hatvalues.ral.md),
[`nobs.ral`](https://www.thomaswiemann.com/ddml/reference/nobs.ral.md),
[`tidy.ddml`](https://www.thomaswiemann.com/ddml/reference/tidy.ddml.md),
[`glance.ddml`](https://www.thomaswiemann.com/ddml/reference/glance.ddml.md),
and
[`diagnostics`](https://www.thomaswiemann.com/ddml/reference/diagnostics.md).

## References

Ahrens A, Chernozhukov V, Hansen C B, Kozbur D, Schaffer M E, Wiemann T
(2026). "An Introduction to Double/Debiased Machine Learning." Journal
of Economic Literature, forthcoming.

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). "Double/debiased machine learning for treatment and
structural parameters." The Econometrics Journal, 21(1), C1-C68.

## See also

Other ddml estimators:
[`ddml_apo()`](https://www.thomaswiemann.com/ddml/reference/ddml_apo.md),
[`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_attgt()`](https://www.thomaswiemann.com/ddml/reference/ddml_attgt.md),
[`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md),
[`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md),
[`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md),
[`ddml_policy()`](https://www.thomaswiemann.com/ddml/reference/ddml_policy.md)
