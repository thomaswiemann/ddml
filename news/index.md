# Changelog

## ddml 0.9.0

- Adds
  [`ddml_attgt()`](https://www.thomaswiemann.com/ddml/reference/ddml_attgt.md)
  for staggered DiD and
  [`ddml_apo()`](https://www.thomaswiemann.com/ddml/reference/ddml_apo.md)
  for average potential outcomes.
- Adds
  [`ddml_policy()`](https://www.thomaswiemann.com/ddml/reference/ddml_policy.md)
  for multi-action policy value estimation.
- Adds [`ddml()`](https://www.thomaswiemann.com/ddml/reference/ddml.md)
  constructor for custom DML estimators with user-supplied scores.
- Adds
  [`lincom()`](https://www.thomaswiemann.com/ddml/reference/lincom.md)
  for inference on linear combinations. Supports computation of dynamic
  average treatment effects via
  [`lincom_weights_did()`](https://www.thomaswiemann.com/ddml/reference/lincom_weights_did.md).
- Influence-function-based inference via the `ral` class; all estimators
  now inherit from `ral`.
- Adds
  [`ddml_rep()`](https://www.thomaswiemann.com/ddml/reference/ddml_rep.md)
  and
  [`ddml_replicate()`](https://www.thomaswiemann.com/ddml/reference/ddml_replicate.md)
  for repeated cross-fitting with median, mean, or spectral-norm
  aggregation.
- Adds
  [`diagnostics()`](https://www.thomaswiemann.com/ddml/reference/diagnostics.md)
  for MSPE, R-squared, stacking weights, and CVC tests.
- Adds `fitted`/`splits` pass-through to all `ddml_*()` estimators.
- New S3 methods:
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`as.list()`](https://rdrr.io/r/base/list.html),
  [`hatvalues()`](https://rdrr.io/r/stats/influence.measures.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html), multi-ensemble
  [`tidy()`](https://generics.r-lib.org/reference/tidy.html)/[`glance()`](https://generics.r-lib.org/reference/glance.html).
- Adds uniform confidence bands via multiplier bootstrap
  (`confint(uniform = TRUE)`).
- Adds HC0/HC3 variance estimators, parallel computation, stratified
  cross-fitting, cluster-aware splitting, and input validation.
- Adds `broom` compatibility.
- Fixes
  [`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md)
  with custom weights.
- Breaking changes:
  - Inference internals use `$inf_func` instead of
    `$scores`/`$J`/`$psi_a`/`$psi_b`.
  - Utility functions (`crosspred`, `crossval`, `ensemble`,
    `ensemble_weights`, `shortstacking`) no longer accept `Z`/`newZ`.
    Pre-concatenate instruments with covariates (e.g., `cbind(X, Z)`).
  - [`crosspred()`](https://www.thomaswiemann.com/ddml/reference/crosspred.md)
    and
    [`shortstacking()`](https://www.thomaswiemann.com/ddml/reference/shortstacking.md)
    drop `compute_insample_predictions` and `insample_fitted` output.
  - [`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md)
    drops the `enforce_LIE` argument.
  - [`shortstacking()`](https://www.thomaswiemann.com/ddml/reference/shortstacking.md)
    drops `shortstack_y`.
  - `ddml_*()` estimators drop `subsamples`, `cv_subsamples`,
    `subsamples_byD`, `cv_subsamples_byD`. Use the new `splits`
    parameter instead.

## ddml 0.3.0

CRAN release: 2024-10-02

- Implements one-way clustered inference.
- Increases defaults for `sample_folds` and `cv_folds` to `10`.
- Fixes typo in `auxiliary_X` arguments.

## ddml 0.2.2

CRAN release: 2024-06-26

- Changes
  [`ddml::ols()`](https://www.thomaswiemann.com/ddml/reference/ols.md)
  default to `const=TRUE`.
- Adds probability forest compatibility to
  [`ddml::mdl_ranger()`](https://www.thomaswiemann.com/ddml/reference/mdl_ranger.md).
- Adds propensity score trimming option to
  [`ddml::ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
  [`ddml::ddml_att()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
  and
  [`ddml::ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md).
- Fixes ATE and LATE scores.
- Fixes output of `ddml::print.summary.ddml_plm` and
  `ddml::print.summary.ddml_ate`
  ([\#57](https://github.com/thomaswiemann/ddml/issues/57)).

## ddml 0.2.1

CRAN release: 2024-05-26

- Fixes permuted residuals returned by
  [`ddml::crossval`](https://www.thomaswiemann.com/ddml/reference/crossval.md)
  ([\#54](https://github.com/thomaswiemann/ddml/issues/54)).

## ddml 0.2.0

CRAN release: 2024-01-09

- Adds support for the average treatment effect on the treated
  estimator.
- Adds support for local average treatment effect estimation with
  perfect compliance or perfect non-compliance.
- Adds support for custom ensemble weights.
- Adds article on integration with the
  [`did`](https://bcallaway11.github.io/did/) package.
- Adds
  [`ddml::mdl_glm`](https://www.thomaswiemann.com/ddml/reference/mdl_glm.md)
  wrapper for [`stats::glm()`](https://rdrr.io/r/stats/glm.html).

## ddml 0.1.0

CRAN release: 2023-08-29

- Initial CRAN submission.
