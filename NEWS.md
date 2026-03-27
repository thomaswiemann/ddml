# ddml 0.9.0

* Adds `ddml_attgt()` for staggered DiD and `ddml_apo()` for average potential outcomes.
* Adds `ddml()` constructor for custom DML estimators with user-supplied scores.
* Adds `lincom()` and `lincom_weights_did()` for inference on linear combinations.
* Influence-function-based inference via the `ral` class; all estimators now inherit from `ral`.
* Adds `ddml_rep()` and `ddml_replicate()` for repeated cross-fitting with median, mean, or spectral-norm aggregation.
* Adds `diagnostics()` for MSPE, R-squared, stacking weights, and CVC tests.
* Adds `fitted`/`splits` pass-through to all `ddml_*()` estimators.
* New S3 methods: `plot()`, `as.list()`, `hatvalues()`, `nobs()`, multi-ensemble `tidy()`/`glance()`.
* Adds uniform confidence bands via multiplier bootstrap (`confint(uniform = TRUE)`).
* Adds HC0/HC3 variance estimators, parallel computation, stratified cross-fitting, cluster-aware splitting, and input validation.
* Adds `broom` compatibility.
* Fixes `ddml_fpliv()` with custom weights.
* Breaking changes:
    - Inference internals use `$inf_func` instead of `$scores`/`$J`/`$psi_a`/`$psi_b`.
    - Utility functions (`crosspred`, `crossval`, `ensemble`, `ensemble_weights`, `shortstacking`) no longer accept `Z`/`newZ`. Pre-concatenate instruments with covariates (e.g., `cbind(X, Z)`).
    - `crosspred()` and `shortstacking()` drop `compute_insample_predictions` and `insample_fitted` output.
    - `ddml_fpliv()` drops the `enforce_LIE` argument.
    - `shortstacking()` drops `shortstack_y`.
    - `ddml_*()` estimators drop `subsamples`, `cv_subsamples`, `subsamples_byD`, `cv_subsamples_byD`. Use the new `splits` parameter instead.

# ddml 0.3.0

* Implements one-way clustered inference.
* Increases defaults for ``sample_folds`` and ``cv_folds`` to ``10``.
* Fixes typo in ``auxiliary_X`` arguments.

# ddml 0.2.2

* Changes ``ddml::ols()`` default to ``const=TRUE``.
* Adds probability forest compatibility to ``ddml::mdl_ranger()``.
* Adds propensity score trimming option to ``ddml::ddml_ate()``, ``ddml::ddml_att()``, and ``ddml::ddml_late()``.
* Fixes ATE and LATE scores.
* Fixes output of ``ddml::print.summary.ddml_plm`` and ``ddml::print.summary.ddml_ate`` (#57).

# ddml 0.2.1

* Fixes permuted residuals returned by ``ddml::crossval`` (#54).

# ddml 0.2.0

* Adds support for the average treatment effect on the treated estimator.
* Adds support for local average treatment effect estimation with perfect compliance or perfect non-compliance.
* Adds support for custom ensemble weights.
* Adds article on integration with the [``did``](https://bcallaway11.github.io/did/) package.
* Adds ``ddml::mdl_glm`` wrapper for ``stats::glm()``.

# ddml 0.1.0

* Initial CRAN submission.
