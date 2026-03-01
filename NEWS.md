# ddml 0.3.1

* Adds `fitted` and `splits` parameters to all `ddml_*()` estimators for re-estimation with different ensemble types without re-fitting base learners.
* Adds `diagnostics()` for per-equation MSPE, R-squared, stacking weights, and cross-validated comparison (CVC) tests.
* Adds HC0 and HC3 variance estimators alongside default HC1 via `type` argument in `vcov()`, `confint()`, `summary()`, and `tidy()`.
* Learner specifications now use `what` instead of `fun`. `fun` is accepted with a deprecation message.
* Sample split arguments consolidated into a single `splits` parameter.
* Adds parallel computation via the `parallel` parameter and progress reporting via `pbapply`.
* Adds stratified cross-fitting to `ddml_ate()`, `ddml_att()`, and `ddml_late()` via `stratify`. Enabled by default.
* Adds `cluster_variable` for cluster-aware sample splitting.
* Adds input validation across all estimators.
* Exports `ensemble()` with documentation and examples.
* Adds `broom` compatibility (`tidy`, `glance`).
* Updates S3 methods (`coef`, `vcov`, `confint`, `summary`, `print`).
* Adds `mdl_bigGLM` sparse-matrix unpenalized regression wrapper.
* `mdl_xgboost()` auto-detects binary factor outcomes.
* Fixes `ddml_fpliv()` with custom weights.

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
