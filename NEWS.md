# ddml 0.3.1

* Adds support for parallel computation in `crosspred()`, `crossval()`, and `shortstacking()` via the `parallel` parameter. Updates progress reporting using `pbapply`.
* Adds stratified cross-fitting to ``ddml_ate()``, ``ddml_att()``, and
  ``ddml_late()`` via the ``stratify`` parameter.
  Stratified splitting is now the default.
* Adds ``cluster_variable`` to ``crossval()``, ``crosspred()``, and
  ``shortstacking()`` for cluster-aware sample splitting.
* Renames ``cv_subsamples_list`` to ``cv_subsamples``. The old name
  still works with a deprecation message.
* Replaces internal ``get_all_indx()`` with ``get_sample_splits()``.
* Warns when any cross-fitting training set has fewer than 100
  observations.
* Updates internals of ``ddml::mdl_xgboost()`` with new ``xgboost`` syntax.
* Fixes ``ddml::ddml_fpliv()`` with custom weights.
* Allows for stacking with no positive stacking weights.
* Fixes ``ddml::mdl_glmnet`` predictions for binomial regression.

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
