# Changelog

## ddml 0.3.1

CRAN release: 2025-12-11

- Updates internals of
  [`ddml::mdl_xgboost()`](https://www.thomaswiemann.com/ddml/reference/mdl_xgboost.md)
  with new `xgboost` syntax.
- Fixes
  [`ddml::ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md)
  with custom weights.
- Allows for stacking with no positive stacking weights.
- Fixes
  [`ddml::mdl_glmnet`](https://www.thomaswiemann.com/ddml/reference/mdl_glmnet.md)
  predictions for binomial regression.

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
- Fixes output of
  [`ddml::print.summary.ddml_plm`](https://www.thomaswiemann.com/ddml/reference/print.summary.ddml_plm.md)
  and
  [`ddml::print.summary.ddml_ate`](https://www.thomaswiemann.com/ddml/reference/print.summary.ddml_ate.md)
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
