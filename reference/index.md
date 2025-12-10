# Package index

## Supported Models

- [`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md)
  : Estimator for the Partially Linear Model.
- [`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md)
  [`ddml_att()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md)
  : Estimators of Average Treatment Effects.
- [`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md)
  : Estimator for the Partially Linear IV Model.
- [`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md)
  : Estimator for the Flexible Partially Linear IV Model.
- [`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md)
  : Estimator of the Local Average Treatment Effect.

## Wrappers for Common (Machine) Learners

- [`ols()`](https://www.thomaswiemann.com/ddml/reference/ols.md) :
  Ordinary least squares.

- [`mdl_glm()`](https://www.thomaswiemann.com/ddml/reference/mdl_glm.md)
  :

  Wrapper for [`stats::glm()`](https://rdrr.io/r/stats/glm.html).

- [`mdl_glmnet()`](https://www.thomaswiemann.com/ddml/reference/mdl_glmnet.md)
  :

  Wrapper for
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).

- [`mdl_ranger()`](https://www.thomaswiemann.com/ddml/reference/mdl_ranger.md)
  :

  Wrapper for
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).

- [`mdl_xgboost()`](https://www.thomaswiemann.com/ddml/reference/mdl_xgboost.md)
  :

  Wrapper for
  [`xgboost::xgboost()`](https://rdrr.io/pkg/xgboost/man/xgboost.html).

## Utilities

- [`crossval()`](https://www.thomaswiemann.com/ddml/reference/crossval.md)
  : Estimator of the Mean Squared Prediction Error using
  Cross-Validation.
- [`crosspred()`](https://www.thomaswiemann.com/ddml/reference/crosspred.md)
  : Cross-Predictions using Stacking.
- [`shortstacking()`](https://www.thomaswiemann.com/ddml/reference/shortstacking.md)
  : Predictions using Short-Stacking.

## Dataset

- [`AE98`](https://www.thomaswiemann.com/ddml/reference/AE98.md) :
  Random subsample from the data of Angrist & Evans (1991).
