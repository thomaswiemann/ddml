# Package index

## Supported DDML Estimators

- [`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md)
  : Estimator for the Partially Linear Regression Coefficient
- [`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md)
  : Estimator for the Partially Linear IV Coefficient
- [`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md)
  : Estimator for the Flexible Partially Linear IV Coefficient
- [`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md)
  [`ddml_att()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md)
  : Estimator for the Average Treatment Effect
- [`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md)
  : Estimator for the Local Average Treatment Effect
- [`ddml_apo()`](https://www.thomaswiemann.com/ddml/reference/ddml_apo.md)
  : Estimator for the Average Potential Outcome
- [`ddml_policy()`](https://www.thomaswiemann.com/ddml/reference/ddml_policy.md)
  : Estimator for the Multi-Action Policy Value
- [`ddml_attgt()`](https://www.thomaswiemann.com/ddml/reference/ddml_attgt.md)
  : Estimator for Group-Time Average Treatment Effects

## Repeated Resampling

- [`ddml_replicate()`](https://www.thomaswiemann.com/ddml/reference/ddml_replicate.md)
  : Replicate a DDML Estimator Across Multiple Resamples

## Linear Combinations and Aggregation

- [`lincom()`](https://www.thomaswiemann.com/ddml/reference/lincom.md)
  [`print(`*`<lincom>`*`)`](https://www.thomaswiemann.com/ddml/reference/lincom.md)
  [`print(`*`<lincom_rep>`*`)`](https://www.thomaswiemann.com/ddml/reference/lincom.md)
  : Linear Combinations of DDML Coefficients
- [`lincom_weights_did()`](https://www.thomaswiemann.com/ddml/reference/lincom_weights_did.md)
  : Difference-in-Differences Aggregation Weights for lincom

## Wrappers for Common (Machine) Learners

- [`ols()`](https://www.thomaswiemann.com/ddml/reference/ols.md) :
  Ordinary Least Squares
- [`mdl_glm()`](https://www.thomaswiemann.com/ddml/reference/mdl_glm.md)
  : Wrapper for stats::glm()
- [`mdl_glmnet()`](https://www.thomaswiemann.com/ddml/reference/mdl_glmnet.md)
  : Wrapper for glmnet::glmnet()
- [`mdl_ranger()`](https://www.thomaswiemann.com/ddml/reference/mdl_ranger.md)
  : Wrapper for ranger::ranger()
- [`mdl_xgboost()`](https://www.thomaswiemann.com/ddml/reference/mdl_xgboost.md)
  : Wrapper for xgboost::xgboost()
- [`mdl_bigGlm()`](https://www.thomaswiemann.com/ddml/reference/mdl_bigGLM.md)
  : Wrapper for glmnet::bigGlm()

## Utilities

- [`crossval()`](https://www.thomaswiemann.com/ddml/reference/crossval.md)
  : Estimator of the Mean Squared Prediction Error Using
  Cross-Validation
- [`crosspred()`](https://www.thomaswiemann.com/ddml/reference/crosspred.md)
  : Cross-Fitted Predictions Using Stacking
- [`shortstacking()`](https://www.thomaswiemann.com/ddml/reference/shortstacking.md)
  : Predictions using Short-Stacking
- [`ensemble()`](https://www.thomaswiemann.com/ddml/reference/ensemble.md)
  : Stacking Estimator Using Combinations of Base Learners
- [`diagnostics()`](https://www.thomaswiemann.com/ddml/reference/diagnostics.md)
  : Stacking Diagnostics for DDML Estimators

## Dataset

- [`AE98`](https://www.thomaswiemann.com/ddml/reference/AE98.md) :
  Random Subsample from the Data of Angrist & Evans (1998)
