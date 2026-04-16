# Predict Method for mdl_xgboost Objects

Predict Method for mdl_xgboost Objects

## Usage

``` r
# S3 method for class 'mdl_xgboost'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted `mdl_xgboost` object.

- newdata:

  A feature matrix for prediction.

- ...:

  Additional arguments passed to
  [`predict.xgb.Booster`](https://rdrr.io/pkg/xgboost/man/predict.xgb.Booster.html).

## Value

A numeric vector of predicted values.
