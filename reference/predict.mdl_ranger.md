# Predict Method for mdl_ranger Objects

Predict Method for mdl_ranger Objects

## Usage

``` r
# S3 method for class 'mdl_ranger'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted `mdl_ranger` object.

- newdata:

  A feature matrix for prediction.

- ...:

  Additional arguments passed to
  [`predict.ranger`](http://imbs-hl.github.io/ranger/reference/predict.ranger.md).

## Value

A numeric vector of predicted values (probabilities for probability
forests, point predictions for regression forests).
