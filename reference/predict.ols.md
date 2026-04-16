# Predict Method for ols Objects

Predict Method for ols Objects

## Usage

``` r
# S3 method for class 'ols'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted `ols` object.

- newdata:

  A feature matrix for prediction. If `NULL`, returns fitted values from
  the training data.

- ...:

  Currently unused.

## Value

A numeric vector of predicted values.
