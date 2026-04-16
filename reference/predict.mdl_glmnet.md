# Predict Method for mdl_glmnet Objects

Predict Method for mdl_glmnet Objects

## Usage

``` r
# S3 method for class 'mdl_glmnet'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted `mdl_glmnet` object.

- newdata:

  A (sparse) feature matrix for prediction.

- ...:

  Additional arguments passed to
  [`predict.glmnet`](https://glmnet.stanford.edu/reference/predict.glmnet.html).

## Value

A numeric vector of predicted values.
