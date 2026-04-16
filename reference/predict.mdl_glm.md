# Predict Method for mdl_glm Objects

Predict Method for mdl_glm Objects

## Usage

``` r
# S3 method for class 'mdl_glm'
predict(object, newdata, ...)
```

## Arguments

- object:

  A fitted `mdl_glm` object.

- newdata:

  A feature matrix for prediction.

- ...:

  Additional arguments passed to
  [`predict.glm`](https://rdrr.io/r/stats/predict.glm.html).

## Value

A numeric vector of predicted response values.
