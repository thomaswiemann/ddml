# Predict Method for `ensemble` Objects

Predict Method for `ensemble` Objects

## Usage

``` r
# S3 method for class 'ensemble'
predict(object, newdata, ..., type = "ensemble")
```

## Arguments

- object:

  A fitted `ensemble` object.

- newdata:

  A feature matrix for prediction.

- ...:

  Currently unused.

- type:

  Character; `"ensemble"` (default) returns weighted ensemble
  predictions, `"bylearner"` returns the raw per-learner prediction
  matrix.

## Value

A matrix of predictions. When `type = "ensemble"`, one column per
ensemble type; when `type = "bylearner"`, one column per base learner.
