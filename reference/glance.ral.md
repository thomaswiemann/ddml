# Glance at a RAL Object

Returns a one-row summary of model-level statistics.

## Usage

``` r
# S3 method for class 'ral'
glance(x, ...)
```

## Arguments

- x:

  An object inheriting from class `ral`.

- ...:

  Currently unused.

## Value

A one-row `data.frame` with columns `nobs` and `estimator_name`.
