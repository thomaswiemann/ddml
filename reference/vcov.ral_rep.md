# Variance-Covariance Matrix for RAL Rep Objects

Variance-Covariance Matrix for RAL Rep Objects

## Usage

``` r
# S3 method for class 'ral_rep'
vcov(
  object,
  fit_idx = 1,
  aggregation = c("median", "mean", "spectral"),
  type = "HC1",
  ...
)
```

## Arguments

- object:

  An object inheriting from class `ral_rep`.

- fit_idx:

  Integer index of the fit. Defaults to 1.

- aggregation:

  Character string. Aggregation rule.

- type:

  Character. HC type. Default `"HC1"`.

- ...:

  Currently unused.

## Value

A \\p \times p\\ variance-covariance matrix.
