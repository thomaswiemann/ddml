# Extract Aggregated Coefficients

Extract Aggregated Coefficients

## Usage

``` r
# S3 method for class 'ral_rep'
coef(object, aggregation = c("median", "mean", "spectral"), ...)
```

## Arguments

- object:

  An object inheriting from class `ral_rep`.

- aggregation:

  Character string: `"median"` (default), `"mean"`, or `"spectral"`.

- ...:

  Currently unused.

## Value

Named vector (single fit) or matrix (multiple).
