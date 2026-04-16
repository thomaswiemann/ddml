# Split a RAL Rep Object by Fit

Returns a named list of single-fit `ral_rep` objects. Each element
aggregates across resamples for a single ensemble type.

## Usage

``` r
# S3 method for class 'ral_rep'
as.list(x, ...)
```

## Arguments

- x:

  An object inheriting from class `ral_rep`.

- ...:

  Currently unused.

## Value

A named list of `ral_rep` objects.
