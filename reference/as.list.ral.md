# Split a RAL Object by Fit

Returns a named list of single-fit `ral` objects, one per column of
`coefficients`. This is the primary mechanism for passing multi-ensemble
results to modelsummary.

## Usage

``` r
# S3 method for class 'ral'
as.list(x, ...)
```

## Arguments

- x:

  An object inheriting from class `ral`.

- ...:

  Currently unused.

## Value

A named list of `ral` objects, each with `nfit = 1`.

## See also

[`ral`](https://www.thomaswiemann.com/ddml/reference/ral.md),
[`lincom`](https://www.thomaswiemann.com/ddml/reference/lincom.md)
