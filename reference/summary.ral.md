# Summary for RAL Estimators

Computes a coefficient table with estimates, standard errors, z-values,
and p-values.

## Usage

``` r
# S3 method for class 'ral'
summary(object, type = "HC1", ...)

# S3 method for class 'summary.ral'
print(x, digits = 3, ...)
```

## Arguments

- object:

  An object inheriting from class `ral`.

- type:

  Character. HC type. Default `"HC1"`.

- ...:

  Currently unused.

- x:

  An object of class `summary.ral`.

- digits:

  Number of significant digits. Default 3.

## Value

An object of class `summary.ral` with:

- `coefficients`:

  A 3-dimensional array (\\p \times 4 \times\\ nfit).

- `type`:

  The HC type used.

- `nobs`:

  Number of observations.
