# Ordinary least squares.

Simple implementation of ordinary least squares that computes with
sparse feature matrices.

## Usage

``` r
ols(y, X, const = TRUE, w = NULL)
```

## Arguments

- y:

  The outcome variable.

- X:

  The feature matrix.

- const:

  Boolean equal to `TRUE` if a constant should be included.

- w:

  A vector of weights for weighted least squares.

## Value

`ols` returns an object of S3 class `ols`. An object of class `ols` is a
list containing the following components:

- `coef`:

  A vector with the regression coefficents.

- `y`, `X`, `const`, `w`:

  Pass-through of the user-provided arguments. See above.

## See also

Other ml_wrapper:
[`mdl_glm()`](https://www.thomaswiemann.com/ddml/reference/mdl_glm.md),
[`mdl_glmnet()`](https://www.thomaswiemann.com/ddml/reference/mdl_glmnet.md),
[`mdl_ranger()`](https://www.thomaswiemann.com/ddml/reference/mdl_ranger.md),
[`mdl_xgboost()`](https://www.thomaswiemann.com/ddml/reference/mdl_xgboost.md)

## Examples

``` r
ols_fit <- ols(rnorm(100), cbind(rnorm(100), rnorm(100)), const = TRUE)
ols_fit$coef
#>             [,1]
#> [1,] -0.06553227
#> [2,] -0.02382648
#> [3,] -0.11071182
```
