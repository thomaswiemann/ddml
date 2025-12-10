# Wrapper for [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).

Simple wrapper for
[`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).
Supports regression (default) and probability forests (set
`probability = TRUE`).

## Usage

``` r
mdl_ranger(y, X, ...)
```

## Arguments

- y:

  The outcome variable.

- X:

  The feature matrix.

- ...:

  Additional arguments passed to `ranger`. See
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md)
  for a complete list of arguments.

## Value

`mdl_ranger` returns an object of S3 class `ranger` as a simple mask of
the return object of
[`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).

## References

Wright M N, Ziegler A (2017). "ranger: A fast implementation of random
forests for high dimensional data in C++ and R." Journal of Statistical
Software 77(1), 1-17.

## See also

[`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md)

Other ml_wrapper:
[`mdl_glm()`](https://www.thomaswiemann.com/ddml/reference/mdl_glm.md),
[`mdl_glmnet()`](https://www.thomaswiemann.com/ddml/reference/mdl_glmnet.md),
[`mdl_xgboost()`](https://www.thomaswiemann.com/ddml/reference/mdl_xgboost.md),
[`ols()`](https://www.thomaswiemann.com/ddml/reference/ols.md)

## Examples

``` r
ranger_fit <- mdl_ranger(rnorm(100), matrix(rnorm(1000), 100, 10))
class(ranger_fit)
#> [1] "mdl_ranger" "ranger"    
```
