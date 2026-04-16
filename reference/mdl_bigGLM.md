# Wrapper for glmnet::bigGlm()

Simple wrapper for
[`glmnet::bigGlm()`](https://glmnet.stanford.edu/reference/bigGlm.html),
designed for sparse matrices.

## Usage

``` r
mdl_bigGlm(y, X, ...)
```

## Arguments

- y:

  The outcome variable.

- X:

  The (sparse) feature matrix.

- ...:

  Additional arguments passed to `bigGlm`. See
  [`glmnet::bigGlm()`](https://glmnet.stanford.edu/reference/bigGlm.html)
  for a complete list of arguments.

## Value

`mdl_bigGlm` returns an object of S3 class `mdl_bigGlm`.

## See also

[`glmnet::bigGlm()`](https://glmnet.stanford.edu/reference/bigGlm.html)

Other ml_wrapper:
[`mdl_glm()`](https://www.thomaswiemann.com/ddml/reference/mdl_glm.md),
[`mdl_glmnet()`](https://www.thomaswiemann.com/ddml/reference/mdl_glmnet.md),
[`mdl_ranger()`](https://www.thomaswiemann.com/ddml/reference/mdl_ranger.md),
[`mdl_xgboost()`](https://www.thomaswiemann.com/ddml/reference/mdl_xgboost.md),
[`ols()`](https://www.thomaswiemann.com/ddml/reference/ols.md)

## Examples

``` r
bigglm_fit <- mdl_bigGlm(rnorm(100), matrix(rnorm(1000), 100, 10))
class(bigglm_fit)
#> [1] "mdl_bigGlm" "list"      
```
