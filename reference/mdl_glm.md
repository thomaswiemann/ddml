# Wrapper for [`stats::glm()`](https://rdrr.io/r/stats/glm.html).

Simple wrapper for [`stats::glm()`](https://rdrr.io/r/stats/glm.html).

## Usage

``` r
mdl_glm(y, X, ...)
```

## Arguments

- y:

  The outcome variable.

- X:

  The feature matrix.

- ...:

  Additional arguments passed to `glm`. See
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) for a complete list
  of arguments.

## Value

`mdl_glm` returns an object of S3 class `mdl_glm` as a simple mask of
the return object of [`stats::glm()`](https://rdrr.io/r/stats/glm.html).

## See also

[`stats::glm()`](https://rdrr.io/r/stats/glm.html)

Other ml_wrapper:
[`mdl_glmnet()`](https://www.thomaswiemann.com/ddml/reference/mdl_glmnet.md),
[`mdl_ranger()`](https://www.thomaswiemann.com/ddml/reference/mdl_ranger.md),
[`mdl_xgboost()`](https://www.thomaswiemann.com/ddml/reference/mdl_xgboost.md),
[`ols()`](https://www.thomaswiemann.com/ddml/reference/ols.md)

## Examples

``` r
glm_fit <- mdl_glm(sample(0:1, 100, replace = TRUE),
                   matrix(rnorm(1000), 100, 10))
class(glm_fit)
#> [1] "mdl_glm" "glm"     "lm"     
```
