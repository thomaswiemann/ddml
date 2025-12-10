# Wrapper for [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html).

Simple wrapper for
[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
and
[`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html).

## Usage

``` r
mdl_glmnet(y, X, cv = TRUE, ...)
```

## Arguments

- y:

  The outcome variable.

- X:

  The (sparse) feature matrix.

- cv:

  Boolean to indicate use of lasso with cross-validated penalty.

- ...:

  Additional arguments passed to `glmnet`. See
  [`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
  and
  [`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)
  for a complete list of arguments.

## Value

`mdl_glmnet` returns an object of S3 class `mdl_glmnet` as a simple mask
of the return object of
[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html)
or
[`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html).

## References

Friedman J, Hastie T, Tibshirani R (2010). "Regularization Paths for
Generalized Linear Models via Coordinate Descent." Journal of
Statistical Software, 33(1), 1–22.

Simon N, Friedman J, Hastie T, Tibshirani R (2011). "Regularization
Paths for Cox's Proportional Hazards Model via Coordinate Descent."
Journal of Statistical Software, 39(5), 1–13.

## See also

[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html),[`glmnet::cv.glmnet()`](https://glmnet.stanford.edu/reference/cv.glmnet.html)

Other ml_wrapper:
[`mdl_glm()`](https://www.thomaswiemann.com/ddml/reference/mdl_glm.md),
[`mdl_ranger()`](https://www.thomaswiemann.com/ddml/reference/mdl_ranger.md),
[`mdl_xgboost()`](https://www.thomaswiemann.com/ddml/reference/mdl_xgboost.md),
[`ols()`](https://www.thomaswiemann.com/ddml/reference/ols.md)

## Examples

``` r
glmnet_fit <- mdl_glmnet(rnorm(100), matrix(rnorm(1000), 100, 10))
class(glmnet_fit)
#> [1] "mdl_glmnet" "cv.glmnet" 
```
