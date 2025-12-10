# Wrapper for [`xgboost::xgboost()`](https://rdrr.io/pkg/xgboost/man/xgboost.html).

Simple wrapper for
[`xgboost::xgboost()`](https://rdrr.io/pkg/xgboost/man/xgboost.html)
with some changes to the default arguments.

## Usage

``` r
mdl_xgboost(y, X, nrounds = 500, verbosity = 0, ...)
```

## Arguments

- y:

  The outcome variable.

- X:

  The (sparse) feature matrix.

- nrounds:

  Number of boosting iterations / rounds.

  Note that the number of default boosting rounds here is not
  automatically tuned, and different problems will have vastly different
  optimal numbers of boosting rounds.

- verbosity:

  Verbosity of printing messages. Valid values of 0 (silent), 1
  (warning), 2 (info), and 3 (debug).

- ...:

  Additional arguments passed to `xgboost`. See
  [`xgboost::xgboost()`](https://rdrr.io/pkg/xgboost/man/xgboost.html)
  for a complete list of arguments.

## Value

`mdl_xgboost` returns an object of S3 class `mdl_xgboost` as a simple
mask to the return object of
[`xgboost::xgboost()`](https://rdrr.io/pkg/xgboost/man/xgboost.html).

## References

Chen T, Guestrin C (2011). "Xgboost: A Scalable Tree Boosting System."
Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge
Discovery and Data Mining, 785–794.

## See also

[`xgboost::xgboost()`](https://rdrr.io/pkg/xgboost/man/xgboost.html)

Other ml_wrapper:
[`mdl_glm()`](https://www.thomaswiemann.com/ddml/reference/mdl_glm.md),
[`mdl_glmnet()`](https://www.thomaswiemann.com/ddml/reference/mdl_glmnet.md),
[`mdl_ranger()`](https://www.thomaswiemann.com/ddml/reference/mdl_ranger.md),
[`ols()`](https://www.thomaswiemann.com/ddml/reference/ols.md)

## Examples

``` r
xgboost_fit <- mdl_xgboost(rnorm(50), matrix(rnorm(150), 50, 3),
                           nrounds = 1)
class(xgboost_fit)
#> [1] "mdl_xgboost" "xgboost"     "xgb.Booster"
```
