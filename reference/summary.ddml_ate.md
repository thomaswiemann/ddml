# Inference Methods for Treatment Effect Estimators.

Inference methods for treatment effect estimators. By default, standard
errors are heteroskedasiticty-robust. If the `ddml` estimator was
computed using a `cluster_variable`, the standard errors are also
cluster-robust by default.

## Usage

``` r
# S3 method for class 'ddml_ate'
summary(object, ...)

# S3 method for class 'ddml_att'
summary(object, ...)

# S3 method for class 'ddml_late'
summary(object, ...)
```

## Arguments

- object:

  An object of class `ddml_ate`, `ddml_att`, and `ddml_late`, as fitted
  by
  [`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
  [`ddml_att()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
  and
  [`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md),
  respectively.

- ...:

  Currently unused.

## Value

A matrix with inference results.

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]

# Estimate the average treatment effect using a single base learner, ridge.
ate_fit <- ddml_ate(y, D, X,
                    learners = list(what = mdl_glmnet,
                                    args = list(alpha = 0)),
                    sample_folds = 2,
                    silent = TRUE)
#> Warning: : 1 propensity scores were trimmed.
summary(ate_fit)
#> ATE estimation results: 
#>  
#>   Estimate Std. Error t value Pr(>|t|)
#>     -0.141     0.0153   -9.23  2.6e-20
```
