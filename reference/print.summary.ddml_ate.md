# Print Methods for Treatment Effect Estimators.

Print methods for treatment effect estimators.

## Usage

``` r
# S3 method for class 'summary.ddml_ate'
print(x, digits = 3, ...)

# S3 method for class 'summary.ddml_att'
print(x, digits = 3, ...)

# S3 method for class 'summary.ddml_late'
print(x, digits = 3, ...)
```

## Arguments

- x:

  An object of class `summary.ddml_ate`, `summary.ddml_att`, and
  `ddml_late`, as returned by
  [`summary.ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_ate.md),
  [`summary.ddml_att()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_ate.md),
  and
  [`summary.ddml_late()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_ate.md),
  respectively.

- digits:

  The number of significant digits used for printing.

- ...:

  Currently unused.

## Value

NULL.

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
#>     -0.142     0.0153   -9.24 2.37e-20
```
