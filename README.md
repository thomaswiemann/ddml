
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ddml

<!-- badges: start -->

[![R-CMD-check](https://github.com/thomaswiemann/ddml/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/thomaswiemann/ddml/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/thomaswiemann/ddml/branch/master/graph/badge.svg?token=PHB9W2TJ6S)](https://codecov.io/gh/thomaswiemann/ddml)
[![CodeFactor](https://www.codefactor.io/repository/github/thomaswiemann/ddml/badge)](https://www.codefactor.io/repository/github/thomaswiemann/ddml)
<!-- badges: end -->

`ddml` is an implementation of double/debiased machine learning
estimators targeting common causal parameters such as the average
treatment effect. The key benefit of the package is the straightforward
estimation of nuisance parameters using (short-)stacking, which
simultaneously leverages multiple base learners to increase robustness
to the underlying data generating process.

`ddml` is the sister R package to our
[Stata](https://github.com/aahrens1/ddml/) package, mirroring its key
features while also leveraging R to simplify estimation with
user-provided learners and/or sparse matrices. See also [Ahrens et
al. (2023)](https://arxiv.org/abs/2301.09397) with additional discussion
of the supported causal models and benefits of (short)-stacking.

## Installation

Install latest development version from GitHub (requires
[devtools](https://github.com/hadley/devtools) package):

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("thomaswiemann/ddml", dependencies = TRUE)
```

## Example using the Data from Angrist & Evans (1998)

To illustrate `ddml` on a simple example, we can use the included random
subsample of 5,000 observations from the data of Angrist & Evans (1998).
The data contains information on the labor supply of mothers, their
children, as well as demographic data. See `?AE98` for details.

The goal of the exercise is to understand how an increase in the number
of children causes a change in the employment status of mothers
(`worked`). For this purpose, we consider the strategy of Angrist &
Evans (1998): To instrument for having more than two children
(`morekids`), the authors use the `samesex` variable that indicates
whether the first two children are either both male or female. The idea
is that families with two same-sex children are more likely to have
preferences for a third child (of the opposite sex). Additional
controls, including the age and years of education of the mother, are
included to address concerns about endogeneity of the instrument.

``` r
# Load ddml and set seet
library(ddml)
set.seed(75523)

# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
Z = AE98[, "samesex"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
```

### LATE Estimation with Lasso

`ddml_late` estimates the local average treatment effect (LATE) using
double/debiased machine learning under the LATE assumptions of Imbens &
Angrist (1994). Double/debiased machine learning uses generic machine
learners for the high dimensional nuisance parameters that arise in
estimation of the LATE. For example, the below snippet estimates the
LATE with a double/debiased machine learning estimator that leverages
gradient boosting to estimate the nuisance parameters.

``` r
# Estimate the local average treatment effect using xgboost.
time_single <- system.time({
  late_fit <- ddml_late(y, D, Z, X,
                        learners = list(what = mdl_xgboost,
                                        args = list(nrounds = 100,
                                                    max_depth = 1)),
                        sample_folds = 10,
                        silent = TRUE)
})

cat("LATE estimate: ", late_fit$late)
#> LATE estimate:  -0.1888416
```

### LATE Estimation with (Short-)Stacking

A key feature of `ddml` is that the nuisance parameter estimates can be
constructed using multiple machine learners simultaneously. Since the
statistical properties of machine learners depend heavily on the
underlying (unknown!) structure of the data, adaptive combination of
multiple machine learners can increase robustness. Below, `ddml_late`
estimates the LATE with shortstacking (a computationally convenient
variant of stacking) based on three base learners: linear regression
(`ols`), lasso (`mdl_glmnet`), and gradient boosting (`mdl_xgboost`).

``` r
# Estimate the local average treatment effect using short-stacking with base
#     learners ols, rlasso, and xgboost.
time_shortstacking <- system.time({
  late_fit_short <- ddml_late(y, D, Z, X,
                              learners = list(list(fun = ols),
                                              list(fun = mdl_glmnet),
                                              list(fun = mdl_xgboost,
                                                   args = list(nrounds = 100,
                                                               max_depth = 1))),
                              ensemble_type = 'nnls1',
                              shortstack = TRUE,
                              sample_folds = 10,
                              silent = TRUE)
})

cat("LATE estimate: ", late_fit_short$late)
#> LATE estimate:  -0.206133
```

Short-stacking combines the considered base learners to minimize the
out-of-sample mean squared prediction error in each reduced form
equation (e.g., to estimate E\[Y\|Z=0,X\]). Setting
`ensemble_type = 'nnls1` ensures that the weights associated with each
base learner are in the unit simplex. The below snippet shows the weight
associated with each of the three base learners (in chronological order)
in the reduced form of E\[Y\|Z=0,X\]: `ols` gets a weight of 0.132,
`mdl_glmnet` gets a weight of 0.564, and `mdl_xgboost` gets a weight of
0.304.

``` r
t(late_fit_short$weights$y_X_Z0)
#>            [,1]    [,2]      [,3]
#> nnls1 0.1315521 0.56452 0.3039279
```

Despite the additional computational burden from considering multiple
base learners and optimally combining their predictions, the
short-stacked LATE estimator does not take substantially longer to
compute.

``` r
cat("Time using boosting: ", time_single[1])
#> Time using boosting:  10.04
cat("Time using shortstacking: ", time_shortstacking[1])
#> Time using shortstacking:  14.29
```

The low additional computational burden is thanks to short-stacking
leveraging the cross-predictions arising in double/debiased machine
learning estimators. A simple stacking estimator that used
cross-validation to construct combinations of base learners in every
sample fold increases computation time by more than 10-fold.

``` r
# Estimate the local average treatment effect using short-stacking with base
#     learners ols, rlasso, and xgboost.
time_stacking <- system.time({
  late_fit_stack <- ddml_late(y, D, Z, X,
                              learners = list(list(fun = ols),
                                              list(fun = mdl_glmnet),
                                              list(fun = mdl_xgboost,
                                                   args = list(nrounds = 100,
                                                               max_depth = 1))),
                              ensemble_type = 'nnls1',
                              shortstack = FALSE,
                              sample_folds = 10,
                              cv_folds = 10,
                              silent = TRUE)
})

cat("LATE estimate: ", late_fit_stack$late)
#> LATE estimate:  -0.213941
cat("Time using stacking: ", time_stacking[1])
#> Time using stacking:  164.34
```

## References

Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2023). “ddml:
Double/debiased machine learning in Stata.”
<https://arxiv.org/abs/2301.09397>

Angrist J, Evans W, (1998). “Children and Their Parents’ Labor Supply:
Evidence from Exogenous Variation in Family Size.” American Economic
Review, 88(3), 450-477.

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). “Double/debiased machine learning for treatment and
structural parameters.” The Econometrics Journal, 21(1), C1-C68.

Imbens G, Angrist J (1004). “Identification and Estimation of Local
Average Treatment Effects.” Econometrica, 62(2), 467-475.

Wolpert D H (1992). “Stacked generalization.” Neural Networks, 5(2),
241-259.
