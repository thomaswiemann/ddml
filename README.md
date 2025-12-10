
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ddml

<!-- badges: start -->

[![R-CMD-check](https://github.com/thomaswiemann/ddml/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/thomaswiemann/ddml/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/thomaswiemann/ddml/branch/master/graph/badge.svg?token=PHB9W2TJ6S)](https://app.codecov.io/gh/thomaswiemann/ddml)
[![CodeFactor](https://www.codefactor.io/repository/github/thomaswiemann/ddml/badge)](https://www.codefactor.io/repository/github/thomaswiemann/ddml)
[![CRAN
Version](https://www.r-pkg.org/badges/version/ddml)](https://cran.r-project.org/package=ddml)
[![cranlogs](https://cranlogs.r-pkg.org/badges/ddml)](https://cran.r-project.org/package=ddml)
<!-- badges: end -->

`ddml` is an implementation of double/debiased machine learning
estimators as proposed by Chernozhukov et al. (2018). The key feature of
`ddml` is the straightforward estimation of nuisance parameters using
(short-)stacking (Wolpert, 1992), which allows for multiple machine
learners to increase robustness to the underlying data generating
process. See also [Ahrens et
al. (2024b)](https://arxiv.org/abs/2401.01645) for a detailed
illustration of the practical benefits of combining DDML with
(short-)stacking.

`ddml` is the sister R package to our
[Stata](https://github.com/aahrens1/ddml/) package, mirroring its key
features while also leveraging R to simplify estimation with
user-provided machine learners and/or sparse matrices. See also [Ahrens
et al. (2024a)](https://arxiv.org/abs/2301.09397) with additional
discussion of the supported causal models and benefits of
(short)-stacking.

## Installation

Install the latest development version from GitHub (requires
[devtools](https://github.com/r-lib/devtools) package):

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("thomaswiemann/ddml", dependencies = TRUE)
```

Install the latest public release from CRAN:

``` r
install.packages("ddml")
```

## Example: LATE Estimation based on (Short-)Stacking

To illustrate `ddml` on a simple example, consider the included random
subsample of 5,000 observations from the data of Angrist & Evans (1998).
The data contains information on the labor supply of mothers, their
children, as well as demographic data. See `?AE98` for details.

``` r
# Load ddml and set seed
library(ddml)
set.seed(75523)

# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
Z = AE98[, "samesex"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
```

`ddml_late` estimates the local average treatment effect (LATE) using
double/debiased machine learning (see `?ddml_late`). Since the
statistical properties of machine learners depend heavily on the
underlying (unknown!) structure of the data, adaptive combination of
multiple machine learners can increase robustness. In the below snippet,
`ddml_late` estimates the LATE with short-stacking based on three base
learners:

- linear regression (see `?ols`)
- lasso (see `?mdl_glmnet`)
- gradient boosting (see `?mdl_xgboost`)

``` r
# Estimate the local average treatment effect using short-stacking with base
#     learners ols, rlasso, and xgboost.
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
summary(late_fit_short)
#> LATE estimation results: 
#>  
#>       Estimate Std. Error t value Pr(>|t|)
#> nnls1   -0.221      0.187   -1.18    0.236
```

## Learn More about `ddml`

Check out our articles to learn more:

- `vignette("ddml")` is a more detailed introduction to `ddml`
- `vignette("stacking")` discusses computational benefits of
  short-stacking
- `vignette("new_ml_wrapper")` shows how to write user-provided base
  learners
- `vignette("sparse")` illustrates support of sparse matrices (see
  `?Matrix`)
- `vignette("did")` discusses integration with the diff-in-diff package
  [`did`](https://bcallaway11.github.io/did/)

For additional applied examples, see our case studies:

- `vignette("example_401k")` revisits the effect of 401k participation
  on retirement savings
- `vignette("example_BLP95")` considers flexible demand estimation with
  endogenous prices

## Other Double/Debiased Machine Learning Packages

`ddml` is built to easily (and quickly) estimate common causal
parameters with multiple machine learners. With its support for
short-stacking, sparse matrices, and easy-to-learn syntax, we hope
`ddml` is a useful complement to
[`DoubleML`](https://docs.doubleml.org/stable/index.html), the expansive
R and Python package.
[`DoubleML`](https://docs.doubleml.org/stable/index.html) supports many
advanced features such as [multiway
clustering](https://docs.doubleml.org/stable/examples/R_double_ml_multiway_cluster.html)
and
[stacking](https://docs.doubleml.org/stable/examples/R_double_ml_pipeline.html).

## References

Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024a). “ddml:
Double/debiased machine learning in Stata.” Stata Journal, 24(1): 3-45.

Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024b). “Model averaging
and double machine learning.” Journal of Applied Econometrics, 40(3): 249-269.

Angrist J, Evans W, (1998). “Children and Their Parents’ Labor Supply:
Evidence from Exogenous Variation in Family Size.” American Economic
Review, 88(3), 450-477.

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). “Double/debiased machine learning for treatment and
structural parameters.” The Econometrics Journal, 21(1), C1-C68.

Wolpert D H (1992). “Stacked generalization.” Neural Networks, 5(2),
241-259.
