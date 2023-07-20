
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
estimation of nuisance parameters using stacking, which simultaneously
leverages multiple base learners to increase robustness to the
underlying data generating process.

`ddml` is the sister R package to our
[Stata](https://github.com/aahrens1/ddml/) package, mirroring its key
features while also leveraging R to simplify estimation with
user-provided learners and/or sparse matrices (see
[Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)).

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
# Load ddml
library(ddml)

# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
Z = AE98[, "samesex"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
```

The causal model under consideration is the interactive model in which
  
![Y = g\_0(D, X) +
U,](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%20%3D%20g_0%28D%2C%20X%29%20%2B%20U%2C
"Y = g_0(D, X) + U,")  
where ![(Y, D, X, Z,
U)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28Y%2C%20D%2C%20X%2C%20Z%2C%20U%29
"(Y, D, X, Z, U)") is a random vector with
![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y
"Y") denoting the outcome,
![D](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;D
"D") denoting the binary variable of interest,
![Z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Z
"Z") denoting the binary instrument,
![X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X
"X") denoting the vector of controls, and
![U](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;U
"U") denoting all other determinants of
![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y
"Y") other than ![(D, X,
Z)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28D%2C%20X%2C%20Z%29
"(D, X, Z)"). In this model, the local average treatment effect (LATE)
is defined as   
![\\theta\_0^{\\textrm{LATE}} \\equiv E\[g\_0(1, X) - g\_0(0, X)\\vert
p\_0(1, X) \> p(0,
X)\],](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_0%5E%7B%5Ctextrm%7BLATE%7D%7D%20%5Cequiv%20%20E%5Bg_0%281%2C%20X%29%20-%20g_0%280%2C%20X%29%5Cvert%20p_0%281%2C%20X%29%20%3E%20p%280%2C%20X%29%5D%2C
"\\theta_0^{\\textrm{LATE}} \\equiv  E[g_0(1, X) - g_0(0, X)\\vert p_0(1, X) \> p(0, X)],")  
where ![p\_0(Z, X) \\equiv \\Pr(D=1\\vert Z,
X)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_0%28Z%2C%20X%29%20%5Cequiv%20%5CPr%28D%3D1%5Cvert%20Z%2C%20X%29
"p_0(Z, X) \\equiv \\Pr(D=1\\vert Z, X)") denotes the propensity score.

`ddml_late` estimates the LATE using double/debiased machine learning
under the LATE assumptions of Imbens & Angrist (1994). Double/debiased
machine learning uses generic machine learners for the high dimensional
nuisance parameters that arise in estimation of the LATE. `ddml` allows
for combination of multiple machine learners using stacking to to
increase robustness to the underlying data generating process.

Below, `ddml_late` estimates the LATE with shortstacking (a
computationally convenient variant of stacking) based on three base
learners: linear regression (`ols`), lasso (`mdl_glmnet`), and gradient
boosting (`mdl_xgboost`).

``` r
# Estimate the local average treatment effect using short-stacking with base
#     learners ols, rlasso, and xgboost.
late_fit <- ddml_late(y, D, Z, X,
                      learners = list(list(fun = ols),
                                      list(fun = mdl_glmnet),
                                      list(fun = mdl_xgboost,
                                         args = list(nrounds = 500,
                                                     max_depth = 3))),
                      ensemble_type = 'nnls',
                      shortstack = TRUE,
                      sample_folds = 5,
                      silent = TRUE)
late_fit$late
#>      nnls 
#> -0.210019
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

Wolpert D H (1992). “Stacked generalization.” Neural Networks, 5(2),
241-259.
