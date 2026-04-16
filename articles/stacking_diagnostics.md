# Stacking Diagnostics and Cross-Validation Criteria

## Introduction

In Double/Debiased Machine Learning, the estimation of nuisance
parameters (conditional expectation functions) can be substantially
improved by adaptively combining multiple machine learners via stacking
or short-stacking.

While stacking is primarily used to minimize the out-of-sample mean
squared prediction error (MSPE) of the nuisance functions, it is often
insightful for researchers to peek under the hood: Which machine learner
performs the best? Are the differences in performance statistically
significant, or just noise from the random cross-fitting splits?

This vignette demonstrates how to extract rich diagnostic information
from `ddml` objects and how to perform statistical inference on base
learner performance using Cross-Validation Criteria (CVC) (Ahrens et
al., 2024).

*(Note: If you are interested in the computational benefits of the
short-stacking approach used here compared to standard stacking, see
`vignette("stacking")`.)*

## Setup and Estimation

For illustration, we use a random subsample of 5,000 observations from
the Angrist & Evans (1998) dataset. We will estimate the local average
treatment effect (LATE) of having more than two children (`morekids`) on
the mother’s labor supply (`worked`), using the sex of the first two
children (`samesex`) as an instrument. See
[`?AE98`](https://www.thomaswiemann.com/ddml/reference/AE98.md) for
details.

``` r
library(ddml)
set.seed(84210)

# Construct variables
y <- AE98[, "worked"]
D <- AE98[, "morekids"]
Z <- AE98[, "samesex"]
X <- AE98[, c("age", "agefst", "black", "hisp", "othrace", "educ")]
```

First, we estimate the LATE using a combination of linear regression,
lasso, and gradient boosting. We use the `ensemble_type = "nnls1"`
scheme, which constructs a convex combination of the base learners
(non-negative weights that sum to one) to minimize the MSPE.

``` r
# Specify multiple base learners
learners <- list(
  list(what = ols),
  list(what = mdl_glmnet),
  list(what = mdl_xgboost, args = list(nrounds = 10, max_depth = 1)),
  list(what = mdl_xgboost, args = list(nrounds = 100, max_depth = 2))
)

# Estimate LATE with short-stacking
late_fit <- ddml_late(
  y, D, Z, X, 
  learners = learners,
  ensemble_type = "nnls1",
  shortstack = TRUE,
  sample_folds = 5,
  silent = TRUE
)
```

## Basic Diagnostics

The
[`diagnostics()`](https://www.thomaswiemann.com/ddml/reference/diagnostics.md)
function extracts the out-of-sample performance of each base learner for
every conditional expectation function (reduced form) estimated during
the procedure.

``` r
# Extract stacking diagnostics
diag <- diagnostics(late_fit)
print(diag)
#> Stacking diagnostics: Local Average Treatment Effect 
#> Obs: 5000 
#> 
#>   y_X_Z0:
#>    learner   mspe     r2 weight_nnls1
#>  learner_1 0.2408 0.0316        0.000
#>  learner_2 0.2407 0.0319        0.901
#>  learner_3 0.2446 0.0165        0.000
#>  learner_4 0.2470 0.0067        0.099
#>      nnls1 0.2407 0.0322           NA
#> 
#>   y_X_Z1:
#>    learner   mspe     r2 weight_nnls1
#>  learner_1 0.2420 0.0309       0.0000
#>  learner_2 0.2419 0.0314       0.8775
#>  learner_3 0.2459 0.0155       0.0874
#>  learner_4 0.2485 0.0049       0.0351
#>      nnls1 0.2419 0.0316           NA
#> 
#>   D_X_Z0:
#>    learner   mspe     r2 weight_nnls1
#>  learner_1 0.2095 0.0776       0.6868
#>  learner_2 0.2095 0.0776       0.0000
#>  learner_3 0.2150 0.0533       0.0693
#>  learner_4 0.2130 0.0621       0.2439
#>      nnls1 0.2090 0.0799           NA
#> 
#>   D_X_Z1:
#>    learner   mspe     r2 weight_nnls1
#>  learner_1 0.2239 0.0797       0.7212
#>  learner_2 0.2239 0.0797       0.0000
#>  learner_3 0.2297 0.0560       0.1049
#>  learner_4 0.2288 0.0594       0.1739
#>      nnls1 0.2235 0.0813           NA
#> 
#>   Z_X:
#>    learner   mspe      r2 weight_nnls1
#>  learner_1 0.2503 -0.0015            0
#>  learner_2 0.2498  0.0002            1
#>  learner_3 0.2501 -0.0010            0
#>  learner_4 0.2551 -0.0209            0
#>      nnls1 0.2498  0.0002           NA
#> 
#> Note: Ensemble MSPE and R2 for short-stacking rely on full-sample weights
#>        and represent in-sample fit over cross-fitted base predictions.
```

For each reduced form equation (e.g.,
$E\left\lbrack Y|Z = 0,X \right\rbrack$), the diagnostics table
reports: - **Learner:** The specific base learner. - **MSPE:** The
out-of-sample mean squared prediction error. - **R2:** The out-of-sample
$R^{2}$, computed as
$1 - \text{MSPE}/\text{Var}\left( \text{Target} \right)$. - **Weight:**
The weight assigned to the learner by the chosen `ensemble_type`.

For workflows integrated with the `broom` package, or when preparing
publication-ready tables with `modelsummary`, the
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) method
converts the rich diagnostic tables into a flat `data.frame`.

``` r
head(tidy(diag))
#>   equation   learner      mspe          r2 weight_nnls1
#> 1   y_X_Z0 learner_1 0.2408175 0.031591434 4.234957e-17
#> 2   y_X_Z0 learner_2 0.2407379 0.031911480 9.010382e-01
#> 3   y_X_Z0 learner_3 0.2445764 0.016475863 0.000000e+00
#> 4   y_X_Z0 learner_4 0.2470071 0.006701109 9.896183e-02
#> 5   y_X_Z0     nnls1 0.2406614 0.032219302           NA
#> 6   y_X_Z1 learner_1 0.2420426 0.030906416 0.000000e+00
```

## Inference on Base Learners (Cross-Validation Criteria)

It is common for researchers to wonder if the base learner assigned the
highest weight is *statistically significantly* better than the
alternatives. Because cross-fitting creates randomly partitioned
hold-out sets, differences in MSPE could simply be due to sample
splitting noise.

`ddml` implements the Cross-Validation Criteria (CVC) methodology to
formally test for equal predictive ability across base learners. Setting
`cvc = TRUE` in the
[`diagnostics()`](https://www.thomaswiemann.com/ddml/reference/diagnostics.md)
function engages a multiplier bootstrap on the out-of-sample fold
residuals to test whether each learner performs equally well as the
single-best learner.

``` r
# Run diagnostics with Cross-Validation Criteria tests
diag_cvc <- diagnostics(late_fit, cvc = TRUE, bootnum = 1000)
print(diag_cvc)
#> Stacking diagnostics: Local Average Treatment Effect 
#> Obs: 5000 
#> 
#>   y_X_Z0:
#>    learner   mspe     r2 weight_nnls1 cvc_pval
#>  learner_1 0.2408 0.0316        0.000    0.505
#>  learner_2 0.2407 0.0319        0.901    0.996
#>  learner_3 0.2446 0.0165        0.000    0.000
#>  learner_4 0.2470 0.0067        0.099    0.001
#>      nnls1 0.2407 0.0322           NA       NA
#> 
#>   y_X_Z1:
#>    learner   mspe     r2 weight_nnls1 cvc_pval
#>  learner_1 0.2420 0.0309       0.0000    0.221
#>  learner_2 0.2419 0.0314       0.8775    1.000
#>  learner_3 0.2459 0.0155       0.0874    0.004
#>  learner_4 0.2485 0.0049       0.0351    0.000
#>      nnls1 0.2419 0.0316           NA       NA
#> 
#>   D_X_Z0:
#>    learner   mspe     r2 weight_nnls1 cvc_pval
#>  learner_1 0.2095 0.0776       0.6868    0.833
#>  learner_2 0.2095 0.0776       0.0000    0.875
#>  learner_3 0.2150 0.0533       0.0693    0.000
#>  learner_4 0.2130 0.0621       0.2439    0.018
#>      nnls1 0.2090 0.0799           NA       NA
#> 
#>   D_X_Z1:
#>    learner   mspe     r2 weight_nnls1 cvc_pval
#>  learner_1 0.2239 0.0797       0.7212    0.879
#>  learner_2 0.2239 0.0797       0.0000    0.809
#>  learner_3 0.2297 0.0560       0.1049    0.000
#>  learner_4 0.2288 0.0594       0.1739    0.006
#>      nnls1 0.2235 0.0813           NA       NA
#> 
#>   Z_X:
#>    learner   mspe      r2 weight_nnls1 cvc_pval
#>  learner_1 0.2503 -0.0015            0    0.076
#>  learner_2 0.2498  0.0002            1    0.994
#>  learner_3 0.2501 -0.0010            0    0.253
#>  learner_4 0.2551 -0.0209            0    0.000
#>      nnls1 0.2498  0.0002           NA       NA
#> 
#> Note: CVC compares individual base learners.
#>        Shortstacked ensemble CVC is not available.
#> Note: Ensemble MSPE and R2 for short-stacking rely on full-sample weights
#>        and represent in-sample fit over cross-fitted base predictions.
```

A new column is added to the table: - **CVC p-value:** The p-value
testing the null hypothesis that the learner has the same true MSPE as
the absolute best-performing learner. A small p-value (e.g., $< 0.05$)
indicates that the learner performs significantly worse than the best
learner. Conversely, a large p-value means the learner belongs to the
“model confidence set” of learners whose performance is not
statistically distinguishable from the best.

## References

Ahrens A, Chernozhukov V, Hansen C B, Kozbur D, Schaffer M E, Wiemann T
(2026). “An Introduction to Double/Debiased Machine Learning.” *Journal
of Economic Literature*, forthcoming.

Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). “Model Averaging
and Double Machine Learning.” *Journal of Applied Econometrics*, 40(3):
249-269.

Angrist J, Evans W (1998). “Children and Their Parents’ Labor Supply:
Evidence from Exogenous Variation in Family Size.” *American Economic
Review*, 88(3), 450-477.

Lei J (2020). “Cross-Validation With Confidence.” *Journal of the
American Statistical Association*, 115(532), 1978-1997.
