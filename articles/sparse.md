# Estimation with Sparse Matrices

## Introduction

`ddml` supports sparse matrices from the `Matrix` package by default.
This article illustrates Double/Debiased Machine Learning estimation
with sparse matrices using the prominent study of Angrist and Krueger
(1991) (AK91, hereafter) on returns to education.

``` r
library(ddml)
library(Matrix) # for sparse matrix operations
set.seed(900837)
```

## Data Construction

One of the coefficients of interest in AK91 is the effect of years of
education on the log weekly wage of American males born between
1930-1939. The authors instrument for years of schooling with quarter of
birth indicators (QOB). This strategy is motivated by mandatory
attendance laws that determine at what age children may drop out of
school. Since children born in later quarters achieve the age cut-off
after more years of schooling, higher QOB should imply more years of
schooling.

Although the data is quite large ($n = 329509$), a need for sparse
matrices only arises when the QOB instrument is interacted with other
variables. Popular control variables are place of birth (POB) and year
of birth (YOB). Depending on whether these are separately or jointly
interacted with QOB, this results in 180 and 1530 instruments,
respectively. The code snippet below constructs these two sets of
instruments as well as the matrix of controls. We use
[`sparse.model.matrix()`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)
to construct sparse matrix objects as supported by the `Matrix` package.

``` r
# Load data
AK91 <- readRDS("data/AK91.rds")

# Obtain instument matrix for 180 IVs
dat_IV180 <- sparse.model.matrix(~ YOB*POB +  QOB*YOB + QOB*POB, data = AK91)
which_instrument <- which(grepl("QOB", colnames(dat_IV180), fixed = TRUE))
Z_IV180 <- dat_IV180[, which_instrument]

# Obtain instument matrix for 1530 IVs
dat_IV1530 <- sparse.model.matrix(~ QOB * YOB * POB, data = AK91)
which_instrument <-  which(grepl("QOB", colnames(dat_IV1530), fixed = TRUE))
Z_IV1530 <- dat_IV1530[, which_instrument]

# Obtain control matrix
X <- dat_IV1530[, -which_instrument]
```

To illustrate the need for sparse matrices when working with the data of
AK91, consider differences in memory needed to hold the control matrix
alone. For the $329509 \times 510$ matrix of control variables, the
sparse version requires only 34.3Mb while the dense version requires
1.3Gb(=1300Mb)! The instrument matrix with 1530 instrument requires even
more space and regularly fails to load into memory.

``` r
# Memory needed for the sparse control matrix
format(object.size(X), units = "Mb")
#> [1] "34.3 Mb"

# Memory needed for the dense control matrix
format(object.size(as.matrix(X)), units = "Mb")
#> [1] "1302.3 Mb"
```

## Estimation with Sparse Matrices

The syntax for estimation with sparse matrices in `ddml` is *exactly*
the same as estimation with dense matrices. In the below, we replicate
AK91 using the simplified set of instruments and controls.

We begin with estimating the returns to schooling using the set of 180
instruments. Following the convention of the returns to
education-literature, our estimator selects only among the instruments
but does regularize the coefficients corresponding to the control
variables. This is achieved by formulating different base learners for
the first and second stage reduced forms (see
[`?ddml_fpliv`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md)),
and by setting the `penalty.factor` of the control variables to zero
(see
[`?mdl_glmnet`](https://www.thomaswiemann.com/ddml/reference/mdl_glmnet.md)).

``` r
learners_XZ <- list(list(what = ols),
                    list(what = mdl_glmnet,
                         args = list(cv = FALSE,
                                     penalty.factor = c(rep(0, 510),
                                                        rep(1, 180)))))

stacking_180IV_fit <- ddml_fpliv(y = AK91$LWKLYWGE, D = AK91$EDUC,
                                 Z = Z_IV180, X = X,
                                 learners = list(list(what = ols)),
                                 learners_DX = list(list(what = ols)),
                                 learners_DXZ = learners_XZ,
                                 ensemble_type = c("nnls1"),
                                 shortstack = TRUE,
                                 sample_folds = 2,
                                 silent = TRUE)
summary(stacking_180IV_fit)
#> DDML estimation: Flexible Partially Linear IV Model 
#> Obs: 329509   Folds: 2  Stacking: short-stack
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> D1           1.12e-01   2.15e-02    5.23  1.7e-07 ***
#> (Intercept) -6.14e-05   1.13e-03   -0.05     0.96    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The exercise can be repeated with the larger set of 1530 instruments as
well. Without support for sparse matrices, this estimation step would
not be possible without very large memory.

``` r
learners_XZ <- list(list(what = ols),
                    list(what = mdl_glmnet,
                         args = list(cv = FALSE,
                                     penalty.factor = c(rep(0, 510),
                                                        rep(1, 1530)))))

stacking_1530IV_fit <- ddml_fpliv(y = AK91$LWKLYWGE, D = AK91$EDUC,
                                 Z = Z_IV1530, X = X,
                                 learners = list(list(what = ols)),
                                 learners_DX = list(list(what = ols)),
                                 learners_DXZ = learners_XZ,
                                 ensemble_type = c("nnls1"),
                                 shortstack = TRUE,
                                 sample_folds = 2,
                                 silent = TRUE)
summary(stacking_1530IV_fit)
#> DDML estimation: Flexible Partially Linear IV Model 
#> Obs: 329509   Folds: 2  Stacking: short-stack
#> 
#>             Estimate Std. Error z value Pr(>|z|)
#> D1          0.046151   0.051943    0.89     0.37
#> (Intercept) 0.000064   0.001112    0.06     0.95
```

The coefficients corresponding to the two sets of instruments are quite
different. Leveraging the `ddml` functionality that allows for
specification of different sets of input variables, we construct a
stacking estimator that considers both first stage regressions
simultaneously.

``` r
# Construct column indices for combined control and instrument sets
Z_c <- cbind(Z_IV180, Z_IV1530); colnames(Z_c) <- 1:(180 + 1530)
set_IV180 <- 1:180; set_IV1530 <- 181:(180 + 1530)

learners_XZ <- list(list(what = ols,
                         assign_Z = set_IV180),
                    list(what = ols,
                          assign_Z = set_IV1530),
                    list(what = mdl_glmnet,
                         args = list(cv = FALSE,
                                     penalty.factor = c(rep(0, 510),
                                                        rep(1, 180))),
                          assign_Z = set_IV180),
                    list(what = mdl_glmnet,
                         args = list(cv = FALSE,
                                     penalty.factor = c(rep(0, 510),
                                                        rep(1, 1530))),
                         assign_Z = set_IV1530))
stacking_fit <- ddml_fpliv(y = AK91$LWKLYWGE, D = AK91$EDUC,
                           Z = Z_c, X = X,
                           learners = list(list(what = ols)),
                           learners_DX = list(list(what = ols)),
                           learners_DXZ = learners_XZ,
                           ensemble_type = c("nnls1"),
                           shortstack = TRUE,
                           sample_folds = 2,
                           silent = TRUE)
summary(stacking_fit)
#> DDML estimation: Flexible Partially Linear IV Model 
#> Obs: 329509   Folds: 2  Stacking: short-stack
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> D1           0.105040   0.015673    6.70  2.1e-11 ***
#> (Intercept) -0.000196   0.001125   -0.17     0.86    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The resulting coefficient is close to the coefficient based on 180
instruments. The stacking diagnostics confirm that the first stage
estimators with 180 instruments contribute almost exclusively to the
final estimate, suggesting that the expansion to 1530 instruments has
little benefit for the ols or lasso-based first stage fits.

``` r
diagnostics(stacking_fit)
#> Stacking diagnostics: Flexible Partially Linear IV Model 
#> Obs: 329509 
#> 
#>   y_X:
#>    learner   mspe     r2 weight_nnls1
#>  learner_1 0.4489 0.0258            1
#> 
#>   D1_X:
#>    learner   mspe     r2 weight_nnls1
#>  learner_1 10.187 0.0538            1
#> 
#>   D1_XZ:
#>    learner    mspe     r2 weight_nnls1
#>  learner_1 10.1851 0.0540       0.5874
#>  learner_2 10.2699 0.0461       0.0000
#>  learner_3 10.1852 0.0540       0.4126
#>  learner_4 10.2673 0.0464       0.0000
#>      nnls1 10.1850 0.0540           NA
#> 
#> Note: Ensemble MSPE and R2 for short-stacking rely on full-sample weights
#>        and represent in-sample fit over cross-fitted base predictions.
```

## References

Angrist J, Krueger A (1991). “Does Compulsory School Attendance Affect
Schooling and Earnings?” Quarterly Journal of Economics, 106(4),
979-1014.
