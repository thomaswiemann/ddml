# Estimator for the Partially Linear Model.

Estimator for the partially linear model.

## Usage

``` r
ddml_plm(
  y,
  D,
  X,
  learners,
  learners_DX = learners,
  sample_folds = 10,
  ensemble_type = "nnls",
  shortstack = FALSE,
  cv_folds = 10,
  custom_ensemble_weights = NULL,
  custom_ensemble_weights_DX = custom_ensemble_weights,
  cluster_variable = seq_along(y),
  subsamples = NULL,
  cv_subsamples_list = NULL,
  silent = FALSE
)
```

## Arguments

- y:

  The outcome variable.

- D:

  A matrix of endogenous variables.

- X:

  A (sparse) matrix of control variables.

- learners:

  May take one of two forms, depending on whether a single learner or
  stacking with multiple learners is used for estimation of the
  conditional expectation functions. If a single learner is used,
  `learners` is a list with two named elements:

  - `what` The base learner function. The function must be such that it
    predicts a named input `y` using a named input `X`.

  - `args` Optional arguments to be passed to `what`.

  If stacking with multiple learners is used, `learners` is a list of
  lists, each containing four named elements:

  - `fun` The base learner function. The function must be such that it
    predicts a named input `y` using a named input `X`.

  - `args` Optional arguments to be passed to `fun`.

  - `assign_X` An optional vector of column indices corresponding to
    control variables in `X` that are passed to the base learner.

  Omission of the `args` element results in default arguments being used
  in `fun`. Omission of `assign_X` results in inclusion of all variables
  in `X`.

- learners_DX:

  Optional argument to allow for different estimators of \\E\[D\|X\]\\.
  Setup is identical to `learners`.

- sample_folds:

  Number of cross-fitting folds.

- ensemble_type:

  Ensemble method to combine base learners into final estimate of the
  conditional expectation functions. Possible values are:

  - `"nnls"` Non-negative least squares.

  - `"nnls1"` Non-negative least squares with the constraint that all
    weights sum to one.

  - `"singlebest"` Select base learner with minimum MSPE.

  - `"ols"` Ordinary least squares.

  - `"average"` Simple average over base learners.

  Multiple ensemble types may be passed as a vector of strings.

- shortstack:

  Boolean to use short-stacking.

- cv_folds:

  Number of folds used for cross-validation in ensemble construction.

- custom_ensemble_weights:

  A numerical matrix with user-specified ensemble weights. Each column
  corresponds to a custom ensemble specification, each row corresponds
  to a base learner in `learners` (in chronological order). Optional
  column names are used to name the estimation results corresponding the
  custom ensemble specification.

- custom_ensemble_weights_DX:

  Optional argument to allow for different custom ensemble weights for
  `learners_DX`. Setup is identical to `custom_ensemble_weights`. Note:
  `custom_ensemble_weights` and `custom_ensemble_weights_DX` must have
  the same number of columns.

- cluster_variable:

  A vector of cluster indices.

- subsamples:

  List of vectors with sample indices for cross-fitting.

- cv_subsamples_list:

  List of lists, each corresponding to a subsample containing vectors
  with subsample indices for cross-validation.

- silent:

  Boolean to silence estimation updates.

## Value

`ddml_plm` returns an object of S3 class `ddml_plm`. An object of class
`ddml_plm` is a list containing the following components:

- `coef`:

  A vector with the \\\theta_0\\ estimates.

- `weights`:

  A list of matrices, providing the weight assigned to each base learner
  (in chronological order) by the ensemble procedure.

- `mspe`:

  A list of matrices, providing the MSPE of each base learner (in
  chronological order) computed by the cross-validation step in the
  ensemble construction.

- `ols_fit`:

  Object of class `lm` from the second stage regression of \\Y -
  \hat{E}\[Y\|X\]\\ on \\D - \hat{E}\[D\|X\]\\.

- `learners`,`learners_DX`,`cluster_variable`, `subsamples`,
  `cv_subsamples_list`, `ensemble_type`:

  Pass-through of selected user-provided arguments. See above.

## Details

`ddml_plm` provides a double/debiased machine learning estimator for the
parameter of interest \\\theta_0\\ in the partially linear model given
by

\\Y = \theta_0D + g_0(X) + U,\\

where \\(Y, D, X, U)\\ is a random vector such that \\E\[Cov(U, D\vert
X)\] = 0\\ and \\E\[Var(D\vert X)\] \neq 0\\, and \\g_0\\ is an unknown
nuisance function.

## References

Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). "Model Averaging
and Double Machine Learning." Journal of Applied Econometrics, 40(3):
249-269.

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). "Double/debiased machine learning for treatment and
structural parameters." The Econometrics Journal, 21(1), C1-C68.

Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2),
241-259.

## See also

[`summary.ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_plm.md)

Other ddml:
[`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md),
[`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md)

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]

# Estimate the partially linear model using a single base learner, ridge.
plm_fit <- ddml_plm(y, D, X,
                    learners = list(what = mdl_glmnet,
                                    args = list(alpha = 0)),
                    sample_folds = 2,
                    silent = TRUE)
summary(plm_fit)
#> PLM estimation results: 
#>  
#> , , single base learner
#> 
#>              Estimate Std. Error t value Pr(>|t|)
#> (Intercept) -0.000825    0.00689   -0.12 9.05e-01
#> D_r         -0.148560    0.01475  -10.08 7.11e-24
#> 

# Estimate the partially linear model using short-stacking with base learners
#     ols, lasso, and ridge. We can also use custom_ensemble_weights
#     to estimate the ATE using every individual base learner.
weights_everylearner <- diag(1, 3)
colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
plm_fit <- ddml_plm(y, D, X,
                    learners = list(list(fun = ols),
                                    list(fun = mdl_glmnet),
                                    list(fun = mdl_glmnet,
                                         args = list(alpha = 0))),
                    ensemble_type = 'nnls',
                    custom_ensemble_weights = weights_everylearner,
                    shortstack = TRUE,
                    sample_folds = 2,
                    silent = TRUE)
summary(plm_fit)
#> PLM estimation results: 
#>  
#> , , nnls
#> 
#>             Estimate Std. Error t value Pr(>|t|)
#> (Intercept)  0.00133    0.00688   0.193 8.47e-01
#> D_r         -0.14952    0.01473 -10.149 3.36e-24
#> 
#> , , mdl:ols
#> 
#>              Estimate Std. Error  t value Pr(>|t|)
#> (Intercept) -0.000266    0.00688  -0.0386 9.69e-01
#> D_r         -0.149559    0.01474 -10.1495 3.33e-24
#> 
#> , , mdl:lasso
#> 
#>              Estimate Std. Error  t value Pr(>|t|)
#> (Intercept) -0.000295    0.00688  -0.0429 9.66e-01
#> D_r         -0.149588    0.01473 -10.1538 3.19e-24
#> 
#> , , mdl:ridge
#> 
#>              Estimate Std. Error  t value Pr(>|t|)
#> (Intercept) -0.000247    0.00688  -0.0359 9.71e-01
#> D_r         -0.149493    0.01473 -10.1464 3.44e-24
#> 
```
