# Estimator of the Local Average Treatment Effect.

Estimator of the local average treatment effect.

## Usage

``` r
ddml_late(
  y,
  D,
  Z,
  X,
  learners,
  learners_DXZ = learners,
  learners_ZX = learners,
  sample_folds = 10,
  ensemble_type = "nnls",
  shortstack = FALSE,
  cv_folds = 10,
  custom_ensemble_weights = NULL,
  custom_ensemble_weights_DXZ = custom_ensemble_weights,
  custom_ensemble_weights_ZX = custom_ensemble_weights,
  cluster_variable = seq_along(y),
  subsamples_byZ = NULL,
  cv_subsamples_byZ = NULL,
  trim = 0.01,
  silent = FALSE
)
```

## Arguments

- y:

  The outcome variable.

- D:

  The binary endogenous variable of interest.

- Z:

  Binary instrumental variable.

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

  - `assign_Z` An optional vector of column indices corresponding to
    instruments in `Z` that are passed to the base learner.

  Omission of the `args` element results in default arguments being used
  in `fun`. Omission of `assign_X` (and/or `assign_Z`) results in
  inclusion of all variables in `X` (and/or `Z`).

- learners_DXZ, learners_ZX:

  Optional arguments to allow for different estimators of \\E\[D \vert
  X, Z\]\\, \\E\[Z \vert X\]\\. Setup is identical to `learners`.

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

- custom_ensemble_weights_DXZ, custom_ensemble_weights_ZX:

  Optional arguments to allow for different custom ensemble weights for
  `learners_DXZ`,`learners_ZX`. Setup is identical to
  `custom_ensemble_weights`. Note: `custom_ensemble_weights` and
  `custom_ensemble_weights_DXZ`,`custom_ensemble_weights_ZX` must have
  the same number of columns.

- cluster_variable:

  A vector of cluster indices.

- subsamples_byZ:

  List of two lists corresponding to the two instrument levels. Each
  list contains vectors with sample indices for cross-fitting.

- cv_subsamples_byZ:

  List of two lists, each corresponding to one of the two instrument
  levels. Each of the two lists contains lists, each corresponding to a
  subsample and contains vectors with subsample indices for
  cross-validation.

- trim:

  Number in (0, 1) for trimming the estimated propensity scores at
  `trim` and `1-trim`.

- silent:

  Boolean to silence estimation updates.

## Value

`ddml_late` returns an object of S3 class `ddml_late`. An object of
class `ddml_late` is a list containing the following components:

- `late`:

  A vector with the average treatment effect estimates.

- `weights`:

  A list of matrices, providing the weight assigned to each base learner
  (in chronological order) by the ensemble procedure.

- `mspe`:

  A list of matrices, providing the MSPE of each base learner (in
  chronological order) computed by the cross-validation step in the
  ensemble construction.

- `psi_a`, `psi_b`:

  Matrices needed for the computation of scores. Used in
  [`summary.ddml_late()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_ate.md).

- `oos_pred`:

  List of matrices, providing the reduced form predicted values.

- `learners`,`learners_DXZ`,`learners_ZX`,
  `cluster_variable`,`subsamples_Z0`,
  `subsamples_Z1`,`cv_subsamples_list_Z0`,
  `cv_subsamples_list_Z1`,`ensemble_type`:

  Pass-through of selected user-provided arguments. See above.

## Details

`ddml_late` provides a double/debiased machine learning estimator for
the local average treatment effect in the interactive model given by

\\Y = g_0(D, X) + U,\\

where \\(Y, D, X, Z, U)\\ is a random vector such that
\\\operatorname{supp} D = \operatorname{supp} Z = \\0,1\\\\, \\E\[U\vert
X, Z\] = 0\\, \\E\[Var(E\[D\vert X, Z\]\vert X)\] \neq 0\\,
\\\Pr(Z=1\vert X) \in (0, 1)\\ with probability 1, \\p_0(1, X) \geq
p_0(0, X)\\ with probability 1 where \\p_0(Z, X) \equiv \Pr(D=1\vert Z,
X)\\, and \\g_0\\ is an unknown nuisance function.

In this model, the local average treatment effect is defined as

\\\theta_0^{\textrm{LATE}} \equiv E\[g_0(1, X) - g_0(0, X)\vert p_0(1,
X) \> p(0, X)\]\\.

## References

Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). "Model Averaging
and Double Machine Learning." Journal of Applied Econometrics, 40(3):
249-269.

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). "Double/debiased machine learning for treatment and
structural parameters." The Econometrics Journal, 21(1), C1-C68.

Imbens G, Angrist J (1004). "Identification and Estimation of Local
Average Treatment Effects." Econometrica, 62(2), 467-475.

Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2),
241-259.

## See also

[`summary.ddml_late()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_ate.md)

Other ddml:
[`ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/ddml_ate.md),
[`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md),
[`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md)

## Examples

``` r
# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
Z = AE98[, "samesex"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]

# Estimate the local average treatment effect using a single base learner,
#     ridge.
late_fit <- ddml_late(y, D, Z, X,
                      learners = list(what = mdl_glmnet,
                                      args = list(alpha = 0)),
                      sample_folds = 2,
                      silent = TRUE)
summary(late_fit)
#> LATE estimation results: 
#>  
#>   Estimate Std. Error t value Pr(>|t|)
#>     -0.223      0.185    -1.2    0.229

# Estimate the local average treatment effect using short-stacking with base
#     learners ols, lasso, and ridge. We can also use custom_ensemble_weights
#     to estimate the ATE using every individual base learner.
weights_everylearner <- diag(1, 3)
colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
late_fit <- ddml_late(y, D, Z, X,
                      learners = list(list(fun = ols),
                                      list(fun = mdl_glmnet),
                                      list(fun = mdl_glmnet,
                                           args = list(alpha = 0))),
                      ensemble_type = 'nnls',
                      custom_ensemble_weights = weights_everylearner,
                      shortstack = TRUE,
                      sample_folds = 2,
                      silent = TRUE)
summary(late_fit)
#> LATE estimation results: 
#>  
#>           Estimate Std. Error t value Pr(>|t|)
#> nnls        -0.207      0.182   -1.14    0.255
#> mdl:ols     -0.183      0.178   -1.03    0.305
#> mdl:lasso   -0.187      0.179   -1.05    0.296
#> mdl:ridge   -0.200      0.181   -1.10    0.270
```
