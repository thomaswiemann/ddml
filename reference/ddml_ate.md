# Estimators of Average Treatment Effects.

Estimators of the average treatment effect and the average treatment
effect on the treated.

## Usage

``` r
ddml_ate(
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
  subsamples_byD = NULL,
  cv_subsamples_byD = NULL,
  trim = 0.01,
  silent = FALSE
)

ddml_att(
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
  subsamples_byD = NULL,
  cv_subsamples_byD = NULL,
  trim = 0.01,
  silent = FALSE
)
```

## Arguments

- y:

  The outcome variable.

- D:

  The binary endogenous variable of interest.

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

- subsamples_byD:

  List of two lists corresponding to the two treatment levels. Each list
  contains vectors with sample indices for cross-fitting.

- cv_subsamples_byD:

  List of two lists, each corresponding to one of the two treatment
  levels. Each of the two lists contains lists, each corresponding to a
  subsample and contains vectors with subsample indices for
  cross-validation.

- trim:

  Number in (0, 1) for trimming the estimated propensity scores at
  `trim` and `1-trim`.

- silent:

  Boolean to silence estimation updates.

## Value

`ddml_ate` and `ddml_att` return an object of S3 class `ddml_ate` and
`ddml_att`, respectively. An object of class `ddml_ate` or `ddml_att` is
a list containing the following components:

- `ate` / `att`:

  A vector with the average treatment effect / average treatment effect
  on the treated estimates.

- `weights`:

  A list of matrices, providing the weight assigned to each base learner
  (in chronological order) by the ensemble procedure.

- `mspe`:

  A list of matrices, providing the MSPE of each base learner (in
  chronological order) computed by the cross-validation step in the
  ensemble construction.

- `psi_a`, `psi_b`:

  Matrices needed for the computation of scores. Used in
  [`summary.ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_ate.md)
  or
  [`summary.ddml_att()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_ate.md).

- `oos_pred`:

  List of matrices, providing the reduced form predicted values.

- `learners`,`learners_DX`,`cluster_variable`,
  `subsamples_D0`,`subsamples_D1`,
  `cv_subsamples_list_D0`,`cv_subsamples_list_D1`, `ensemble_type`:

  Pass-through of selected user-provided arguments. See above.

## Details

`ddml_ate` and `ddml_att` provide double/debiased machine learning
estimators for the average treatment effect and the average treatment
effect on the treated, respectively, in the interactive model given by

\\Y = g_0(D, X) + U,\\

where \\(Y, D, X, U)\\ is a random vector such that
\\\operatorname{supp} D = \\0,1\\\\, \\E\[U\vert D, X\] = 0\\, and
\\\Pr(D=1\vert X) \in (0, 1)\\ with probability 1, and \\g_0\\ is an
unknown nuisance function.

In this model, the average treatment effect is defined as

\\\theta_0^{\textrm{ATE}} \equiv E\[g_0(1, X) - g_0(0, X)\]\\.

and the average treatment effect on the treated is defined as

\\\theta_0^{\textrm{ATT}} \equiv E\[g_0(1, X) - g_0(0, X)\vert D =
1\]\\.

## References

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). "Double/debiased machine learning for treatment and
structural parameters." The Econometrics Journal, 21(1), C1-C68.

Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2),
241-259.

## See also

[`summary.ddml_ate()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_ate.md),
[`summary.ddml_att()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml_ate.md)

Other ddml:
[`ddml_fpliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_fpliv.md),
[`ddml_late()`](https://www.thomaswiemann.com/ddml/reference/ddml_late.md),
[`ddml_pliv()`](https://www.thomaswiemann.com/ddml/reference/ddml_pliv.md),
[`ddml_plm()`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md)

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
#>     -0.141     0.0152   -9.28 1.65e-20

# Estimate the average treatment effect using short-stacking with base
#     learners ols, lasso, and ridge. We can also use custom_ensemble_weights
#     to estimate the ATE using every individual base learner.
weights_everylearner <- diag(1, 3)
colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
ate_fit <- ddml_ate(y, D, X,
                    learners = list(list(fun = ols),
                                    list(fun = mdl_glmnet),
                                    list(fun = mdl_glmnet,
                                         args = list(alpha = 0))),
                    ensemble_type = 'nnls',
                    custom_ensemble_weights = weights_everylearner,
                    shortstack = TRUE,
                    sample_folds = 2,
                    silent = TRUE)
#> Warning: nnls: 1 propensity scores were trimmed.
summary(ate_fit)
#> ATE estimation results: 
#>  
#>           Estimate Std. Error t value Pr(>|t|)
#> nnls        -0.143     0.0153   -9.36 7.96e-21
#> mdl:ols     -0.143     0.0155   -9.22 2.92e-20
#> mdl:lasso   -0.143     0.0153   -9.34 1.01e-20
#> mdl:ridge   -0.143     0.0153   -9.39 6.16e-21
```
