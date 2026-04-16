# Construct a `ddml` Object.

Build a `"ddml"` object from user-supplied score components. The
resulting object inherits all S3 methods available for `ddml` objects,
including
[`summary.ddml`](https://www.thomaswiemann.com/ddml/reference/summary.ddml.md),
[`confint.ral`](https://www.thomaswiemann.com/ddml/reference/confint.ral.md),
[`vcov.ral`](https://www.thomaswiemann.com/ddml/reference/vcov.ral.md),
and
[`tidy.ddml`](https://www.thomaswiemann.com/ddml/reference/tidy.ddml.md).

## Usage

``` r
ddml(
  coefficients,
  scores,
  J,
  inf_func,
  nobs,
  coef_names,
  estimator_name,
  ensemble_type = colnames(coefficients),
  cluster_variable = seq_len(nobs),
  sample_folds = NULL,
  cv_folds = NULL,
  shortstack = FALSE,
  ensemble_weights = NULL,
  mspe = NULL,
  r2 = NULL,
  fitted = NULL,
  splits = NULL,
  call = match.call(),
  subclass = NULL,
  dinf_dtheta = NULL,
  ...
)
```

## Arguments

- coefficients:

  A `(p x nensb)` matrix of estimated target parameters. Rows correspond
  to components of \\\theta\\, columns to ensemble types.

- scores:

  A 3D array of evaluated Neyman orthogonal scores with dimensions
  `(n x p x nensb)`.

- J:

  A 3D array of evaluated Jacobians with dimensions `(p x p x nensb)`.

- inf_func:

  A 3D array of evaluated influence functions with dimensions
  `(n x p x nensb)`.

- nobs:

  Number of observations.

- coef_names:

  Character vector of coefficient names (length `p`).

- estimator_name:

  Character string identifying the estimator (e.g.,
  `"My Custom Estimator"`).

- ensemble_type:

  Character vector of ensemble types. Defaults to
  `colnames(coefficients)`.

- cluster_variable:

  A vector of cluster indices. Defaults to `seq_len(nobs)`.

- sample_folds:

  Number of cross-fitting folds used. Optional.

- cv_folds:

  Number of cross-validation folds used. Optional.

- shortstack:

  Logical indicating whether short-stacking was used. Default `FALSE`.

- ensemble_weights:

  A named list of ensemble weight matrices. Optional.

- mspe:

  A named list of per-learner MSPEs. Optional.

- r2:

  A named list of per-learner R-squared values. Optional.

- fitted:

  A named list of per-equation cross-fitted prediction objects.
  Optional.

- splits:

  A list of sample split objects. Optional.

- call:

  The matched call. Defaults to
  [`match.call()`](https://rdrr.io/r/base/match.call.html).

- subclass:

  Optional character string for a subclass name. If provided, the object
  will have class `c(subclass, "ddml")`.

- dinf_dtheta:

  An optional 4D array of dimensions `(nobs x p x p x nensb)` containing
  the derivatives of the influence functions.

- ...:

  Additional named components to include in the object.

## Value

An object of S3 class `"ddml"` (or `c(subclass, "ddml")` if `subclass`
is specified). See
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)
for the output structure.

## See also

Other utilities:
[`crosspred()`](https://www.thomaswiemann.com/ddml/reference/crosspred.md),
[`crossval()`](https://www.thomaswiemann.com/ddml/reference/crossval.md),
[`diagnostics()`](https://www.thomaswiemann.com/ddml/reference/diagnostics.md),
[`ensemble()`](https://www.thomaswiemann.com/ddml/reference/ensemble.md),
[`ensemble_weights()`](https://www.thomaswiemann.com/ddml/reference/ensemble_weights.md),
[`shortstacking()`](https://www.thomaswiemann.com/ddml/reference/shortstacking.md)

## Examples

``` r
# A minimal example: construct a ddml object from pre-computed
#     score components for a simple mean estimator.
n <- 100
y <- rnorm(n)
theta <- mean(y)

scores <- array(y - theta, dim = c(n, 1, 1))
J <- array(-1, dim = c(1, 1, 1))
psi_b <- list(matrix(y, ncol = 1))
psi_a <- list(array(-1, dim = c(n, 1, 1)))
inf_func <- array(y - theta, dim = c(n, 1, 1))
dinf_dtheta <- array(1, dim = c(n, 1, 1, 1))
coef <- matrix(theta, 1, 1, dimnames = list("mean", "custom"))

fit <- ddml(coefficients = coef, scores = scores, J = J,
        inf_func = inf_func, nobs = n, coef_names = "mean",
        dinf_dtheta = dinf_dtheta,
        estimator_name = "Sample Mean")
summary(fit)
#> DDML estimation: Sample Mean 
#> Obs: 100   Folds:
#> 
#>      Estimate Std. Error z value Pr(>|z|)  
#> mean  -0.2161     0.0916   -2.36    0.018 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
