# Linear Combinations of DDML Coefficients

Computes linear combinations \\R'\hat\theta\\ of DDML coefficient
estimates.

## Usage

``` r
lincom(fit, R, ...)

# S3 method for class 'ddml'
lincom(
  fit,
  R,
  fit_idx = NULL,
  labels = NULL,
  inf_func_R = NULL,
  dinf_dR = NULL,
  ...
)

# S3 method for class 'ddml_rep'
lincom(
  fit,
  R,
  inf_func_R = NULL,
  dinf_dR = NULL,
  fit_idx = NULL,
  labels = NULL,
  ...
)

# S3 method for class 'lincom'
print(x, ...)

# S3 method for class 'lincom_rep'
print(x, ...)
```

## Arguments

- fit:

  A `ddml` or `ddml_rep` object.

- R:

  A \\(p \times q)\\ contrast matrix. Each column defines one linear
  combination.

- ...:

  Currently unused.

- fit_idx:

  Integer index of the fit to use, or `NULL` (default) for all ensemble
  types. When `NULL`, the output carries all ensembles from the parent
  fit.

- labels:

  Optional character vector of length \\q\\ naming the linear
  combinations. Defaults to column names of `R`, or `"lc1"`, `"lc2"`,
  etc.

- inf_func_R:

  An optional \\(n \times p \times q)\\ n x p x q array of influence
  functions \\\Phi\_{R,i}\\ for the contrast matrix \\R\\. Slice `[,,k]`
  contains the IFs for column \\k\\ of \\R\\. When supplied, a
  delta-method correction is applied to the variance. When `NULL`
  (default), \\R\\ is treated as fixed and only the first term
  contributes.

- dinf_dR:

  An optional \\(n \times q \times q)\\ array of observation-level
  derivatives \\\partial \phi\_{R,i} / \partial R\\, representing the
  Weighting Leverage. When supplied, it is added to the Structural
  Leverage to form the total leverage used by HC3.

- x:

  A `lincom` or `lincom_rep` object.

## Value

An object of class `"lincom"` (inheriting from `"ral"`) for `ddml`
input, or `"lincom_rep"` (inheriting from `"ral_rep"`) for `ddml_rep`
input.

## Details

For a \\p\\-dimensional coefficient vector \\\hat\theta\\ and a \\(p
\times q)\\ contrast matrix \\R\\, the linear combination is \\\gamma =
R'\hat\theta\\.

The influence function for \\\gamma\\ is

\$\$\phi\_\gamma(W_i; \theta, R) = R'\\ \phi\_\theta(W_i; \theta) +
\Phi\_{R}(W_i)\\ \theta,\$\$

where \\\phi\_\theta\\ is the influence function of \\\hat\theta\\ (see
[`ral`](https://www.thomaswiemann.com/ddml/reference/ral.md) and
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md)),
\\\Phi_R(W_i)\\ is the \\(p \times q)\\ matrix of influence functions
for the contrast matrix \\R\\ (the \\i\\-th slice of `inf_func_R`), and
the second term vanishes when \\R\\ is fixed. The estimated influence
function is

\$\$\hat\phi\_{\gamma,i} = R'\\ \hat\phi\_{\theta,i} + \hat\Phi\_{R,i}\\
\hat\theta.\$\$

The leverage for \\\gamma\\ is

\$\$h\_\gamma(W_i; \theta, R) = \mathrm{tr}\\ \left( R'\\h\_\theta(W_i;
\theta)\\R + h_R(W_i)\right),\$\$

where \\h\_\theta\\ is the structural leverage from the parent estimator
(see
[`hatvalues.ral`](https://www.thomaswiemann.com/ddml/reference/hatvalues.ral.md))
and \\h_R\\ is the weighting leverage. The sample analog is

\$\$\hat{h}\_{\gamma,i} = \mathrm{tr}\\ \left(
R'\\\hat{h}\_{\theta,i}\\R + \hat{h}\_{R,i}\right),\$\$

where \\\hat{h}\_{\theta,i}\\ is mapped from the parent's `dinf_dtheta`
and \\\hat{h}\_{R,i}\\ from the optional `dinf_dR` argument.

The resulting `lincom` object inherits from `ral` and supports all
standard inference methods: `vcov`, `confint`, `summary`, `tidy`, and
`hatvalues`. For `ddml_rep` objects, `lincom` returns a `lincom_rep`
inheriting from `ral_rep`.

Note that `inf_func_R` is needed for inference when \\R\\ is estimated.
Leverage computation further requires `dinf_dR`. See
[`vcov.ral`](https://www.thomaswiemann.com/ddml/reference/vcov.ral.md)
and
[`hatvalues.ral`](https://www.thomaswiemann.com/ddml/reference/hatvalues.ral.md)
for more details.

## References

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). "Double/debiased machine learning for treatment and
structural parameters." The Econometrics Journal, 21(1), C1-C68.

## See also

[`lincom_weights_did`](https://www.thomaswiemann.com/ddml/reference/lincom_weights_did.md)
for constructing DiD aggregation weights.

Other ddml inference:
[`summary.ddml()`](https://www.thomaswiemann.com/ddml/reference/summary.ddml.md)

## Examples

``` r
# \donttest{
set.seed(42)
n <- 200; T_ <- 4
X <- matrix(rnorm(n * 2), n, 2)
G <- sample(c(3, 4, Inf), n, replace = TRUE,
            prob = c(0.3, 0.3, 0.4))
y <- matrix(rnorm(n * T_), n, T_)
for (i in seq_len(n)) {
  if (is.finite(G[i])) {
    for (j in seq_len(T_)) {
      if (j >= G[i]) y[i, j] <- y[i, j] + 1
    }
  }
}
fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                learners = list(what = ols),
                sample_folds = 2,
                silent = TRUE)
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
#> Warning: One of the crossfitting subsamples only uses 28 observations for training. Consider increasing ``sample_folds`` if possible.
# Simple contrast: first cell minus second
p <- nrow(fit$coefficients)
R <- matrix(0, p, 1)
R[1, 1] <- 1; R[2, 1] <- -1
lc <- lincom(fit, R = R, labels = "ATT1-ATT2")
summary(lc)
#> RAL estimation: Linear Combination 
#> Obs: 200
#> 
#>           Estimate Std. Error z value Pr(>|z|)    
#> ATT1-ATT2   -1.165      0.238    -4.9  9.6e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
