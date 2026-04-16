# Variance-Covariance Matrix for RAL Estimators

Computes a heteroskedasticity-robust variance-covariance matrix.

## Usage

``` r
# S3 method for class 'ral'
vcov(object, fit_idx = 1, type = "HC1", ...)
```

## Arguments

- object:

  An object inheriting from class `ral`.

- fit_idx:

  Integer index of the fit. Defaults to 1.

- type:

  Character. One of `"HC1"` (default), `"HC0"`, or `"HC3"`.

- ...:

  Currently unused.

## Value

A \\p \times p\\ variance-covariance matrix.

## Details

Let \\\hat\phi_i\\ denote the estimated influence function at
observation \\i\\. Three variance estimators are available:

**HC0**: \$\$V\_{\mathrm{HC0}} = \frac{1}{n^2}\sum_i
\hat\phi_i\\\hat\phi_i'\$\$

**HC1** (default): \$\$V\_{\mathrm{HC1}} = V\_{\mathrm{HC0}} \times
\frac{n}{n - p}\$\$

**HC3**: \$\$V\_{\mathrm{HC3}} = \frac{1}{n^2}\sum_i
\frac{\hat\phi_i\\\hat\phi_i'} {(1 - \hat{h}\_{\theta,i})^2}\$\$

where \\\hat{h}\_{\theta,i}\\ is the leverage; see
[`hatvalues.ral`](https://www.thomaswiemann.com/ddml/reference/hatvalues.ral.md).

**Cluster-robust inference.** When `cluster_variable` is non-`NULL` and
identifies fewer groups than observations, the observation-level
influence functions are aggregated to cluster-level influence functions

\$\$\hat\Phi_g = \frac{G}{n} \sum\_{i \in C_g} \hat\phi_i\$\$

and the variance is computed as

\$\$V\_{\mathrm{HC0}} = \frac{1}{G^2} \sum\_{g=1}^{G}
\hat\Phi_g\\\hat\Phi_g'.\$\$

## See also

[`hatvalues.ral`](https://www.thomaswiemann.com/ddml/reference/hatvalues.ral.md),
[`confint.ral`](https://www.thomaswiemann.com/ddml/reference/confint.ral.md)
