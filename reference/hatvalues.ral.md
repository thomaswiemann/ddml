# Extract leverage (Hat Values)

Computes the leverage (hat values) for a RAL estimator. Used internally
for HC3 standard errors.

## Usage

``` r
# S3 method for class 'ral'
hatvalues(model, fit_idx = 1, ...)
```

## Arguments

- model:

  An object inheriting from class `ral`.

- fit_idx:

  Integer index of the fit to extract leverage values for. Defaults to
  1.

- ...:

  Currently unused.

## Value

A numeric vector of leverage values.

## Details

The leverage for observation \\i\\ is

\$\$h_i(\theta) = -\mathrm{tr}\\\left( \frac{1}{n} \frac{\partial
\phi_i(\theta)} {\partial \theta} \right).\$\$

The sample analog replaces \\\phi_i(\theta)\\ with its estimate
\\\hat\phi_i\\: \\\hat{h}\_i = h_i(\hat\theta).\\

The derivative \\\partial \hat\phi_i / \partial \theta\\ is stored in
the `dinf_dtheta` slot. For the specific form of this derivative in the
DML context, see
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md).
For the leverage of linear combinations, see
[`lincom`](https://www.thomaswiemann.com/ddml/reference/lincom.md).
