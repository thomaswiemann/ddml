# Construct a RAL Inference Object

Creates a regular asymptotically linear (RAL) inference object from
pre-computed influence functions. This is the base class for all
influence-function-based inference in ddml.

## Usage

``` r
ral(
  coefficients,
  inf_func,
  dinf_dtheta = NULL,
  nobs,
  coef_names,
  cluster_variable = NULL,
  estimator_name = "RAL estimator",
  subclass = NULL,
  ...
)
```

## Arguments

- coefficients:

  A \\p \times\\ `nfit` matrix of estimated coefficients. Rows are
  parameters, columns are fits (e.g., ensemble types).

- inf_func:

  A 3D array of dimension \\n \times p \times\\ `nfit`. The influence
  function evaluated at each observation.

- dinf_dtheta:

  Optional 4D array of dimension \\n \times p \times p \times\\ `nfit`.
  The derivative of the influence function with respect to \\\theta\\,
  used for HC3 leverage. If `NULL`, HC3 is unavailable.

- nobs:

  Integer number of observations.

- coef_names:

  Character vector of parameter names (length \\p\\).

- cluster_variable:

  Optional vector of cluster identifiers (length \\n\\). If non-`NULL`,
  cluster-robust inference is used.

- estimator_name:

  Character string for display.

- subclass:

  Optional character string prepended to the class vector.

- ...:

  Additional named elements stored in the object.

## Value

An object of class `ral` (or `c(subclass, "ral")`).

## Details

A regular asymptotically linear (RAL) estimator \\\hat\theta\\ satisfies

\$\$\hat\theta - \theta_0 = \frac{1}{n} \sum\_{i=1}^{n} \phi(W_i;
\theta_0) + o_p(n^{-1/2}),\$\$

where \\\phi(W_i; \theta_0)\\ is the *influence function*. This package
stores the estimated influence function \\\hat\phi_i \equiv \phi(W_i;
\hat\theta)\\ in the `inf_func` slot.

When an observation-level derivative \\\partial \hat\phi_i / \partial
\theta\\ is available (stored in `dinf_dtheta`), the estimator supports
HC3 inference via leverage; see
[`hatvalues.ral`](https://www.thomaswiemann.com/ddml/reference/hatvalues.ral.md).

The RAL framework is estimator-agnostic: it consumes pre-computed
influence functions and does not prescribe how they are obtained. For
the specific construction under cross-fitting and Neyman-orthogonal
scores, see
[`ddml-intro`](https://www.thomaswiemann.com/ddml/reference/ddml-intro.md).
For linear combinations of `ddml` estimators, see
[`lincom`](https://www.thomaswiemann.com/ddml/reference/lincom.md).
