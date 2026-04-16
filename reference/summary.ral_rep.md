# Summary for RAL Rep Objects

Aggregates coefficient estimates and covariance matrices across
independent replications.

## Usage

``` r
# S3 method for class 'ral_rep'
summary(
  object,
  aggregation = c("median", "mean", "spectral"),
  type = "HC1",
  ...
)

# S3 method for class 'summary.ral_rep'
print(x, digits = 3, ...)
```

## Arguments

- object:

  An object inheriting from class `ral_rep`.

- aggregation:

  Character string. Aggregation rule.

- type:

  Character. HC type. Default `"HC1"`.

- ...:

  Currently unused.

- x:

  An object of class `summary.ral_rep`.

- digits:

  Number of significant digits. Default 3.

## Value

An object of class `"summary.ral_rep"`.

## Details

Let \\\hat\theta_s\\ and \\\hat\Sigma_s\\ denote the coefficient vector
and sandwich covariance matrix from replication \\s\\.

**Coefficient aggregation.** For `"mean"`: \\\tilde\theta = S^{-1}
\sum\_{s=1}^{S} \hat\theta_s\\. For `"median"` and `"spectral"`:
\\\tilde\theta_j = \mathrm{median}\_{s}(\hat\theta\_{s,j})\\.

**Covariance aggregation.** Define the inflated per-replication
covariance as \$\$V_s = \hat\Sigma_s + (\hat\theta_s - \tilde\theta)
(\hat\theta_s - \tilde\theta)^\top.\$\$ For `"mean"`: \\\tilde\Sigma =
S^{-1} \sum\_{s=1}^{S} V_s\\. For `"median"`: \\\tilde\Sigma\_{s,ij} =
\mathrm{median}\_{s}(V\_{s,ij})\\. For `"spectral"`: solved via CVXR,
guaranteeing PSD.

## References

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). "Double/debiased machine learning for treatment and
structural parameters." The Econometrics Journal, 21(1), C1-C68.
