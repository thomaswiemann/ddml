# Difference-in-Differences Aggregation Weights for lincom

Constructs the contrast matrix \\R\\ and its influence function matrix
`inf_func_R` for standard DiD aggregation types. The output is designed
to be passed directly to
[`lincom`](https://www.thomaswiemann.com/ddml/reference/lincom.md).

## Usage

``` r
lincom_weights_did(
  fit,
  type = c("dynamic", "group", "simple", "calendar"),
  min_e = -Inf,
  max_e = Inf,
  fit_idx = NULL
)
```

## Arguments

- fit:

  A `ddml_attgt` or `ddml_rep` object whose underlying fits are
  `ddml_attgt`.

- type:

  Aggregation type: `"dynamic"` (default), `"group"`, `"simple"`, or
  `"calendar"`.

- min_e, max_e:

  Event-time range filter (dynamic only). Cells with event time outside
  `[min_e, max_e]` are excluded.

- fit_idx:

  Integer index of the fit (ensemble type) to use for computing the
  weighting leverage, or `NULL` (default) for all ensemble types. When
  `NULL`, `dinf_dR` is a 4D array.

## Value

A list with elements:

- `R`:

  A \\(C \times q)\\ contrast matrix where \\C\\ is the number of GT
  cells.

- `inf_func_R`:

  An \\(n \times C \times q)\\ n x C x q array of influence functions
  for the contrast matrix \\R\\. Slice `[,,k]` contains the IFs for
  column \\k\\ of \\R\\.

- `dinf_dR`:

  An \\(n \times q \times q)\\ n x q x q array of weighting leverage.
  Each slice is the constant matrix \\V'V\\ where \\V\_{g,k} = \sum\_{c
  \in \mathcal{K}\_k,\\ g_c = g} (\theta_c - \gamma_k) / S_k\\. See D89
  §5.3.

- `labels`:

  Character vector of length \\q\\ naming the aggregated quantities.

## Details

Let \\\theta^{(g,t)}\_0\\ denote the GT-ATT from
[`ddml_attgt`](https://www.thomaswiemann.com/ddml/reference/ddml_attgt.md).
Each aggregation type defines a summary parameter as a weighted average
of GT-ATTs over a subset of post-treatment cells (\\t \geq g\\).

**Dynamic** (`type = "dynamic"`): aggregates by event time \\e = t - g\\
(Callaway and Sant'Anna, 2021, eq. 9). For each \\e\\:

\$\$\tau_0(e) = \sum\_{g \in \mathcal{G}} \mathbf{1}\\g + e \leq T\\\\
\Pr(G = g \mid G + e \leq T)\\ \theta^{(g,\\ g+e)}\_0.\$\$

**Group** (`type = "group"`): aggregates by cohort \\g\\. For each
\\g\\:

\$\$\theta(g) = \frac{1}{\|\mathcal{T}\_g\|} \sum\_{t \in
\mathcal{T}\_g} \theta^{(g,t)}\_0,\$\$

where \\\mathcal{T}\_g = \\t : t \geq g\\\\ and the weights reduce to
uniform within each cohort.

**Calendar** (`type = "calendar"`): aggregates by time period \\t\\. For
each \\t\\:

\$\$\theta(t) = \sum\_{g:\\ g \leq t} \frac{P(G = g)}{\sum\_{g':\\ g'
\leq t} P(G = g')}\\ \theta^{(g,t)}\_0.\$\$

**Simple** (`type = "simple"`): a single weighted average across all
post-treatment cells:

\$\$\theta\_{ATT} = \sum\_{(g,t):\\ t \geq g} \frac{P(G =
g)}{\sum\_{(g',t'):\\ t' \geq g'} P(G = g')}\\ \theta^{(g,t)}\_0.\$\$

The influence function for the estimated weights is derived via the
quotient rule and passed to
[`lincom`](https://www.thomaswiemann.com/ddml/reference/lincom.md) as
`inf_func_R`.

## References

Callaway B, Sant'Anna P H C (2021). "Difference-in-Differences with
multiple time periods." Journal of Econometrics, 225(2), 200-230.

## See also

[`lincom`](https://www.thomaswiemann.com/ddml/reference/lincom.md)

## Examples

``` r
# \donttest{
set.seed(42)
n <- 200; T_ <- 4
X <- matrix(rnorm(n * 2), n, 2)
G <- sample(c(3, 4, Inf), n, replace = TRUE,
            prob = c(0.3, 0.3, 0.4))
y <- matrix(rnorm(n * T_), n, T_)
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
w <- lincom_weights_did(fit, type = "dynamic")
dyn <- lincom(fit, R = w$R,
              inf_func_R = w$inf_func_R,
              dinf_dR = w$dinf_dR,
              labels = w$labels)
summary(dyn)
#> RAL estimation: Linear Combination 
#> Obs: 200
#> 
#>      Estimate Std. Error z value Pr(>|z|)
#> e=-3  -0.0144     0.3091   -0.05     0.96
#> e=-2  -0.1743     0.1809   -0.96     0.34
#> e=0    0.0443     0.1879    0.24     0.81
#> e=1    0.0890     0.2664    0.33     0.74
# }
```
