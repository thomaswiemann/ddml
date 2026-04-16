# Plot Coefficients from a RAL Estimator

Plots point estimates with confidence intervals from an object
inheriting from class `ral`.

## Usage

``` r
# S3 method for class 'ral'
plot(
  x,
  parm = NULL,
  level = 0.95,
  uniform = TRUE,
  fit_idx = 1,
  type = "HC1",
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  col = "black",
  pch = 19,
  lwd = 1.5,
  ...
)

# S3 method for class 'ral_rep'
plot(
  x,
  parm = NULL,
  level = 0.95,
  uniform = TRUE,
  type = "HC1",
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  col = "black",
  pch = 19,
  lwd = 1.5,
  ...
)
```

## Arguments

- x:

  An object inheriting from class `ral`.

- parm:

  A specification of which parameters to plot. Either a vector of names
  or indices. Default: all.

- level:

  Numeric. Confidence level. Default `0.95`.

- uniform:

  Logical. If `TRUE`, uses uniform confidence bands via the multiplier
  bootstrap. Default `TRUE`.

- fit_idx:

  Integer. Which fit to plot (column index of `coefficients`). Default
  `1`.

- type:

  Character. HC type for standard errors. Default `"HC1"`.

- xlab:

  Character. Label for the x-axis.

- ylab:

  Character. Label for the y-axis.

- main:

  Character. Title for the plot.

- col:

  Color for points and segments. Default `"black"`.

- pch:

  Point character. Default `19` (solid dot).

- lwd:

  Line width for confidence interval segments. Default `1.5`.

- ...:

  Additional arguments passed to
  [`plot.default`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns a list with components `coefficients`, `ci`, and
`labels`.

## See also

[`confint.ral`](https://www.thomaswiemann.com/ddml/reference/confint.ral.md),
[`summary.ral`](https://www.thomaswiemann.com/ddml/reference/summary.ral.md)

## Examples

``` r
# Simulate a simple example
n <- 200
X <- cbind(1, stats::rnorm(n))
theta <- c(0.5, -0.3)
inf <- matrix(stats::rnorm(n * 2), n, 2)
obj <- ral(matrix(theta, 2, 1),
           array(inf, c(n, 2, 1)),
           nobs = n,
           coef_names = c("b1", "b2"))
plot(obj)

```
