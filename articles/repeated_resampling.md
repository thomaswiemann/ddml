# Robust Inference and Repeated Resampling

## Introduction

Double/Debiased Machine Learning revolves around the use of
**cross-fitting** (sample splitting) to estimate nuisance parameters and
structural effects without inducing severe over-fitting bias.

Because cross-fitting relies on random sample splits, structural
estimates inherently contain *splitting variance*. If you run a `ddml`
estimator twice with two different random seeds, you will likely get two
slightly different point estimates.

To stabilize these estimates and eliminate the influence of a “lucky” or
“unlucky” sample split, researchers recommend repeating the
cross-fitting split multiple times and aggregating the results (see
Chernozhukov et al. 2018 and Ahrens et al. 2024).

`ddml` introduces dedicated functionality for repeated resampling via
[`ddml_replicate()`](https://www.thomaswiemann.com/ddml/reference/ddml_replicate.md)
and the `ddml_rep` class.

``` r
library(ddml)
set.seed(847291)
```

## Single vs. Repeated Fits

Let’s use a subsample of the `AE98` dataset to show how a single fit
compares to a repeated fit.

``` r
# Use a subset of data
sub_idx = sample(1:nrow(AE98), 1000)
y = AE98[sub_idx, "worked"]
D = AE98[sub_idx, "morekids"]
X = AE98[sub_idx, c("age", "agefst", "black", "hisp", "othrace", "educ")]

# Define a simple learner list
learners <- list(
  list(what = ols),
  list(what = mdl_glmnet)
)
```

First, let’s look at the result of a *single* `ddml_plm` fit:

``` r
# Single fit
fit_single <- ddml_plm(
  y = y, D = D, X = X,
  learners = learners,
  shortstack = TRUE,
  sample_folds = 5,
  silent = TRUE
)

summary(fit_single)
#> DDML estimation: Partially Linear Model 
#> Obs: 1000   Folds: 5  Stacking: short-stack
#> 
#>             Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.11986    0.03306   -3.63  0.00029 ***
#> (Intercept)  0.00422    0.01541    0.27  0.78425    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Next, let’s use the
[`ddml_replicate()`](https://www.thomaswiemann.com/ddml/reference/ddml_replicate.md)
wrapper. We pass it the target estimator (`ddml_plm`) and the desired
number of independent resamples.

*Note: Repeated resampling takes proportionally longer to compute, as
the estimator is run entirely to completion `resamples` times in
sequence.*

``` r
# Replicated fit (e.g., 5 independent cross-fitting splits)
fit_rep <- ddml_replicate(
  ddml_plm,
  y = y, D = D, X = X,
  learners = learners,
  shortstack = TRUE,
  sample_folds = 5,
  resamples = 5,
  silent = TRUE
)

# Observe the class of the returned object
class(fit_rep)
#> [1] "ddml_rep" "ral_rep"
```

## Interacting with `ddml_rep` Objects

[`ddml_replicate()`](https://www.thomaswiemann.com/ddml/reference/ddml_replicate.md)
returns an object of class `ddml_rep`. This object natively supports all
the standard generics you would expect from a single `ddml` model.

When summary statistics or coefficients are requested, `ddml_rep`
automatically aggregates the results from all individual fits. By
default, it applies the **median** aggregation approach (Chernozhukov et
al., 2018), which computes the median of the parameter estimates across
all splits and suitably inflates the standard errors to account for
between-split variation.

``` r
# Print the object details
fit_rep
#> DDML replicated fits: Partially Linear Model 
#>   Resamples: 5   Obs: 1000   Folds: 5 
#> 
#> Use summary() for aggregated inference.
#> Use x[[i]] to access individual fits.

# Aggregate inference (Median Aggregation)
summary(fit_rep)
#> DDML estimation: Partially Linear Model 
#> Obs: 1000   Folds: 5   Resamples: 5   Aggregation: median  Stacking: short-stack
#> 
#>             Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.11409    0.03318   -3.44  0.00058 ***
#> (Intercept)  0.00384    0.01544    0.25  0.80372    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Extract only the coefficients
coef(fit_rep)
#>           D1  (Intercept) 
#> -0.114085768  0.003836949

# Extract confidence intervals (with HC3 standard errors passed to the VCoV estimator)
confint(fit_rep, type = "HC3")
#>                   2.5 %      97.5 %
#> D1          -0.17920978 -0.04896176
#> (Intercept) -0.02645093  0.03412483
#> attr(,"crit_val")
#> [1] 1.959964
```

You can alternatively switch to **mean** aggregation (Ahrens et al.,
2024), which computes the empirical mean across split estimates:

``` r
summary(fit_rep, aggregation = "mean")
#> DDML estimation: Partially Linear Model 
#> Obs: 1000   Folds: 5   Resamples: 5   Aggregation: mean  Stacking: short-stack
#> 
#>             Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.11553    0.03323   -3.48  0.00051 ***
#> (Intercept)  0.00415    0.01548    0.27  0.78843    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Visualizing Splitting Variance

We can inspect the exact variation across the five sample splits by
extracting the point estimates from the underlying list of fits in the
`ddml_rep` object:

``` r
# Individual point estimates for the primary ensemble across the 5 independent splits
split_estimates <- vapply(fit_rep$fits, function(fit) coef(fit)[1], numeric(1))
print(split_estimates)
#> [1] -0.1117597 -0.1187530 -0.1133071 -0.1197519 -0.1140858

# Variance of the point estimates across splits
var(split_estimates)
#> [1] 1.23637e-05
```

### Multiple Treatment Variables and Joint Variance

If your model includes multiple treatment variables, `ddml_rep`
seamlessly aggregates the full, joint variance-covariance matrix across
all the sample splits, preserving the covariance between different
parameter estimates.

``` r
# Fit a model with two treatments: 'morekids' and 'boy2nd'
fit_rep_multi <- ddml_replicate(
  ddml_plm,
  y = y, 
  D = AE98[sub_idx, c("morekids", "boy2nd")], 
  X = X,
  learners = learners,
  shortstack = TRUE,
  sample_folds = 5,
  resamples = 5,
  silent = TRUE
)

# Extract the full aggregated Variance-Covariance matrix
vcov(fit_rep_multi)
#>                  morekids        boy2nd   (Intercept)
#> morekids     1.090179e-03 -2.156755e-05 -7.584066e-06
#> boy2nd      -2.156755e-05  9.450183e-04  1.153763e-05
#> (Intercept) -7.584066e-06  1.153763e-05  2.368640e-04
```

### Extracting Individual Fits

Because a `ddml_rep` object is essentially a structured list of pure
`ddml` outputs, you can always extract the specific estimates and
diagnostics of any single resample using double brackets `[[ ]]`:

``` r
# Extract the 3rd replicate fit
run_3 <- fit_rep[[3]]
class(run_3)
#> [1] "ddml_plm" "ddml"     "ral"

# View its specific diagnostics
diagnostics(run_3)
#> Stacking diagnostics: Partially Linear Model 
#> Obs: 1000 
#> 
#>   y_X:
#>    learner   mspe     r2 weight_nnls
#>  learner_1 0.2402 0.0394      0.0000
#>  learner_2 0.2396 0.0415      0.9986
#>       nnls 0.2396 0.0415          NA
#> 
#>   D1_X:
#>    learner   mspe     r2 weight_nnls
#>  learner_1 0.2205 0.0659      0.0000
#>  learner_2 0.2202 0.0671      0.9926
#>       nnls 0.2202 0.0672          NA
#> 
#> Note: Ensemble MSPE and R2 for short-stacking rely on full-sample weights
#>        and represent in-sample fit over cross-fitted base predictions.
```

## Manual Construction and HPC Compatibility

[`ddml_replicate()`](https://www.thomaswiemann.com/ddml/reference/ddml_replicate.md)
computes sequential cross-fits using a simple `for` loop.

If you are running very large regressions (e.g., millions of
observations or complex machine learners), you might want to parallelize
these resamples across a high-performance computing (HPC) cluster.

Because `ddml` fits are completely independent of one another, you can
easily compute them in parallel (using `mclapply`, `future`, or job
arrays) and combine the raw list of returned objects using the
structural
[`ddml_rep()`](https://www.thomaswiemann.com/ddml/reference/ddml_rep.md)
constructor.

``` r
# 1. Generate an independent list of 3 ddml fits (this could be from parallel::mclapply)
fits_list <- lapply(1:3, function(r) {
  ddml_plm(
    y = y, D = D, X = X,
    learners = learners,
    shortstack = TRUE,
    sample_folds = 5,
    silent = TRUE
  )
})

# 2. Combine them into a ddml_rep object
manual_rep <- ddml_rep(fits_list)
class(manual_rep)
#> [1] "ddml_rep" "ral_rep"

# 3. Aggregate results
summary(manual_rep)
#> DDML estimation: Partially Linear Model 
#> Obs: 1000   Folds: 5   Resamples: 3   Aggregation: median  Stacking: short-stack
#> 
#>             Estimate Std. Error z value Pr(>|z|)    
#> D1          -0.11819    0.03308   -3.57  0.00035 ***
#> (Intercept)  0.00309    0.01543    0.20  0.84144    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

This flexibility allows `ddml` to scale efficiently from local laptops
to massive distributed compute environments, ensuring robust inference
no matter your hardware setup.

## References

Ahrens A, Chernozhukov V, Hansen C B, Kozbur D, Schaffer M E, Wiemann T
(2026). “An Introduction to Double/Debiased Machine Learning.” *Journal
of Economic Literature*, forthcoming.

Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). “Model Averaging
and Double Machine Learning.” *Journal of Applied Econometrics*, 40(3):
249-269.

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C, Newey W,
Robins J (2018). “Double/debiased machine learning for treatment and
structural parameters.” *The Econometrics Journal*, 21(1), C1-C68.
