# Neyman Orthogonality in Linear Regression

## Introduction

This article builds intuition for *Neyman orthogonality*—a key component
of Double/Debiased Machine Learning (DML). Throughout, we focus on
ordinary least squares.

Consider estimating the coefficient $\theta_{0}$ on a variable of
interest $D$ in a linear model with controls $X$. Introductory
econometrics courses typically discuss *partialling out* $X$ before
estimating $\theta_{0}$. Two variants of this idea lead to the same
point estimate:

- **Procedure A** (partial $D$ only): residualize $D$ with respect to
  $X$, then regress $Y$ on $\widetilde{D}$.
- **Procedure B** (partial both): residualize *both* $Y$ and $D$ with
  respect to $X$, then regress $\widetilde{Y}$ on $\widetilde{D}$.

It is tempting to take the first-step coefficients from partialling out
$X$ from $D$ and $Y$ as given, and compute standard errors only for the
second-step regression coefficient. As it turns out, this naive approach
results in valid inference in Procedure B, but is incorrect for
Procedure A. Here, we aim to explain why one procedure is immune to
first-step estimation error while the other is not—a property
conveniently summarized by the notion of *Neyman orthogonality*.

The article covers:

1.  Setting up the two procedures and comparing their output to full OLS
2.  Deriving how first-step estimation error propagates—or doesn’t—into
    inference about $\theta_{0}$
3.  Introducing Neyman orthogonality and its implications
4.  Discussing the bridge from OLS to nonparametric/ML first steps and
    `ddml_plm`

## Two-Step Procedures with Identical Coefficient Estimates

To fix ideas, we generate data from a linear model with a single
control:

$$Y_{i} = \theta_{0}D_{i} + \beta_{0}X_{i} + U_{i},$$

where $D_{i} = \alpha X_{i} + V_{i}$ so that $X$ is correlated with both
$D$ and $Y$.

``` r
library(ddml)
set.seed(48105)

n <- 500
theta_0 <- 1
beta_0 <- 2
alpha <- 0.8

X <- rnorm(n)
V <- rnorm(n)
D <- alpha * X + V
U <- rnorm(n)
Y <- theta_0 * D + beta_0 * X + U
```

### Procedure A: Residualize $D$ Only

For our first procedure, we partial out $X$ from $D$, then regress $Y$
on the residual:

``` r
D_tilde <- residuals(lm(D ~ X))
fit_A <- lm(Y ~ D_tilde)
```

### Procedure B: Residualize Both $Y$ and $D$

For our second procedure, we partial out $X$ from *both* $D$ and $Y$,
then regress the residuals on each other:

``` r
Y_tilde <- residuals(lm(Y ~ X))
fit_B <- lm(Y_tilde ~ D_tilde - 1)
```

### Comparing to Full OLS

It is useful to compare both two-step procedures with conventional
multiple regression.

``` r
library(sandwich)
library(lmtest)

fit_full <- lm(Y ~ D + X)

theta_A <- coef(fit_A)["D_tilde"]
theta_B <- coef(fit_B)["D_tilde"]
theta_full <- coef(fit_full)["D"]

se_A <- sqrt(diag(vcovHC(fit_A, type = "HC0")))["D_tilde"]
se_B <- sqrt(diag(vcovHC(fit_B, type = "HC0")))["D_tilde"]
se_full <- sqrt(diag(vcovHC(fit_full, type = "HC0")))["D"]

cat("Full OLS:      lm(Y ~ D + X)\n")
#> Full OLS:      lm(Y ~ D + X)
cat("  theta_hat:", round(theta_full, 4),
    "  SE:", round(se_full, 4), "\n\n")
#>   theta_hat: 0.9988   SE: 0.047

cat("Procedure A:   lm(Y ~ D_tilde)\n")
#> Procedure A:   lm(Y ~ D_tilde)
cat("  theta_hat:", round(theta_A, 4),
    "  SE:", round(se_A, 4), "\n\n")
#>   theta_hat: 0.9988   SE: 0.1188

cat("Procedure B:   lm(Y_tilde ~ D_tilde - 1)\n")
#> Procedure B:   lm(Y_tilde ~ D_tilde - 1)
cat("  theta_hat:", round(theta_B, 4),
    "  SE:", round(se_B, 4), "\n")
#>   theta_hat: 0.9988   SE: 0.047
```

All three point estimates are identical. This is no accident:
$\widetilde{D}$ is orthogonal to $X$ by construction, so omitting $X$
does not bias the coefficient on $\widetilde{D}$.

But the standard errors tell a different story. **Procedure B** matches
full OLS exactly—both point estimate and standard error. **Procedure A**
gives a standard error that is far too large.

In both procedures, the software computes the sandwich standard error
treating first-step coefficients as fixed. **If both procedures ignore
first-step estimation uncertainty, why does Procedure A fail while
Procedure B succeeds?**

## The First-Order Impact of Nuisance Estimation

Estimation of $\theta_{0}$ in each procedure is based on solving a
sample average of a score function
$m\left( W_{i};\theta,\eta \right) = 0$, where $\eta$ denotes nuisance
parameters.

For **Procedure A**, the score is

$$m_{A}\left( W_{i};\theta,\eta_{D} \right) = {\widetilde{D}}_{i}\left( Y_{i} - {\widetilde{D}}_{i}\theta \right),$$

where ${\widetilde{D}}_{i} = D_{i} - X_{i}\prime\eta_{D}$ and the
nuisance parameter is $\eta = \eta_{D}$ (the coefficient from the
regression of $D$ on $X$).

For **Procedure B**, the score is

$$m_{B}\left( W_{i};\theta,\eta \right) = {\widetilde{D}}_{i}\left( {\widetilde{Y}}_{i} - {\widetilde{D}}_{i}\theta \right),$$

where ${\widetilde{Y}}_{i} = Y_{i} - X_{i}\prime\eta_{Y}$ and the
nuisance parameter is $\eta = \left( \eta_{Y},\eta_{D} \right)$ (the
coefficients from the regressions of $Y$ and $D$ on $X$, respectively).

Given the score, software computes the sandwich standard error by
treating $\widehat{\eta}$ as known:

$$\widehat{\text{SE}}\left( \widehat{\theta} \right) = \sqrt{{\widehat{J}}^{- 2} \cdot \frac{1}{n^{2}}\sum\limits_{i = 1}^{n}m\left( W_{i};\,\widehat{\theta},\,\widehat{\eta} \right)^{2}}\,,$$

where
$\widehat{J} = \frac{1}{n}\sum_{i}\frac{\partial}{\partial\theta}m\left( W_{i};\,\widehat{\theta},\,\widehat{\eta} \right)$
is the sample Jacobian. Both procedures share the same Jacobian
$\widehat{J} \approx - E\left\lbrack {\widetilde{D}}^{2} \right\rbrack$,
so the SE difference traces entirely to the squared scores in the
numerator. At $\theta_{0}$, Procedure A’s score residual
$Y_{i} - {\widetilde{D}}_{i}\theta_{0}$ still contains terms involving
$X_{i}$, inflating the numerator. Procedure B’s score residual
${\widetilde{Y}}_{i} - {\widetilde{D}}_{i}\theta_{0} = U_{i}$ contains
only the structural error—hence its smaller SE.

Whether the asymptotic distribution of $\widehat{\theta}$ depends on
$\widehat{\eta}$ can be assessed via a standard Taylor expansion of the
sample moment around the true parameters
$\left( \theta_{0},\eta_{0} \right)$ (see, e.g., Ahrens et al., 2026):

$$\sqrt{n}\left( \widehat{\theta} - \theta_{0} \right) = - J_{\theta}^{- 1}\lbrack\underset{\text{CLT}}{\underbrace{\frac{1}{\sqrt{n}}\sum\limits_{i}m\left( W_{i};\theta_{0},\eta_{0} \right)}} + \underset{{( \star )}:{\mspace{6mu}\text{first-order impact of nuisance estimation}}}{\underbrace{\frac{1}{\sqrt{n}}\sum\limits_{i}\frac{\partial}{\partial\eta}m\left( W_{i};\theta_{0},\eta_{0} \right)\left( \widehat{\eta} - \eta_{0} \right)}}\rbrack$$

plus higher order terms. The term $( \star )$ captures the first-order
impact of estimating the nuisance parameter $\eta_{0}$. If $( \star )$
vanishes, inference about $\theta_{0}$ can proceed as if $\eta_{0}$ were
known.

### Procedure A: $( \star )$ Does Not Vanish

Procedure A’s score depends on $\eta_{D}$ through
${\widetilde{D}}_{i} = D_{i} - X_{i}\prime\eta_{D}$. The pathwise
derivative of the expected score with respect to the nuisance parameter
is:

$$\frac{\partial}{\partial\eta_{D}\prime}\, E\left\lbrack m_{A}\left( W;\,\theta_{0},\eta_{D} \right) \right\rbrack|_{\eta_{D,0}} = - E\lbrack XX\prime\rbrack\,\beta_{Y|X}$$

where $\beta_{Y|X} = \theta_{0}\gamma_{D|X} + \beta_{0}$ is the
coefficient from the regression of $Y$ on $X$. This is **non-zero**
whenever $X$ is related to $Y$.

As a consequence, estimation error in ${\widehat{\eta}}_{D}$ propagates
directly into ${\widehat{\theta}}_{A}$. The term $( \star )$ does not
vanish, and the naive standard error—which ignores $( \star )$—is wrong.

### Procedure B: $( \star )$ Vanishes

Procedure B’s score depends on two nuisance parameters: $\eta_{Y}$ (from
$Y$ on $X$) and $\eta_{D}$ (from $D$ on $X$). The pathwise derivatives
of the expected score are:

$$\frac{\partial}{\partial\eta_{Y}\prime}E\left\lbrack m_{B} \right\rbrack = - E\left\lbrack \widetilde{D}\, X\prime \right\rbrack = 0$$

because $\widetilde{D}\bot X$ by the projection of $D$ on $X$, and

$$\frac{\partial}{\partial\eta_{D}\prime}E\left\lbrack m_{B} \right\rbrack = - E\left\lbrack X\widetilde{Y} \right\rbrack + 2\theta_{0}E\left\lbrack X\widetilde{D} \right\rbrack = - E\lbrack XU\rbrack + \theta_{0}E\left\lbrack X\widetilde{D} \right\rbrack = 0$$

because $\widetilde{Y} = \widetilde{D}\theta_{0} + U$ and both $U$ and
$\widetilde{D}$ are orthogonal to $X$.

**Both pathwise derivatives vanish.** The term $( \star )$ in the
expansion is zero, so first-step estimation of ${\widehat{\eta}}_{Y}$
and ${\widehat{\eta}}_{D}$ does not affect the asymptotic distribution
of ${\widehat{\theta}}_{B}$—even to first order. The naive standard
error is correct as-is.

### Numerical Verification

We can verify the pathwise derivatives numerically:

``` r
J_gamma_A <- mean(X * (Y - D_tilde * theta_0))
J_gamma_B <- mean(X * (Y_tilde - D_tilde * theta_0))

cat("Pathwise derivatives (should vanish for B only):\n")
#> Pathwise derivatives (should vanish for B only):
cat("  Procedure A: E[X*(Y  - Dt*theta)] =",
    round(J_gamma_A, 4), "\n")
#>   Procedure A: E[X*(Y  - Dt*theta)] = 2.672
cat("  Procedure B: E[X*(Yt - Dt*theta)] =",
    round(J_gamma_B, 4), "\n")
#>   Procedure B: E[X*(Yt - Dt*theta)] = 0
```

## Neyman Orthogonality

The property that makes Procedure B’s standard error correct is called
**Neyman orthogonality**. Formally, a score $m(W;\,\theta,\eta)$ is
Neyman orthogonal at $\left( \theta_{0},\eta_{0} \right)$ if small
perturbations of the nuisance parameter away from $\eta_{0}$ do not
create first-order changes in the expected score:

$$\left. \frac{\partial}{\partial\lambda}E\lbrack m\left( W;\,\theta_{0},\,\eta_{0} + \lambda\left( \eta - \eta_{0} \right) \right)\rbrack \right|_{\lambda = 0} = 0,\quad\forall\,\eta.$$

Procedure B’s score satisfies this condition: the expected score is
insensitive to perturbations in *both* $\eta_{Y}$ and $\eta_{D}$. As a
result, replacing the true nuisance parameters with first-step estimates
does not distort inference about $\theta_{0}$.

Procedure A’s score does not satisfy Neyman orthogonality: it is
sensitive to perturbations in $\eta_{D}$, and the standard error
computed by treating ${\widehat{\eta}}_{D}$ as known is incorrect.

An important practical implication is that estimators based on Neyman
orthogonal scores yield inference about $\theta_{0}$ that does not
depend on the detailed statistical properties of the nuisance estimator
$\widehat{\eta}$. This is useful even in classical low-dimensional
settings, where it avoids cumbersome variance adjustments to account for
first-step estimation. It becomes particularly important when modern
flexible methods—including machine learning—are used for nuisance
estimation, as only coarse convergence rates are currently available for
many promising ML methods. By alleviating the first-order impact of
nuisance estimation, Neyman orthogonality makes it possible to combine
such methods with standard asymptotic approximations.

## From OLS to Machine Learning

In our linear example, both nuisance parameters—
$\eta_{Y} = \text{argmin}_{\eta}\, E\left\lbrack (Y - X\prime\eta)^{2} \right\rbrack$
and
$\eta_{D} = \text{argmin}_{\eta}\, E\left\lbrack (D - X\prime\eta)^{2} \right\rbrack$—are
estimated by OLS. Partialling out both $Y$ and $D$ yields the same
estimate and standard error as full OLS of $Y$ on $(D,X)$.

The partially linear regression model generalizes this. Instead of
assuming linearity, we allow $X$ to enter flexibly:

$$Y_{i} = \theta_{0}D_{i} + g_{0}\left( X_{i} \right) + \varepsilon_{i}$$

where $g_{0}( \cdot ) = E\left\lbrack Y - \theta_{0}D|X \right\rbrack$
is an unknown, potentially complex, function of $X$. The nuisance
parameters become $\ell_{0}(X) = E\left\lbrack Y|X \right\rbrack$ and
$r_{0}(X) = E\left\lbrack D|X \right\rbrack$—conditional expectation
functions that can be estimated by machine learning.

The Neyman orthogonal score for the partially linear model is

$$m_{PLM}(W;\,\theta,\eta) = \left\lbrack \left( Y - \ell(X) \right) - \theta\left( D - r(X) \right) \right\rbrack\left( D - r(X) \right),$$

which is the flexible analog of Procedure B’s partialling-out score. As
in our linear example, a naive score that adjusts only $D$ (analogous to
Procedure A) is not Neyman orthogonal and leads to invalid inference.

DML combines the Neyman orthogonal score with *cross-fitting*—a form of
sample splitting that addresses a second source of bias (overfitting
bias) arising when the nuisance estimator $\widehat{\eta}$ and the
observations used in the moment condition are statistically dependent.
Together, these two ingredients form the core of the DML framework.

The `ddml` package implements this approach via `ddml_plm`:

``` r
# Example: DML with gradient boosting (not run)
fit_dml <- ddml_plm(Y, D, X,
                    learners = list(what = mdl_xgboost),
                    sample_folds = 5)
summary(fit_dml)
```

See
[`?ddml_plm`](https://www.thomaswiemann.com/ddml/reference/ddml_plm.md)
for details and
[`vignette("ddml")`](https://www.thomaswiemann.com/ddml/articles/ddml.md)
for a full introduction.

## References

Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
Robins J (2018). “Double/debiased machine learning for treatment and
structural parameters.” The Econometrics Journal, 21(1), C1-C68.

Ahrens A, Chernozhukov V, Hansen C B, Kozbur D, Schaffer M E, Wiemann T
(2025). “An Introduction to Double/Debiased Machine Learning.”
Forthcoming.
