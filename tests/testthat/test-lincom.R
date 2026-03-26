# Tests for lincom =============================================================

# Shared DGP: ddml_attgt fit with known structure
make_attgt_fit <- function(seed = 42, n = 800, T_ = 4) {
  set.seed(seed)
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)
  for (i in seq_len(n)) {
    if (is.finite(G[i])) {
      for (tt in seq_len(T_)) {
        if (tt >= G[i]) y[i, tt] <- y[i, tt] + 1
      }
    }
  }
  ddml_attgt(y, X, t = 1:T_, G = G,
            learners = list(what = ols),
            sample_folds = 2, cv_folds = 3,
            silent = TRUE)
}#MAKE_ATTGT_FIT

# lincom ======================================================================

test_that("lincom identity R recovers original", {
  fit <- make_attgt_fit()
  p <- nrow(fit$coefficients)
  R <- diag(p)
  lc <- lincom(fit, R = R)

  # Class: lincom > ral
  expect_s3_class(lc, "lincom")
  expect_true(inherits(lc, "ral"))

  # Coefficients match
  expect_equal(as.numeric(coef(lc)),
               as.numeric(coef(fit)),
               tolerance = 1e-10)

  # Vcov match (dimnames differ, use unname)
  V_lc <- vcov(lc)
  V_fit <- vcov(fit)
  expect_equal(unname(V_lc), unname(V_fit),
               tolerance = 1e-8)

  # nobs
  expect_equal(lc$nobs, fit$nobs)

  # fixed_R = TRUE (no delta-method)
  expect_true(lc$fixed_R)
})#TEST_THAT

test_that("lincom scalar contrast works", {
  fit <- make_attgt_fit()
  p <- nrow(fit$coefficients)

  # Contrast: first cell minus second
  R <- matrix(0, p, 1)
  R[1, 1] <- 1
  R[2, 1] <- -1
  lc <- lincom(fit, R = R, labels = "ATT1-ATT2")

  expect_length(coef(lc), 1)
  expect_equal(names(coef(lc)), "ATT1-ATT2")

  # Point estimate = difference of first two coefs
  expected <- coef(fit)[1] - coef(fit)[2]
  expect_equal(as.numeric(coef(lc)),
               as.numeric(expected),
               tolerance = 1e-10)

  # Vcov is 1x1, positive
  V <- vcov(lc)
  expect_equal(dim(V), c(1, 1))
  expect_true(V[1, 1] > 0)
})#TEST_THAT

test_that("vcov of lincom is PSD with positive diagonal", {
  fit <- make_attgt_fit()
  p <- nrow(fit$coefficients)
  R <- diag(p)
  lc <- lincom(fit, R = R)

  V <- vcov(lc)
  expect_true(all(diag(V) > 0))

  # PSD: eigenvalues >= 0
  evals <- eigen(V, only.values = TRUE)$values
  expect_true(all(evals >= -1e-10))
})#TEST_THAT

test_that("confint lower < upper", {
  fit <- make_attgt_fit()
  p <- nrow(fit$coefficients)
  R <- diag(p)
  lc <- lincom(fit, R = R)

  ci <- confint(lc)
  expect_equal(ncol(ci), 2)
  expect_equal(nrow(ci), p)
  expect_true(all(ci[, 1] < ci[, 2]))

  # Pointwise has Gaussian quantile
  expect_equal(attr(ci, "crit_val"), qnorm(0.975), tolerance = 1e-6)
})#TEST_THAT

test_that("confint uniform wider than pointwise", {
  fit <- make_attgt_fit()
  p <- nrow(fit$coefficients)
  R <- diag(p)
  lc <- lincom(fit, R = R)

  ci_pw <- confint(lc)

  set.seed(1)
  ci_uf <- confint(lc, uniform = TRUE, bootstraps = 499)

  # Uniform bands >= pointwise (p > 1)
  width_pw <- ci_pw[, 2] - ci_pw[, 1]
  width_uf <- ci_uf[, 2] - ci_uf[, 1]
  expect_true(all(width_uf >= width_pw))

  # Critical value attribute
  expect_false(is.null(attr(ci_uf, "crit_val")))
  expect_true(attr(ci_uf, "crit_val") >
    qnorm(0.975))
})#TEST_THAT

test_that("lincom dimension mismatch errors", {
  fit <- make_attgt_fit()
  p <- nrow(fit$coefficients)

  # Wrong number of rows
  R_bad <- matrix(1, p + 1, 1)
  expect_error(lincom(fit, R = R_bad))

  # inf_func_R wrong dimensions
  R <- diag(p)
  inf_func_R_bad <- matrix(0, fit$nobs, p + 1)
  expect_error(lincom(fit, R = R,
    inf_func_R = inf_func_R_bad))
})#TEST_THAT

test_that("summary returns summary.ddml via inheritance", {
  fit <- make_attgt_fit()
  p <- nrow(fit$coefficients)
  R <- diag(p)
  lc <- lincom(fit, R = R)

  s <- summary(lc)
  # Returns summary.ral via inheritance
  expect_true(inherits(s, "summary.ral"))
  expect_true(is.array(s$coefficients))
  expect_equal(s$nobs, fit$nobs)

  # Print does not error
  capture_output({print(s)}, print = FALSE)
})#TEST_THAT

test_that("hatvalues with dinf_dtheta works for lincom", {
  fit <- make_attgt_fit()
  p <- nrow(fit$coefficients)
  R <- diag(p)
  lc <- lincom(fit, R = R)

  # hatvalues should work (dinf_dtheta computed in lincom)
  h <- hatvalues(lc)
  expect_length(h, fit$nobs)
})#TEST_THAT

test_that("tidy works via inheritance", {
  fit <- make_attgt_fit()
  p <- nrow(fit$coefficients)
  R <- diag(p)
  lc <- lincom(fit, R = R)

  td <- tidy(lc)
  expect_true(is.data.frame(td))
  expect_equal(nrow(td), p)
  expect_true("estimate" %in% names(td))
  expect_true("std.error" %in% names(td))
})#TEST_THAT

test_that("print.lincom works", {
  fit <- make_attgt_fit()
  p <- nrow(fit$coefficients)
  R <- diag(p)
  lc <- lincom(fit, R = R)

  out <- capture_output({print(lc)}, print = FALSE)
  expect_true(grepl("Linear Combination", out))
  expect_true(grepl("Obs:", out))
})#TEST_THAT

# lincom_rep ===================================================================

test_that("lincom_rep inherits from ral_rep", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)
  for (i in seq_len(n)) {
    if (is.finite(G[i])) {
      for (tt in seq_len(T_)) {
        if (tt >= G[i]) y[i, tt] <- y[i, tt] + 1
      }
    }
  }

  reps <- ddml_replicate(
    ddml_attgt, y = y, X = X,
    t = 1:T_, G = G,
    learners = list(what = ols),
    sample_folds = 2,
    resamples = 3, silent = TRUE)

  p <- nrow(reps$fits[[1]]$coefficients)
  R <- diag(p)
  lc <- lincom(reps, R = R)

  # Class hierarchy
  expect_s3_class(lc, "lincom_rep")
  expect_true(inherits(lc, "ral_rep"))
  expect_equal(lc$nresamples, 3)

  expect_s3_class(lc$fits[[1]], "lincom")
  expect_true(inherits(lc$fits[[1]], "ral"))

  # coef, vcov, confint via inheritance
  expect_length(coef(lc), p)
  V <- vcov(lc)
  expect_equal(dim(V), c(p, p))
  expect_true(all(diag(V) > 0))

  ci <- confint(lc)
  expect_equal(nrow(ci), p)
  expect_true(all(ci[, 1] < ci[, 2]))
})#TEST_THAT

test_that("rep vcov matches inflate-then-median", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)

  reps <- ddml_replicate(
    ddml_attgt, y = y, X = X,
    t = 1:T_, G = G,
    learners = list(what = ols),
    sample_folds = 2,
    resamples = 3, silent = TRUE)

  p <- nrow(reps$fits[[1]]$coefficients)
  R <- diag(p)
  lc <- lincom(reps, R = R)

  # Manual inflate-then-median
  agg_coef <- coef(lc)
  V_arr <- array(NA, dim = c(p, p, 3))
  for (r in seq_len(3)) {
    Sigma_r <- vcov(lc$fits[[r]])
    bdiff <- coef(lc$fits[[r]]) - agg_coef
    V_arr[, , r] <- Sigma_r + tcrossprod(bdiff)
  }
  V_manual <- apply(V_arr, c(1, 2), median)
  dimnames(V_manual) <- dimnames(vcov(lc))

  expect_equal(vcov(lc), V_manual, tolerance = 1e-10)
})#TEST_THAT

test_that("uniform CI for lincom_rep", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)

  reps <- ddml_replicate(
    ddml_attgt, y = y, X = X,
    t = 1:T_, G = G,
    learners = list(what = ols),
    sample_folds = 2,
    resamples = 3, silent = TRUE)

  p <- nrow(reps$fits[[1]]$coefficients)
  R <- diag(p)
  lc <- lincom(reps, R = R)

  ci_pw <- confint(lc)

  set.seed(1)
  ci_uf <- confint(lc, uniform = TRUE, bootstraps = 499)

  width_pw <- ci_pw[, 2] - ci_pw[, 1]
  width_uf <- ci_uf[, 2] - ci_uf[, 1]
  expect_true(all(width_uf >= width_pw))

  expect_false(is.null(attr(ci_uf, "crit_val")))
})#TEST_THAT

test_that("summary/print for rep via inheritance", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)

  reps <- ddml_replicate(
    ddml_attgt, y = y, X = X,
    t = 1:T_, G = G,
    learners = list(what = ols),
    sample_folds = 2,
    resamples = 3, silent = TRUE)

  p <- nrow(reps$fits[[1]]$coefficients)
  R <- diag(p)
  lc <- lincom(reps, R = R)

  s <- summary(lc)
  expect_true(inherits(s, "summary.ral_rep"))

  out <- capture_output({print(lc)}, print = FALSE)
  expect_true(grepl("Resamples:", out))

  capture_output({print(s)}, print = FALSE)
})#TEST_THAT

test_that("lincom.ddml_rep rejects non-rep fit", {
  fit <- make_attgt_fit()
  R <- diag(nrow(fit$coefficients))
  expect_error(lincom(fit, R = R), NA)  # ddml works
})#TEST_THAT

test_that("lincom_rep HC3 vcov works now", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)

  reps <- ddml_replicate(
    ddml_attgt, y = y, X = X,
    t = 1:T_, G = G,
    learners = list(what = ols),
    sample_folds = 2,
    resamples = 3, silent = TRUE)

  p <- nrow(reps$fits[[1]]$coefficients)
  R <- diag(p)
  lc <- lincom(reps, R = R)

  # HC3 should now work for lincom via ral
  V_hc3 <- vcov(lc, type = "HC3")
  expect_equal(dim(V_hc3), c(p, p))
  expect_true(all(diag(V_hc3) > 0))
})#TEST_THAT

# Multi-ensemble lincom ======================================================

test_that("lincom all ensembles with fit_idx = NULL", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)
  for (i in seq_len(n)) {
    if (is.finite(G[i])) {
      for (tt in seq_len(T_)) {
        if (tt >= G[i]) y[i, tt] <- y[i, tt] + 1
      }
    }
  }
  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                   learners = list(
                     list(what = ols),
                     list(what = ols,
                          args = list(const = FALSE))),
                   ensemble_type = c("nnls", "singlebest"),
                   sample_folds = 2, cv_folds = 3,
                   silent = TRUE)

  p <- nrow(fit$coefficients)
  nensb <- ncol(fit$coefficients)
  R <- diag(p)

  # Default fit_idx = NULL -> all ensembles
  lc <- lincom(fit, R = R)
  expect_equal(ncol(lc$coefficients), nensb)
  expect_equal(dim(lc$inf_func)[3], nensb)

  # Each ensemble matches single-ensemble lincom
  for (j in seq_len(nensb)) {
    lc_j <- lincom(fit, R = R, fit_idx = j)
    expect_equal(lc$coefficients[, j],
                 lc_j$coefficients[, 1],
                 tolerance = 1e-10)
  }

  # summary works with multiple ensembles
  s <- summary(lc)
  expect_equal(dim(s$coefficients)[3], nensb)
})#TEST_THAT

test_that("lincom_weights_did multi-ensemble dinf_dR", {
  fit <- make_attgt_fit()

  # Single-ensemble fit: only 1 ensemble, dinf_dR is 3D
  w <- lincom_weights_did(fit, type = "dynamic")
  expect_equal(length(dim(w$dinf_dR)), 3L)
  q <- ncol(w$R)
  expect_equal(dim(w$dinf_dR), c(fit$nobs, q, q))

  # Explicit fit_idx = 1: same shape
  w1 <- lincom_weights_did(fit, type = "dynamic", fit_idx = 1)
  expect_equal(dim(w1$dinf_dR), c(fit$nobs, q, q))
  expect_equal(w$dinf_dR, w1$dinf_dR, tolerance = 1e-12)
})#TEST_THAT

test_that("as.list on lincom with multi-ensemble", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)
  for (i in seq_len(n)) {
    if (is.finite(G[i])) {
      for (tt in seq_len(T_)) {
        if (tt >= G[i]) y[i, tt] <- y[i, tt] + 1
      }
    }
  }
  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                   learners = list(
                     list(what = ols),
                     list(what = ols,
                          args = list(const = FALSE))),
                   ensemble_type = c("nnls", "singlebest"),
                   sample_folds = 2, cv_folds = 3,
                   silent = TRUE)

  p <- nrow(fit$coefficients)
  lc <- lincom(fit, R = diag(p))

  L <- as.list(lc)
  expect_length(L, 2)
  for (j in 1:2) {
    expect_s3_class(L[[j]], "lincom")
    expect_equal(ncol(L[[j]]$coefficients), 1L)
  }
})#TEST_THAT
