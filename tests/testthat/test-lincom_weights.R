# Tests for lincom_weights_did ---------------------------------------------

# Shared DGP
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

test_that("did_weights dynamic structure", {
  fit <- make_attgt_fit()
  w <- lincom_weights_did(fit, type = "dynamic")

  C <- nrow(fit$coefficients)
  n <- fit$nobs

  # R is C x q
  expect_equal(nrow(w$R), C)
  expect_true(ncol(w$R) > 0)

  # inf_func_R is n x C x q
  q <- ncol(w$R)
  expect_equal(dim(w$inf_func_R), c(n, C, q))

  # Labels exist and match R columns
  expect_equal(length(w$labels), ncol(w$R))
  expect_equal(colnames(w$R), w$labels)

  # Columns of R sum to 1 (or 0 for excluded)
  col_sums <- colSums(w$R)
  for (s in col_sums) {
    expect_true(abs(s - 1) < 1e-10 || abs(s) < 1e-10)
  }

  # Labels start with "e="
  expect_true(all(grepl("^e=", w$labels)))
})#TEST_THAT

test_that("did_weights group structure", {
  fit <- make_attgt_fit()
  w <- lincom_weights_did(fit, type = "group")

  # Labels start with "g="
  expect_true(all(grepl("^g=", w$labels)))

  # Columns of R sum to 1
  col_sums <- colSums(w$R)
  expect_true(all(abs(col_sums - 1) < 1e-10 |
    abs(col_sums) < 1e-10))
})#TEST_THAT

test_that("did_weights calendar structure", {
  fit <- make_attgt_fit()
  w <- lincom_weights_did(fit, type = "calendar")

  # Labels start with "t="
  expect_true(all(grepl("^t=", w$labels)))
})#TEST_THAT

test_that("did_weights simple structure", {
  fit <- make_attgt_fit()
  w <- lincom_weights_did(fit, type = "simple")

  # Single column
  expect_equal(ncol(w$R), 1)
  expect_equal(w$labels, "ATT")

  # Weights sum to 1
  expect_equal(sum(w$R), 1, tolerance = 1e-10)
})#TEST_THAT

test_that("did_weights min_e/max_e filter", {
  fit <- make_attgt_fit()
  w_all <- lincom_weights_did(fit, type = "dynamic")
  w_sub <- lincom_weights_did(fit, type = "dynamic",
                               min_e = 0, max_e = 1)

  # Filtered has fewer columns
  expect_true(ncol(w_sub$R) <= ncol(w_all$R))

  # All e= labels within range
  evals <- as.numeric(sub("^e=", "", w_sub$labels))
  expect_true(all(evals >= 0))
  expect_true(all(evals <= 1))
})#TEST_THAT

test_that("did_weights rejects non-attgt fit", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- 1 * (rnorm(nobs) > 0)
  y <- D + rnorm(nobs)

  fit <- ddml_ate(y, D, X,
                   learners = list(what = ols),
                   sample_folds = 2, silent = TRUE)

  expect_error(lincom_weights_did(fit),
               "ddml_attgt")
})#TEST_THAT

test_that("did_weights dynamic point estimates match manual", {
  fit <- make_attgt_fit()
  w <- lincom_weights_did(fit, type = "dynamic")
  lc <- lincom(fit, R = w$R,
                     inf_func_R = w$inf_func_R,
                     labels = w$labels)

  # Manual: R'theta
  manual <- as.numeric(crossprod(w$R, coef(fit)))
  expect_equal(as.numeric(coef(lc)), manual,
               tolerance = 1e-10)

  # delta-method was applied
  expect_false(lc$fixed_R)
})#TEST_THAT

test_that("delta method matters for SEs", {
  fit <- make_attgt_fit()
  w <- lincom_weights_did(fit, type = "dynamic")

  # With delta-method correction
  lc_dm <- lincom(fit, R = w$R,
                   inf_func_R = w$inf_func_R,
                   labels = w$labels)

  # Without delta-method (fixed R)
  lc_fixed <- lincom(fit, R = w$R,
                      labels = w$labels)

  # SEs should differ (delta-method adds weight IF)
  se_dm <- sqrt(diag(vcov(lc_dm)))
  se_fixed <- sqrt(diag(vcov(lc_fixed)))
  expect_false(all(abs(se_dm - se_fixed) < 1e-10))
})#TEST_THAT

test_that("simple aggregation equals group-share weighted average", {
  fit <- make_attgt_fit()

  # Native
  w <- lincom_weights_did(fit, type = "simple")
  lc <- lincom(fit, R = w$R,
                     inf_func_R = w$inf_func_R,
                     labels = w$labels)

  # Manual: group-share weighted average of GT-ATTs
  ci <- fit$cell_info
  post <- which(ci$time >= ci$group)
  pg <- ci$n_treated[post] / sum(ci$n_treated[post])
  theta <- coef(fit)[post]
  manual_att <- sum(pg * theta)

  expect_equal(as.numeric(coef(lc)), manual_att, tolerance = 1e-10)
})#TEST_THAT

test_that("dinf_dR has correct dimensions", {
  fit <- make_attgt_fit()
  w <- lincom_weights_did(fit, type = "dynamic")

  n <- fit$nobs
  q <- ncol(w$R)
  expect_equal(dim(w$dinf_dR), c(n, q, q))

  # Constant across observations
  for (i in 2:min(5, n)) {
    expect_equal(w$dinf_dR[i, , ], w$dinf_dR[1, , ],
                 tolerance = 1e-12)
  }
})#TEST_THAT

test_that("dinf_dR is zero for group aggregation", {
  fit <- make_attgt_fit()
  w <- lincom_weights_did(fit, type = "group")

  expect_equal(max(abs(w$dinf_dR)), 0, tolerance = 1e-12)
})#TEST_THAT

test_that("dinf_dR matches manual V'V computation", {
  fit <- make_attgt_fit()
  w <- lincom_weights_did(fit, type = "dynamic")

  ci <- fit$cell_info
  theta <- coef(fit)
  gamma <- as.numeric(crossprod(w$R, theta))
  q <- ncol(w$R)

  glist <- sort(unique(ci$group))
  nG <- length(glist)
  pg <- sapply(glist, function(g) mean(fit$G == g))

  et <- ci$time - ci$group
  all_et <- sort(unique(et))
  V_manual <- matrix(0, nG, q)
  for (k in seq_along(all_et)) {
    e <- all_et[k]
    keepers <- which(et == e)
    S_k <- sum(pg[match(ci$group[keepers], glist)])
    for (j in seq_len(nG)) {
      g <- glist[j]
      cells <- keepers[ci$group[keepers] == g]
      if (length(cells) > 0) {
        V_manual[j, k] <- sum(theta[cells] - gamma[k]) / S_k
      }
    }
  }
  VtV_manual <- crossprod(V_manual)

  expect_equal(w$dinf_dR[1, , ], -VtV_manual, tolerance = 1e-10)
})#TEST_THAT

test_that("HC3 with dinf_dR differs from HC3 without", {
  fit <- make_attgt_fit()
  w <- lincom_weights_did(fit, type = "dynamic")

  lc_with <- lincom(fit, R = w$R,
                    inf_func_R = w$inf_func_R,
                    dinf_dR = w$dinf_dR,
                    labels = w$labels)
  lc_without <- lincom(fit, R = w$R,
                       inf_func_R = w$inf_func_R,
                       labels = w$labels)

  se_with <- sqrt(diag(vcov(lc_with, type = "HC3")))
  se_without <- sqrt(diag(vcov(lc_without, type = "HC3")))

  # Should differ for dynamic (V'V != 0)
  expect_false(all(abs(se_with - se_without) < 1e-10))

  # With dinf_dR should give larger hat values (more leverage)
  h_with <- hatvalues(lc_with)
  h_without <- hatvalues(lc_without)
  expect_true(all(h_with >= h_without - 1e-10))
})#TEST_THAT
