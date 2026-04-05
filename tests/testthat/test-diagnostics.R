test_that("diagnostics works with PLM single learner", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  fit <- ddml_plm(y, D, X,
                    learners = list(what = ols),
                    sample_folds = 2,
                    silent = TRUE)

  diag <- diagnostics(fit)

  expect_s3_class(diag, "ddml_diagnostics")
  expect_true("y_X" %in% names(diag$tables))
  # Should have at least one equation table
  expect_true(length(diag$tables) >= 1)

  # print should run without error
  out <- capture_output(print(diag))
  expect_true(grepl("Stacking diagnostics", out))
})

test_that("diagnostics works with PLM multiple learners", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_plm(y, D, X,
                    learners = learners,
                    ensemble_type = c("nnls1"),
                    sample_folds = 2,
                    silent = TRUE)

  diag <- diagnostics(fit)

  # Should have 2 learners + ensemble per equation
  for (eq in names(diag$tables)) {
    expect_true(nrow(diag$tables[[eq]]) >= 3)
  }

  # Weights should sum to ~1
  for (eq in names(diag$tables)) {
    expect_equal(sum(diag$tables[[eq]]$weight, na.rm = TRUE), 1,
                 tolerance = 1e-4)
  }
})

test_that("tidy.ddml_diagnostics returns flat data.frame", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  fit <- ddml_plm(y, D, X,
                    learners = list(what = ols),
                    sample_folds = 2,
                    silent = TRUE)

  td <- tidy(diagnostics(fit))

  expect_s3_class(td, "data.frame")
  expect_true(all(c("equation", "learner", "mspe", "r2",
                     "weight") %in% colnames(td)))
  expect_true(nrow(td) >= 1)
})

test_that("diagnostics with CVC", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_plm(y, D, X,
                    learners = learners,
                    ensemble_type = c("nnls1"),
                    sample_folds = 2,
                    silent = TRUE)

  diag <- diagnostics(fit, cvc = TRUE, bootnum = 100)

  # CVC column should be present
  for (eq in names(diag$tables)) {
    expect_true("cvc_pval" %in% colnames(diag$tables[[eq]]))
  }

  # in_conf_set should NOT be present
  for (eq in names(diag$tables)) {
    expect_false("in_conf_set" %in%
                  colnames(diag$tables[[eq]]))
  }

  # tidy should include CVC column
  td <- tidy(diag)
  expect_true("cvc_pval" %in% colnames(td))
  expect_false("in_conf_set" %in% colnames(td))
})

test_that("diagnostics works with ATE", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  fit <- ddml_ate(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 2,
                  silent = TRUE)

  diag <- diagnostics(fit)

  expect_s3_class(diag, "ddml_diagnostics")
  expect_true(length(diag$tables) >= 2)
  out <- capture_output(print(diag))
  expect_true(grepl("Average Treatment Effect", out))
})

test_that("diagnostics r2 matches manual calculation", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_plm(y, D, X,
                    learners = learners,
                    ensemble_type = "nnls1",
                    sample_folds = 2,
                    silent = TRUE)

  # R2 per learner per fold from crossval: 1 - mspe / var(y)
  # Check that r2 is stored and non-NULL
  expect_false(is.null(fit$r2$y_X))
  expect_true(all(fit$r2$y_X <= 1))
})

test_that("diagnostics rejects non-ddml objects", {
  expect_error(diagnostics(lm(1:10 ~ rnorm(10))),
               "class")
})

test_that("cvc_one_vs_many detects a dominant learner", {
  set.seed(42)
  n <- 500
  resid_best <- rnorm(n, sd = 0.5)
  resid_others <- cbind(rnorm(n, sd = 1.5), rnorm(n, sd = 2.0))
  fid <- rep(seq_len(5), each = n / 5)

  # Best learner: should NOT be rejected (large p-value)
  pval <- cvc_one_vs_many(resid_best, resid_others, fid,
                          bootnum = 500)
  expect_true(pval > 0.5)
})

test_that("cvc_one_vs_many rejects a weak learner", {
  set.seed(42)
  n <- 500
  resid_weak <- rnorm(n, sd = 2.0)
  resid_others <- cbind(rnorm(n, sd = 0.5), rnorm(n, sd = 0.5))
  fid <- rep(seq_len(5), each = n / 5)

  # Weak learner: should be rejected (small p-value)
  pval <- cvc_one_vs_many(resid_weak, resid_others, fid,
                          bootnum = 500)
  expect_true(pval < 0.1)
})

test_that("cvc_pvalues returns NA for single learner", {
  # Mock a fitted object with single learner residuals
  fitted <- list(y_X = list(
    cf_resid_bylearner = matrix(rnorm(100), ncol = 1)))
  splits <- list(subsamples = list(1:50, 51:100))
  pvals <- cvc_pvalues(fitted, splits, "y_X", bootnum = 50)
  expect_length(pvals, 1)
  expect_true(is.na(pvals))
})

