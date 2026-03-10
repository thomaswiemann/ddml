test_that("get_CEF dispatches to crosspred for standard stacking", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  y <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  learners <- list(list(what = ols), list(what = ols))
  splits <- get_sample_splits(seq_len(nobs),
                               sample_folds = 2,
                               cv_folds = 2)
  res <- get_CEF(y, X,
                 learners = learners,
                 ensemble_type = "ols",
                 shortstack = FALSE,
                 subsamples = splits$subsamples,
                 cv_subsamples = splits$cv_subsamples,
                 silent = TRUE)
  # Check standard output components
  expect_true(is.matrix(res$cf_fitted))
  expect_equal(nrow(res$cf_fitted), nobs)
  expect_true(!is.null(res$weights))
  expect_true(!is.null(res$mspe))
})#TEST_THAT

test_that("get_CEF dispatches to shortstacking", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  y <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  learners <- list(list(what = ols), list(what = ols))
  splits <- get_sample_splits(seq_len(nobs), sample_folds = 2)
  res <- get_CEF(y, X,
                 learners = learners,
                 ensemble_type = "ols",
                 shortstack = TRUE,
                 subsamples = splits$subsamples,
                 cv_subsamples = splits$cv_subsamples,
                 silent = TRUE)
  expect_true(is.matrix(res$cf_fitted))
  expect_equal(nrow(res$cf_fitted), nobs)
  expect_true(!is.null(res$cf_resid_bylearner))
})#TEST_THAT

test_that("get_CEF handles constant y gracefully", {
  nobs <- 100
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  y <- rep(5, nobs)
  learners <- list(list(what = ols), list(what = ols))
  splits <- get_sample_splits(seq_len(nobs),
                               sample_folds = 2,
                               cv_folds = 2)
  res <- get_CEF(y, X,
                 learners = learners,
                 ensemble_type = "ols",
                 shortstack = FALSE,
                 subsamples = splits$subsamples,
                 cv_subsamples = splits$cv_subsamples,
                 silent = TRUE)
  # All predictions should be the constant value
  expect_true(all(res$cf_fitted == 5))
  expect_true(all(res$cf_resid_bylearner == 0))
})#TEST_THAT

test_that("get_CEF pass-through with pre-computed fitted works", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  y <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  learners <- list(list(what = ols), list(what = ols))
  splits <- get_sample_splits(seq_len(nobs),
                               sample_folds = 2,
                               cv_folds = 2)
  # First call: compute from scratch
  res1 <- get_CEF(y, X,
                  learners = learners,
                  ensemble_type = "ols",
                  shortstack = FALSE,
                  subsamples = splits$subsamples,
                  cv_subsamples = splits$cv_subsamples,
                  silent = TRUE)
  # Pass-through: Rule 1 path with pre-ensembled cf_fitted
  fitted_preens <- list(cf_fitted = res1$cf_fitted)
  res2 <- get_CEF(y, X,
                  learners = learners,
                  ensemble_type = "ols",
                  shortstack = FALSE,
                  subsamples = splits$subsamples,
                  cv_subsamples = splits$cv_subsamples,
                  silent = TRUE,
                  fitted = fitted_preens)
  expect_equal(res2$cf_fitted, res1$cf_fitted)
  expect_null(res2$weights)
})#TEST_THAT

test_that("get_CEF pass-through with pre-ensembled fitted uses Rule 1", {
  set.seed(42)
  nobs <- 100
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  y <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  learners <- list(what = ols)
  splits <- get_sample_splits(seq_len(nobs),
                               sample_folds = 2,
                               cv_folds = 2)
  # Single-learner call: no cf_fitted_bylearner
  res <- get_CEF(y, X,
                 learners = learners,
                 ensemble_type = "ols",
                 shortstack = FALSE,
                 subsamples = splits$subsamples,
                 cv_subsamples = splits$cv_subsamples,
                 silent = TRUE)
  # Pass-through with pre-ensembled cf_fitted (Rule 1)
  fitted_preensembled <- list(cf_fitted = res$cf_fitted)
  res2 <- get_CEF(y, X,
                  learners = learners,
                  ensemble_type = "ols",
                  shortstack = FALSE,
                  subsamples = splits$subsamples,
                  cv_subsamples = splits$cv_subsamples,
                  silent = TRUE,
                  fitted = fitted_preensembled)
  expect_equal(res2$cf_fitted, res$cf_fitted)
  expect_null(res2$weights)
})#TEST_THAT

test_that("get_CEF Rule 2 with cv_resid_byfold reproduces stacking exactly", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  y <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  learners <- list(list(what = ols), list(what = ols))
  splits <- get_sample_splits(seq_len(nobs),
                               sample_folds = 2,
                               cv_folds = 2)

  for (ens in c("nnls", "nnls1")) {
    res_fresh <- get_CEF(y, X,
                         learners = learners,
                         ensemble_type = ens,
                         shortstack = FALSE,
                         subsamples = splits$subsamples,
                         cv_subsamples = splits$cv_subsamples,
                         silent = TRUE)
    expect_true(!is.null(res_fresh$cv_resid_byfold),
                info = paste(ens, "cv_resid_byfold stored"))

    fitted_r2 <- list(
      cf_fitted_bylearner = res_fresh$cf_fitted_bylearner,
      cv_resid_byfold = res_fresh$cv_resid_byfold)
    res_pt <- get_CEF(y, X,
                      learners = learners,
                      ensemble_type = ens,
                      shortstack = FALSE,
                      subsamples = splits$subsamples,
                      cv_subsamples = splits$cv_subsamples,
                      silent = TRUE,
                      fitted = fitted_r2)
    expect_equal(res_pt$cf_fitted, res_fresh$cf_fitted,
                 tolerance = 1e-10,
                 info = paste(ens, "exact reproduction"))
    expect_true(!is.null(res_pt$weights),
                info = paste(ens, "weights returned"))
  }
})#TEST_THAT

test_that("get_CEF Rule 2 without cv_resid_byfold uses approximate path", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  y <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  learners <- list(list(what = ols), list(what = ols))
  splits <- get_sample_splits(seq_len(nobs),
                               sample_folds = 2,
                               cv_folds = 2)

  # average: exact even without cv_resid_byfold
  res_avg <- get_CEF(y, X,
                     learners = learners,
                     ensemble_type = "average",
                     shortstack = FALSE,
                     subsamples = splits$subsamples,
                     cv_subsamples = splits$cv_subsamples,
                     silent = TRUE)
  fitted_no_cv <- list(
    cf_fitted_bylearner = res_avg$cf_fitted_bylearner)
  res_avg_pt <- get_CEF(y, X,
                        learners = learners,
                        ensemble_type = "average",
                        shortstack = FALSE,
                        subsamples = splits$subsamples,
                        cv_subsamples = splits$cv_subsamples,
                        silent = TRUE,
                        fitted = fitted_no_cv)
  expect_equal(res_avg_pt$cf_fitted, res_avg$cf_fitted,
               tolerance = 1e-10)
  # Branch B produces 2D weights (not per-fold 3D)
  expect_equal(length(dim(res_avg_pt$weights)), 2)

  # nnls1: valid output but approximate (no inner-CV residuals)
  res_nnls <- get_CEF(y, X,
                      learners = learners,
                      ensemble_type = "nnls1",
                      shortstack = FALSE,
                      subsamples = splits$subsamples,
                      cv_subsamples = splits$cv_subsamples,
                      silent = TRUE)
  fitted_nnls_no_cv <- list(
    cf_fitted_bylearner = res_nnls$cf_fitted_bylearner)
  res_nnls_pt <- get_CEF(y, X,
                         learners = learners,
                         ensemble_type = "nnls1",
                         shortstack = FALSE,
                         subsamples = splits$subsamples,
                         cv_subsamples = splits$cv_subsamples,
                         silent = TRUE,
                         fitted = fitted_nnls_no_cv)
  expect_equal(nrow(res_nnls_pt$cf_fitted), nobs)
  expect_true(is.numeric(res_nnls_pt$cf_fitted))
  expect_equal(length(dim(res_nnls_pt$weights)), 2)
})#TEST_THAT

test_that("get_CEF Rule 2 reproduces shortstacking exactly", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  y <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  learners <- list(list(what = ols), list(what = ols))
  splits <- get_sample_splits(seq_len(nobs), sample_folds = 2)

  for (ens in c("average", "nnls1")) {
    res_fresh <- get_CEF(y, X,
                         learners = learners,
                         ensemble_type = ens,
                         shortstack = TRUE,
                         subsamples = splits$subsamples,
                         cv_subsamples = NULL,
                         silent = TRUE)
    expect_null(res_fresh$cv_resid_byfold,
                info = paste(ens, "no cv_resid_byfold"))

    fitted_ss <- list(
      cf_fitted_bylearner = res_fresh$cf_fitted_bylearner)
    res_pt <- get_CEF(y, X,
                      learners = learners,
                      ensemble_type = ens,
                      shortstack = TRUE,
                      subsamples = splits$subsamples,
                      cv_subsamples = NULL,
                      silent = TRUE,
                      fitted = fitted_ss)
    expect_equal(res_pt$cf_fitted, res_fresh$cf_fitted,
                 tolerance = 1e-10,
                 info = paste(ens, "exact reproduction"))
  }
})#TEST_THAT
