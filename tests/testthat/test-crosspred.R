test_that("crosspred computes with ensemble procedures & custom weights", {
  # generate test data
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  y <- X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols),
                   list(what = ols))
  # Define custom weights
  custom_ensemble_weights <- diag(1, length(learners))
  colnames(custom_ensemble_weights) <- c("mdl_ols1", "mdl_ols2", "mdl_ols3")
  # Compute cross-sample predictions
  crosspred_res <- crosspred(y, X,
                             learners,
                             ensemble_type = c("average", "ols",
                                               "nnls1", "nnls",
                                               "singlebest"),
                             cv_folds = 3,
                             sample_folds = 3,
                             custom_ensemble_weights = custom_ensemble_weights,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(dim(crosspred_res$cf_fitted), c(length(y), 8))
})#TEST_THAT

test_that("crosspred computes with ensemble procedures and sparse matrices", {
  # generate test data
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  y <- X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                 list(what = ols))
  # Compute cross-sample predictions
  crosspred_res <- crosspred(y, as(X, "sparseMatrix"),
                             learners,
                             ensemble_type = c("average", "ols",
                                          "nnls1", "nnls",
                                          "singlebest"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(dim(crosspred_res$cf_fitted), c(length(y), 5))
})#TEST_THAT

test_that("crosspred computes auxilliary predictions", {
  # generate test data
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  y <- X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols),
                   list(what = ols))
  # Compute cross-sample and auxilliary predictions
  crosspred_res <- crosspred(y, X,
                             learners = learners,
                             ensemble_type = c("average", "ols",
                                               "nnls1", "nnls",
                                               "singlebest"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE,
                             auxiliary_X = list(X, X, X))
  # Check output with expectations
  expect_equal(dim(crosspred_res$auxiliary_fitted[[1]]), c(length(y), 5))
})#TEST_THAT

test_that("crosspred returns identical results with parallel", {
  skip_on_cran()
  skip_if_not_installed("parallel")
  set.seed(42)
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs * 39), nobs, 39))
  y <- X %*% runif(40) + rnorm(nobs)
  learners <- list(list(what = ols),
                   list(what = ols))
  splits <- get_sample_splits(seq_len(nobs),
                              sample_folds = 3, cv_folds = 3)
  # Sequential
  res_seq <- crosspred(y, X, learners = learners,
                       ensemble_type = "average",
                       sample_folds = 3, cv_folds = 3,
                       subsamples = splits$subsamples,
                       cv_subsamples = splits$cv_subsamples,
                       silent = TRUE)
  # Parallel
  res_par <- crosspred(y, X, learners = learners,
                       ensemble_type = "average",
                       sample_folds = 3, cv_folds = 3,
                       subsamples = splits$subsamples,
                       cv_subsamples = splits$cv_subsamples,
                       silent = TRUE,
                       parallel = list(cores = 2))
  expect_equal(res_par$cf_fitted, res_seq$cf_fitted)
  expect_equal(res_par$weights, res_seq$weights)
})#TEST_THAT
