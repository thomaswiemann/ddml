test_that("crossval_compute returns residuals", {
  # Simulate small dataset
  X <- matrix(rnorm(100*100), 100, 100) # Simulate features
  y <- 1 + X %*% (10*runif(100) * (runif(100) < 0.05)) + rnorm(100)
  # Define arguments
  test_sample <- sample(1:length(y), 33)
  learner <- list(what = ols)
  # Compute cross-validation instance
  cv_resid <- crossval_compute(test_sample, learner,
                                y, X)
  # Check output with expectations
  expect_equal(length(cv_resid), 33)
})#TEST_THAT

test_that("crossval returns residuals by learner", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*99), 100, 99)) # Simulate features
  nonzero_X <- (runif(100) < 0.05)
  y <- X %*% (10*runif(100) * nonzero_X) + rnorm(100)
  # Define arguments
  learners <- list(list(what = ols),
                 list(what = ols),
                 list(what = ols,
                      assign_X = which(nonzero_X)))
  # Compute cross-validation instance
  cv_res <- crossval(y, X,
                     learners,
                     cv_folds = 3,
                     silent = TRUE)
  # Check output with expectations
  expect_equal(dim(cv_res$cv_resid), c(length(y), length(learners)))
})#TEST_THAT

test_that("crossval returns residuals by learner in correct order", {
  # Simulate small dataset
  n <- 147
  X <- matrix(rnorm(n * 10), n, 10)
  y <- rowSums(X[, 1:3]) + rnorm(n)
  # split data to two folds and compute residuals manually
  subsample_list <- generate_subsamples(n, 2)
  cv_resid_manual <- matrix(0, n, 2)
  for(i in seq_along(subsample_list)) {
    idx_i <- subsample_list[[i]]
    # ols 1
    ols_fit <- ols(y[-idx_i], X[-idx_i, 1:5])
    cv_resid_manual[idx_i, 1] <- y[idx_i] -
      ddml:::predict.ols(ols_fit, X[idx_i, 1:5])
    # ols 1
    ols_fit <- ols(y[-idx_i], X[-idx_i, 1:10])
    cv_resid_manual[idx_i, 2] <- y[idx_i] -
      ddml:::predict.ols(ols_fit, X[idx_i, 1:10])
  }#FOR
  # Compute cross-validation with crossval using the same subsamples
  cv_res <- crossval(y, X,
                     learners = list(list(what = ols,
                                          assign_X = 1:5),
                                     list(what = ols,
                                          assign_X = 1:10)),
                     cv_subsamples = subsample_list,
                     silent = TRUE)

  # Check output with expectations
  expect_equal(round(cv_res$cv_resid[, 1], 3), round(cv_resid_manual[, 1], 3))
  expect_equal(round(cv_res$cv_resid[, 2], 3), round(cv_resid_manual[, 2], 3))
})#TEST_THAT

test_that("crossval returns identical results with parallel", {
  skip_on_cran()
  skip_if_not_installed("parallel")
  set.seed(42)
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs * 39), nobs, 39))
  y <- X %*% runif(40) + rnorm(nobs)
  learners <- list(list(what = ols),
                   list(what = ols))
  cv_subs <- generate_subsamples(nobs, 3)
  # Sequential
  res_seq <- crossval(y, X, learners = learners,
                      cv_subsamples = cv_subs, silent = TRUE)
  # Parallel
  res_par <- crossval(y, X, learners = learners,
                      cv_subsamples = cv_subs, silent = TRUE,
                      parallel = list(cores = 2))
  expect_equal(res_par$mspe, res_seq$mspe)
  expect_equal(res_par$cv_resid, res_seq$cv_resid)
})#TEST_THAT
