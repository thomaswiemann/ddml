test_that("shortstacking computes with ensemble procedures & custom weights", {
  # generate test data
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10) # overidentified
  y <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols),
                   list(fun = ols))
  # Compute cross-sample predictions
  shortstacking_res <- shortstacking(y, X, Z,
                             learners,
                             ensemble_type = c("average", "ols",
                                               "nnls1", "nnls",
                                               "singlebest"),
                             custom_ensemble_weights = diag(1, 3),
                             sample_folds = 3,
                             compute_insample_predictions = T,
                             silent = F)
  # Check output with expectations
  expect_equal(dim(shortstacking_res$oos_fitted), c(length(y), 8))
  expect_equal(length(shortstacking_res$is_fitted), 8)
})#TEST_THAT
