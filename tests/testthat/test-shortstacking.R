test_that("shortstacking computes with ensemble procedures & custom weights", {
  # generate test data
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  y <- X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols),
                   list(what = ols))
  # Compute cross-sample predictions
  shortstacking_res <- shortstacking(y, X,
                             learners,
                             ensemble_type = c("average", "ols",
                                               "nnls1", "nnls",
                                               "singlebest"),
                             custom_ensemble_weights = diag(1, 3),
                             sample_folds = 3,
                             silent = FALSE)
  # Check output with expectations
  expect_equal(dim(shortstacking_res$cf_fitted), c(length(y), 8))
})#TEST_THAT
