library(ddml)
context("Testing lsml objects.")

test_that("lsml converges w/ easy DGP", {
  # Simulate small dataset
  X <- matrix(rnorm(100*20), 100, 20) # Simulate features
  y <- X %*% (10*runif(20)) + rnorm(100)
  mdl_fit <- lsml(y, X,
                  model = mdl_xgboost,
                  args = list(colsample_bytree = 1,
                              subsample = 1,
                              nrounds = 2),
                  split_X = replicate(2, c(1:ncol(X)), F),
                  max_iter = 50, eps = 1e-2,
                  silent = T)
  # Check output with expectations
  expect_is(mdl_fit$ls_fit, "ols")
  expect_is(mdl_fit$ml_fit, "mdl_xgboost")
})#TEST_THAT

test_that("predict.lsml returns fitted values", {
  # Simulate small dataset
  X <- matrix(rnorm(100*20), 100, 20) # Simulate features
  y <- X %*% (10*runif(20)) + rnorm(100)
  mdl_fit <- lsml(y, X,
                  model = mdl_xgboost,
                  args = list(colsample_bytree = 1,
                              subsample = 1,
                              nrounds = 2),
                  split_X = replicate(2, c(1:ncol(X)), F),
                  max_iter = 50, eps = 1e-2,
                  silent = T)
  fitted <- predict(mdl_fit)
  # Check output with expectations
  expect_equal(length(fitted), 100)
})#TEST_THAT
