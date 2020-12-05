library(ddml)
context("Testing wrapper functions for imported ML procedures.")

sim_data <- function() {
  X <- matrix(rnorm(100*50), 100, 50) # Simulate features
  y <- 1 + X %*% (10*runif(50) * (runif(50) < 0.1)) + rnorm(100)
  return(list(y = y, X = X))
}#SIM_DATA

test_that("mdl_glmnet is working", {
  # Simulate small dataset and fit the model
  dat <- sim_data()
  mdl_fit <- mdl_glmnet(dat$y, dat$X)
  # Check methods predict() and any_iv()
  fitted <- predict(mdl_fit, newdata = dat$X)
  iv_selected <- any_iv(mdl_fit, index_iv = c(11:20))
  # Check output with expectations
  expect_is(mdl_fit, "mdl_glmnet")
  expect_equal(length(fitted), 100)
  expect_is(iv_selected, "logical")
})#TEST_THAT

test_that("mdl_xgboost is working", {
  # Simulate small dataset and fit the model
  dat <- sim_data()
  mdl_fit <- mdl_xgboost(dat$y, dat$X)
  # Check methods predict() and any_iv()
  fitted <- predict(mdl_fit, newdata = dat$X)
  iv_selected <- any_iv(mdl_fit, index_iv = c(11:20),
                        names_iv = colnames(dat$X))
  # Check output with expectations
  expect_is(mdl_fit, "mdl_xgboost")
  expect_equal(length(fitted), 100)
  expect_is(iv_selected, "logical")
})#TEST_THAT

test_that("mdl_randomForest is working", {
  # Simulate small dataset and fit the model
  dat <- sim_data()
  mdl_fit <- mdl_randomForest(dat$y, dat$X)
  # Check methods predict() and any_iv()
  fitted <- predict(mdl_fit, newdata = dat$X)
  iv_selected <- any_iv(mdl_fit, index_iv = c(11:20))
  # Check output with expectations
  expect_is(mdl_fit, "mdl_randomForest")
  expect_equal(length(fitted), 100)
  expect_is(iv_selected, "logical")
})#TEST_THAT

test_that("mdl_grf is working", {
  # Simulate small dataset and fit the model
  dat <- sim_data()
  mdl_fit <- mdl_grf(dat$y, dat$X)
    # Check methods predict() and any_iv()
  fitted <- predict(mdl_fit, newdata = dat$X)
  iv_selected <- any_iv(mdl_fit, index_iv = c(11:20))
  # Check output with expectations
  expect_is(mdl_fit, "mdl_grf")
  expect_equal(length(fitted), 100)
  expect_is(iv_selected, "logical")
})#TEST_THAT
