library(ddml)
context("Testing ols objects.")

# ols ==========================================================================
test_that("ols returns output of correct types", {
  # Simulate small dataset and fit ols
  X <- matrix(rnorm(100*3), 100, 3)
  y <- 1 + X %*% c(-1, 1, 0) + rnorm(100)
  mdl_fit <- ols(y, X, const = T)
  # Check output with expectations
  expect_is(mdl_fit, "ols")
  expect_is(mdl_fit$coef, "matrix") # Return of type double
  expect_equal(dim(mdl_fit$coef), c(1+ncol(X), 1)) # Return a vector
  expect_is(mdl_fit$w, "NULL") # Return NULL weights
})#TEST_THAT

test_that("predict.ols returns fitted values", {
  # Simulate small dataset and fit ols
  X1 <- matrix(rnorm(100*3), 100, 3); X2 <- matrix(rnorm(100*3), 100, 3)
  y <- 1 + X1 %*% c(-1, 1, 0) + rnorm(100)
  mdl_fit <- ols(y, X1, const = T, w = runif(100))
  # Check in- and out-of-sample fitted values
  expect_equal(predict(mdl_fit), cbind(1, X1) %*% mdl_fit$coef)
  expect_equal(predict(mdl_fit, newdata = X2), cbind(1, X2) %*% mdl_fit$coef)
})#TEST_THAT

test_that("summary.ols returns numerical se, t-stats, and p-values", {
  # Simulate small dataset and fit ols
  X <- matrix(rnorm(100*3), 100, 3)
  y <- 1 + X %*% c(-1, 1, 0) + rnorm(100)
  mdl_fit <- ols(y, X, const = T)
  sum_res_const <- summary(mdl_fit, type = "const")
  sum_res_HC1 <- summary(mdl_fit, type = "HC1")
  # Check in- and out-of-sample fitted values
  expect_is(sum_res_const$res, "matrix")
  expect_is(sum_res_HC1$res, "matrix")
  expect_equal(dim(sum_res_const$res), c(1+ncol(X), 4))
  expect_equal(dim(sum_res_HC1$res), c(1+ncol(X), 4))
})#TEST_THAT

# wls ==========================================================================
test_that("wls returns output of correct types", {
  # Simulate small dataset and fit ols
  X <- matrix(rnorm(100*3), 100, 3)
  y <- 1 + X %*% c(-1, 1, 0) + rnorm(100)
  mdl_fit <- ols(y, X, const = T, w = runif(100))
  # Check output with expectations
  expect_is(mdl_fit, "ols")
  expect_is(mdl_fit$coef, "matrix") # Return of type double
  expect_equal(dim(mdl_fit$coef), c(1+ncol(X), 1)) # Return a vector
  expect_is(mdl_fit$w, "numeric") # Return passthrough weights
})#TEST_THAT
