library(ddml)
context("Testing liml objects.")

test_that("liml returns output of correct types", {
  # Simulate small dataset and fit tsls
  X <- cbind(1, matrix(rnorm(100*3), 100, 3))
  Z <- matrix(rnorm(100*2), 100, 2) # overidentified
  UV <- matrix(rnorm(2*100), 100, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(4) + Z %*% c(1, runif(1)) + UV[, 1]
  y <- D + X %*% runif(4) + UV[, 2]
  mdl_fit <- liml(y, D, Z, X)
  # Check output with expectations
  expect_is(mdl_fit, "liml")
  expect_is(mdl_fit$coef, "matrix") # Return of type double
  expect_equal(dim(mdl_fit$coef), c(1+ncol(X), 1)) # Return a vector
})#TEST_THAT

test_that("predict.liml returns fitted values", {
  # Simulate small dataset and fit tsls
  X <- cbind(1, matrix(rnorm(100*3), 100, 3));
  Z <- matrix(rnorm(100*2), 100, 2) # overidentified
  UV <- matrix(rnorm(2*100), 100, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(4) + Z %*% c(1, runif(1)) + UV[, 1]
  y <- D + X %*% runif(4) + UV[, 2]
  mdl_fit <- liml(y, D, Z, X)
  # Check in- and out-of-sample fitted values
  expect_equal(predict(mdl_fit),  cbind(D, X) %*% mdl_fit$coef)
})#TEST_THAT

test_that("summary.liml returns numerical se, t-stats, and p-values", {
  # Simulate small dataset and fit tsls
  nobs <- 10000
  X <- cbind(1, matrix(rnorm(nobs*3), nobs, 3));
  Z <- matrix(rnorm(nobs*2), nobs, 2) # overidentified
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(4) + Z %*% c(1, runif(1)) + UV[, 1]
  y <- D + X %*% runif(4) + UV[, 2]
  mdl_fit <- liml(y, D, Z, X)
  sum_res_const <- summary(mdl_fit, type = "const")
  sum_res_CSE <- summary(mdl_fit, type = "CSE")
  # Check in- and out-of-sample fitted values
  expect_is(sum_res_const$res, "matrix")
  expect_is(sum_res_CSE$res, "matrix")
  expect_equal(dim(sum_res_const$res), c(1+ncol(X), 4))
  expect_equal(dim(sum_res_CSE$res), c(1+ncol(X), 4))
})#TEST_THAT
