library(ddml)
context("Testing supplementary helper functions.")

test_that("csolve returns a (generalized) inverse", {
  # Simulate small matrices, one of them singular
  X <- matrix(rnorm(30*30), 30, 30)
  X_sing <- rbind(cbind(X[1:15, 1:15], X[1:15, 1:15]),
                  cbind(X[1:15, 1:15], X[1:15, 1:15]))
  # Check combuted inverses
  expect_equal(csolve(X) %*% X, diag(1, 30))
  expect_equal((csolve(X_sing) %*% X_sing) %*% X_sing, X_sing)
})#TEST_THAT

test_that("ccov returns a covariance matrix", {
  # Simulate small matrices, one of them singular
  X <- matrix(rnorm(6*100), 100, 6)
  # Check combuted inverses
  expect_equal(ccov(X), cov(X))
})#TEST_THAT

test_that("ccor returns a correlation matrix", {
  # Simulate small matrices, one of them singular
  X <- matrix(rnorm(6*100), 100, 6)
  # Check combuted inverses
  expect_equal(ccor(X), cor(X))
  expect_equal(ccor(X, X[, 1:2]), cor(X, X[, 1:2]))
})#TEST_THAT
