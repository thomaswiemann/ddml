library(ddml)
context("Testing kcmeans objects.")

sim_data <- function(N = 1000, J = 10) {
  X <- sample(c(1:J), N, replace = T)
  y <- 2 * X + rnorm(N)
  return(list(y = y, X = X))
}#SIM_DATA

test_that("init_kcmeans computes initial clustering", {
  # Simulate small dataset
  dat <- sim_data()
  # Calculate initial clustering
  alpha_0 <- init_kcmeans(dat$y, dat$X, K = 4)
  # Check output with expectations
  expect_is(alpha_0, "numeric")
  expect_equal(length(alpha_0), 4)
})#TEST_THAT

test_that("kcmeans computes k conditional means estimator", {
  # Simulate small dataset
  dat <- sim_data()
  # Estimate kcmeans
  mdl_kcmeans <- kcmeans(dat$y, dat$X, K = 4)
  # Check output with expectations
  expect_is(mdl_kcmeans$alpha, "numeric")
  expect_equal(length(mdl_kcmeans$alpha), 4)
})#TEST_THAT

test_that("kcmeans predicts out of sample", {
  # Simulate small dataset
  dat <- sim_data()
  # Estimate kcmeans
  mdl_kcmeans <- kcmeans(dat$y, dat$X, K = 4)
  # Predict out of sample
  fitted <- predict(mdl_kcmeans, sample(c(1:10), 100, replace = T))
  # Check output with expectations
  expect_is(fitted, "matrix")
  expect_equal(length(fitted), 100)
})#TEST_THAT
