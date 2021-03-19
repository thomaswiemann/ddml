library(ddml)
context("Testing pooled kcmeans objects.")

sim_data <- function(N = 1000, J = 10) {
  X <- sample(c(1:J), N, replace = T)
  y <- 2 * X + rnorm(N)
  return(list(y = y, X = X))
}#SIM_DATA

test_that("pooled kcmeans computes k conditional means estimator", {
  # Simulate small dataset
  dat <- sim_data()
  # Estimate kcmeans
  mdl_pooled_kcmeans <- pooled_kcmeans(dat$y, dat$X, K = 4, ncluster = 10)
  # Check output with expectations
  expect_equal(length(mdl_pooled_kcmeans$cluster_maps), 10)
  expect_equal(dim(mdl_pooled_kcmeans$alpha_mat), c(10, 4))
})#TEST_THAT

test_that("kcmeans predicts out of sample", {
  # Simulate small dataset
  dat <- sim_data()
  # Estimate kcmeans
  mdl_pooled_kcmeans <- pooled_kcmeans(dat$y, dat$X, K = 4, ncluster = 10)
  # Predict out of sample
  fitted <- predict(mdl_pooled_kcmeans, sample(c(1:10), 100, replace = T))
  # Check output with expectations
  expect_is(fitted, "matrix")
  expect_equal(length(fitted), 100)
})#TEST_THAT
