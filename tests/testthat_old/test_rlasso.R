library(ddml)
context("Testing rlasso objects.")

test_that("rlasso_penalty computes penalty terms", {
  # Simulate small dataset
  X <- matrix(rnorm(100*100), 100, 100) # Simulate features
  y <- 1 + X %*% (10*runif(100) * (runif(100) < 0.05)) + rnorm(100)
  penalty_const <- rlasso_penalty(y, X,
                 include = c(1:10),
                 iter_resid = 5, d = 5,
                 HC_robust = FALSE)
  penalty_HC_robust <- rlasso_penalty(y, X,
                              include = c(1:10),
                              iter_resid = 5, d = 5,
                              HC_robust = T)
  # Check output with expectations
  expect_is(penalty_const$lambda, "numeric")
  expect_is(penalty_const$psi, "matrix")
  expect_is(penalty_const$sigma, "numeric")
  expect_equal(penalty_const$psi[1:10], rep(0, 10))
  expect_is(penalty_HC_robust$lambda, "numeric")
  expect_is(penalty_HC_robust$psi, "matrix")
  expect_is(penalty_HC_robust$sigma, "numeric")
  expect_equal(penalty_HC_robust$psi[1:10], rep(0, 10))
})#TEST_THAT

test_that("rlasso computes lasso and post lasso coefficients", {
  # Simulate small dataset
  X <- matrix(rnorm(100*100), 100, 100) # Simulate features
  y <- 1 + X %*% (10*runif(100) * (runif(100) < 0.05)) + rnorm(100)
  mdl_rlasso <- rlasso(y, X,
                    include = c(1:10),
                    post = TRUE,
                    partial = FALSE,
                    iter_resid = 1, d = 5) # Compute rlasso
  # Check output with expectations
  expect_equal(all(mdl_rlasso$coef_lasso[1:10] != 0), T)
  expect_equal(all(mdl_rlasso$coef_postlasso[1:10] != 0), T)
  expect_is(mdl_rlasso$coef_lasso, c("dgCMatrix", "matrix"))
  expect_is(mdl_rlasso$coef_postlasso, c("dgCMatrix", "matrix"))
})#TEST_THAT
