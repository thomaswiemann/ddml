library(ddml)
context("Testing crossval.")

test_that("crossval_compute returns residuals (w/o instruments)", {
  # Simulate small dataset
  X <- matrix(rnorm(100*100), 100, 100) # Simulate features
  y <- 1 + X %*% (10*runif(100) * (runif(100) < 0.05)) + rnorm(100)
  # Define arguments
  test_sample <- sample(c(1:length(y)), 33)
  model <- list(fun = rlasso,
                args = list(include = c(1:10),
                            iter_resid = 1, d = 5))
  # Compute cross-validation instance
  cv_res <- crossval_compute(test_sample, model,
                             y, X, Z = NULL)
  # Check output with expectations
  expect_equal(length(cv_res$oos_resid), 33)
  expect_equal(cv_res$cv_Z, TRUE)
})#TEST_THAT

test_that("crossval returns residuals by model (w/o instruments)", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*99), 100, 99)) # Simulate features
  nonzero_X <- (runif(100) < 0.05)
  y <- X %*% (10*runif(100) * nonzero_X) + rnorm(100)
  # Define arguments
  models <- list(list(fun = rlasso,
                      args = list(iter_resid = 1, d = 5)), # rlasso
                 list(fun = ols), # ols w/ all features
                 list(fun = ols,
                      assign_X = which(nonzero_X))) # ols w/ important features
  # Compute cross-validation instance
  cv_res <- crossval(y, X, Z = NULL,
                     models,
                     cv_folds = 3,
                     silent = T)
  # Check output with expectations
  expect_equal(dim(cv_res$oos_resid), c(length(y), length(models)))
})#TEST_THAT

test_that("crossval returns residuals by model (w/ instruments)", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  Z <- matrix(rnorm(100*10), 100, 10)
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(100)
  # Define arguments
  models <- list(list(fun = rlasso,
                      args = list(include = c(35:45),
                                  iter_resid = 1, d = 5)), # rlasso
                 list(fun = ols)) # ols w/ important features
  # Compute cross-validation instance
  cv_res <- crossval(D, X, Z,
                     models,
                     cv_folds = 3,
                     silent = T)
  # Check output with expectations
  expect_equal(length(cv_res), 2)
})#TEST_THAT
