library(ddml)
context("Testing parallel computation of crossval.")

# Skip tests if parallel processing not possible
check_cores <- function() {
  cores <- tryCatch(expr = {parallel::detectCores(logical = FALSE)},
                    error = function(e) {0})
  if (cores < 2)
    skip("Multiple CPU cores not available")
}

test_that("crossval works w/ parallelization w/ static job scheduling", {
  check_cores()
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*99), 100, 99)) # Simulate features
  nonzero_X <- (runif(100) < 0.05)
  y <- X %*% (10*runif(100) * nonzero_X) + rnorm(100)
  # Define arguments
  models <- list(list(fun = rlasso,
                      args = list(include = NULL,
                                  iter_resid = 1, d = 5),
                      assign_X = c(1:ncol(X))), # rlasso
                 list(fun = ols,
                      args = list(),
                      assign_X = c(1:ncol(X))), # ols w/ all features
                 list(fun = ols,
                      args = list(),
                      assign_X = which(nonzero_X))) # ols w/ important features
  # Compute cross-validation instance
  cv_res <- crossval(y, X, Z = NULL,
                     models,
                     cv_folds = 3,
                     setup_parallel = list(type = "static",
                                           cores = 2),
                     silent = T)
  # Check output with expectations
  expect_equal(dim(cv_res$oos_resid), c(length(y), length(models)))
})#TEST_THAT

test_that("crossval works w/ parallelization w/ dynamic job scheduling", {
  check_cores()
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*99), 100, 99)) # Simulate features
  nonzero_X <- (runif(100) < 0.05)
  y <- X %*% (10*runif(100) * nonzero_X) + rnorm(100)
  # Define arguments
  models <- list(list(fun = rlasso,
                      args = list(include = NULL,
                                  iter_resid = 1, d = 5),
                      assign_X = c(1:ncol(X))), # rlasso
                 list(fun = ols,
                      args = list(),
                      assign_X = c(1:ncol(X))), # ols w/ all features
                 list(fun = ols,
                      args = list(),
                      assign_X = which(nonzero_X))) # ols w/ important features
  # Compute cross-validation instance
  cv_res <- crossval(y, X, Z = NULL,
                     models,
                     cv_folds = 3,
                     setup_parallel = list(type = "dynamic",
                                           cores = 2),
                     silent = T)
  # Check output with expectations
  expect_equal(dim(cv_res$oos_resid), c(length(y), length(models)))
})#TEST_THAT
