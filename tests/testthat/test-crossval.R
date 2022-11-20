test_that("crossval_compute returns residuals (w/o instruments)", {
  # Simulate small dataset
  X <- matrix(rnorm(100*100), 100, 100) # Simulate features
  y <- 1 + X %*% (10*runif(100) * (runif(100) < 0.05)) + rnorm(100)
  # Define arguments
  test_sample <- sample(c(1:length(y)), 33)
  learner <- list(fun = mdl_glmnet)
  # Compute cross-validation instance
  oos_resid <- crossval_compute(test_sample, learner,
                             y, X, Z = NULL)
  # Check output with expectations
  expect_equal(length(oos_resid), 33)
})#TEST_THAT

test_that("crossval returns residuals by learner (w/o instruments)", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*99), 100, 99)) # Simulate features
  nonzero_X <- (runif(100) < 0.05)
  y <- X %*% (10*runif(100) * nonzero_X) + rnorm(100)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet), # lasso
                 list(fun = ols), # ols w/ all features
                 list(fun = ols,
                      assign_X = which(nonzero_X))) # ols w/ important features
  # Compute cross-validation instance
  cv_res <- crossval(y, X, Z = NULL,
                     learners,
                     cv_folds = 3,
                     silent = T)
  # Check output with expectations
  expect_equal(dim(cv_res$oos_resid), c(length(y), length(learners)))
})#TEST_THAT

test_that("crossval returns residuals by learner (w/ instruments)", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  Z <- matrix(rnorm(100*10), 100, 10)
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(100)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet), # lasso
                 list(fun = ols), # ols
                 list(fun = ols)) # ols again

  # Compute cross-validation instance
  cv_res <- crossval(D, X, Z,
                     learners,
                     cv_folds = 3,
                     silent = T)
  # Check output with expectations
  expect_equal(all(cv_res$oos_resid[,2] == cv_res$oos_resid[,3]), TRUE)
})#TEST_THAT
