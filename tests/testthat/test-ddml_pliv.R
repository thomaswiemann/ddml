test_that("ddml_pliv computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(what = mdl_glmnet,
                   args = list(alpha = 0.5))
  # Compute DDML PLIV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             sample_folds = 2,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_pliv computes with an ensemble procedure", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             ensemble_type = c("stacking"),
                             cv_folds = 2,
                             sample_folds = 2,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_pliv computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             ensemble_type = c("stacking", "stacking_nn",
                                               "stacking_01",
                                               "stacking_best", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 5)
})#TEST_THAT


test_that("ddml_pliv computes with different sets of learners", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols),
                   list(fun = ols))
  learners_ZX <- list(list(fun = mdl_glmnet,
                           args = list(alpha = 0.5)),
                      list(fun = mdl_glmnet,
                           args = list(alpha = 1)))
  learners_DX <- list(list(fun = ols),
                      list(fun = mdl_glmnet,
                           args = list(alpha = 0.5)),
                      list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             learners_ZX = learners_ZX,
                             learners_DX = learners_DX,
                             ensemble_type = c("stacking", "stacking_nn",
                                               "stacking_01",
                                               "stacking_best", "average"),
                             cv_folds = 2,
                             sample_folds = 2,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 5)
})#TEST_THAT