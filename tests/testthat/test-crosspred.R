sim_dat <- function(nobs) {
  # generate test data
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10) # overidentified
  y <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(nobs)
  # Organize and return output
  output <- list(D = D, Z = Z, X = X)
  return(output)
}#SIM_DAT

test_that("crosspred computes with a single model", {
  # generate test data
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10) # overidentified
  y <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(nobs)
  # Define arguments
  learners <- list(what = mdl_xgboost)
  # Compute cross-sample predictions
  crosspred_res <- crosspred(y, X, Z,
                             learners,
                             sample_folds = 3,
                             compute_insample_predictions = T,
                             silent = T)
  # Check output with expectations
  expect_equal(length(crosspred_res$oos_fitted), length(y))
  expect_equal(length(crosspred_res$is_fitted), 3)
})#TEST_THAT

test_that("crosspred computes with ensemble procedures", {
  # generate test data
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10) # overidentified
  y <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet),
                 list(fun = mdl_randomForest),
                 list(fun = ols))
  # Compute cross-sample predictions
  crosspred_res <- crosspred(y, X, Z,
                             learners,
                             ensemble_type = c("average", "stacking",
                                          "stacking_01", "stacking_nn",
                                          "stacking_best"),
                             cv_folds = 3,
                             sample_folds = 3,
                             compute_insample_predictions = T,
                             silent = T)
  # Check output with expectations
  expect_equal(dim(crosspred_res$oos_fitted), c(length(y), 5))
  expect_equal(length(crosspred_res$is_fitted), 5)
})#TEST_THAT

test_that("crosspred computes with ensemble procedures and sparse matrices", {
  # generate test data
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10) # overidentified
  y <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet),
                 list(fun = ols))
  # Compute cross-sample predictions
  crosspred_res <- crosspred(y, as(X, "sparseMatrix"),
                             as(Z, "sparseMatrix"),
                             learners,
                             ensemble_type = c("average", "stacking",
                                          "stacking_01", "stacking_nn",
                                          "stacking_best"),
                             cv_folds = 3,
                             sample_folds = 3,
                             compute_insample_predictions = T,
                             silent = T)
  # Check output with expectations
  expect_equal(dim(crosspred_res$oos_fitted), c(length(y), 5))
  expect_equal(length(crosspred_res$is_fitted), 5)
})#TEST_THAT

test_that("crosspred computes auxilliary predictions", {
  # generate test data
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  y <-  X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet),
                   list(fun = mdl_randomForest),
                   list(fun = ols))
  # Compute cross-sample and auxilliary predictions
  crosspred_res <- crosspred(y, X,
                             learners = learners,
                             ensemble_type = c("average", "stacking",
                                               "stacking_01", "stacking_nn",
                                               "stacking_best"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T,
                             auxilliary_X = list(X, X, X))
  # Check output with expectations
  expect_equal(dim(crosspred_res$auxilliary_fitted[[1]]), c(length(y), 5))
})#TEST_THAT
