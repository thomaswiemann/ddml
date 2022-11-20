library(ddml)
context("Testing crosspred.")

sim_dat <- function(nobs) {
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10) # overidentified
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(nobs)
  # Organize and return output
  output <- list(D = D, Z = Z, X = X)
  return(output)
}#SIM_DAT

test_that("crosspred computes with a single model", {
  # Simulate small dataset
  dat <- sim_dat(100)
  # Define arguments
  learners <- list(what = mdl_xgboost)
  # Compute cross-sample predictions
  y = dat$D; X = dat$X;  Z = dat$Z
  crosspred_res <- crosspred(y, X, Z,
                             learners,
                             sample_folds = 3,
                             compute_insample_predictions = T,
                             silent = T)
  # Check output with expectations
  expect_equal(length(crosspred_res$oos_fitted), length(dat$D))
  expect_equal(length(crosspred_res$is_fitted), 3)
})#TEST_THAT

test_that("crosspred computes with ensemble procedures", {
  # Simulate small dataset
  dat <- sim_dat(100)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet),
                 list(fun = mdl_randomForest),
                 list(fun = ols))
  # Compute cross-sample predictions
  crosspred_res <- crosspred(dat$D, dat$X, dat$Z,
                             learners,
                             ensemble_type = c("average", "stacking",
                                          "stacking_01", "stacking_nn",
                                          "cv"),
                             cv_folds = 3,
                             sample_folds = 3,
                             compute_insample_predictions = T,
                             silent = T)
  # Check output with expectations
  expect_equal(dim(crosspred_res$oos_fitted), c(length(dat$D), 5))
  expect_equal(length(crosspred_res$is_fitted), 5)
})#TEST_THAT

test_that("crosspred computes with ensemble procedures and sparse matrices", {
  # Simulate small dataset
  dat <- sim_dat(100)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet),
                 list(fun = ols))
  # Compute cross-sample predictions
  crosspred_res <- crosspred(dat$D, as(dat$X, "sparseMatrix"),
                             as(dat$Z, "sparseMatrix"),
                             learners,
                             ensemble_type = c("average", "stacking",
                                          "stacking_01", "stacking_nn",
                                          "cv"),
                             cv_folds = 3,
                             sample_folds = 3,
                             compute_insample_predictions = T,
                             silent = T)
  # Check output with expectations
  expect_equal(dim(crosspred_res$oos_fitted), c(length(dat$D), 5))
  expect_equal(length(crosspred_res$is_fitted), 5)
})#TEST_THAT
