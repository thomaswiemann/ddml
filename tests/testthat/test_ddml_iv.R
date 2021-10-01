library(ddml)
context("Testing ddml_iv.")

sim_dat <- function(nobs) {
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10) # overidentified
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Organize and return output
  output <- list(y = y, D = D, Z = Z, X = X)
  return(output)
}#SIM_DAT

test_that("ddml_iv computes with a single model", {
  # Simulate small dataset
  dat <- sim_dat(100)
  # Define arguments
  models <- list(what = mdl_glmnet,
                 args = list(alpha = 0.5))
  ddml_iv_fit <- ddml_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         cv_folds = 3,
                         sample_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_iv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_iv computes with an ensemble procedure", {
  # Simulate small dataset
  dat <- sim_dat(100)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  models <- list(list(fun = mdl_xgboost,
                      args = list(num_parallel_tree = 3)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_iv_fit <- ddml_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         ens_type = c("stacking"),
                         cv_folds = 3,
                         sample_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_iv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_iv computes with an ensemble procedure w/o enforcing the LIE", {
  # Simulate small dataset
  dat <- sim_dat(100)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  models <- list(list(fun = mdl_xgboost,
                      args = list(num_parallel_tree = 3)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_iv_fit <- ddml_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         ens_type = c("stacking"),
                         cv_folds = 3,
                         sample_folds = 3,
                         enforce_LIE = FALSE,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_iv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_iv computes with multiple ensemble procedures", {
  # Simulate small dataset
  dat <- sim_dat(100)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  models <- list(list(fun = mdl_xgboost,
                      args = list(num_parallel_tree = 3)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_iv_fit <- ddml_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         ens_type = c("stacking", "stacking_nn", "stacking_01",
                                      "cv", "average"),
                         cv_folds = 3,
                         sample_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_iv_fit$coef), 5)
})#TEST_THAT

test_that("ddml_iv computes with multiple ensemble procedures w/o the LIE", {
  # Simulate small dataset
  dat <- sim_dat(100)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  models <- list(list(fun = mdl_xgboost,
                      args = list(num_parallel_tree = 3)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_iv_fit <- ddml_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         ens_type = c("stacking", "stacking_nn", "stacking_01",
                                      "cv", "average"),
                         cv_folds = 3,
                         sample_folds = 3,
                         enforce_LIE = FALSE,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_iv_fit$coef), 5)
})#TEST_THAT

test_that("ddml_iv computes with multiple ensemble procedures and
          sparse matrices", {
  # Simulate small dataset
  dat <- sim_dat(100)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  models <- list(list(fun = mdl_xgboost,
                      args = list(num_parallel_tree = 3)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_iv_fit <- ddml_iv(dat$y, dat$D, as(dat$Z, "sparseMatrix"),
                         as(dat$X, "sparseMatrix"),
                         models,
                         ens_type = c("stacking", "stacking_nn", "stacking_01",
                                      "cv", "average"),
                         cv_folds = 3,
                         sample_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_iv_fit$coef), 5)
})#TEST_THAT

test_that("ddml_iv computes with two sets of models", {
  # Simulate small dataset
  dat <- sim_dat(100)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  models <- list(list(fun = ols),
                 list(fun = ols),
                 list(fun = ols))
  models_FS <- list(list(fun = mdl_xgboost,
                      args = list(num_parallel_tree = 3)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 1)))
  # Compute LIE-conform DDML IV estimator
  ddml_iv_fit <- ddml_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         models_FS = models_FS,
                         ens_type = c("stacking", "stacking_nn", "stacking_01",
                                      "cv", "average"),
                         cv_folds = 3,
                         sample_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_iv_fit$coef), 5)
})#TEST_THAT
