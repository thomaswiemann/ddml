library(ddml)
context("Testing pl_iv.")

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

test_that("pl_iv computes with a single model", {
  # Simulate small dataset
  dat <- sim_dat(100)
  # Define arguments
  models <- list(what = mdl_glmnet,
                 args = list(alpha = 0.5))
  y = dat$y; D = dat$D; Z = dat$Z; X = dat$X
  # Estimate IV coefficient from partially linear model
  pl_iv_fit <- pl_iv(y, D, Z, X,
                     models,
                     cv_folds = 3,
                     silent = T)
  # Check output with expectations
  expect_equal(length(pl_iv_fit$coef), 1)
})#TEST_THAT

test_that("pl_iv computes with an ensemble procedure", {
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
  # Estimate IV coefficient from partially linear model
  pl_iv_fit <- pl_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         ens_type = c("stacking"),
                         cv_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(pl_iv_fit$coef), 1)
})#TEST_THAT

test_that("pl_iv computes with an ensemble procedure w/o enforcing the LIE", {
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
  # Estimate IV coefficient from partially linear model
  pl_iv_fit <- pl_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         ens_type = c("stacking"),
                         cv_folds = 3,
                         enforce_LIE = FALSE,
                         silent = T)
  # Check output with expectations
  expect_equal(length(pl_iv_fit$coef), 1)
})#TEST_THAT

test_that("pl_iv computes with multiple ensemble procedures", {
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
  # Estimate IV coefficient from partially linear model
  pl_iv_fit <- pl_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         ens_type = c("stacking", "stacking_nn", "stacking_01",
                                      "cv", "average"),
                         cv_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(pl_iv_fit$coef), 5)
})#TEST_THAT

test_that("pl_iv computes with multiple ensemble procedures w/o the LIE", {
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
  # Estimate IV coefficient from partially linear model
  pl_iv_fit <- pl_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         ens_type = c("stacking", "stacking_nn", "stacking_01",
                                      "cv", "average"),
                         cv_folds = 3,
                         enforce_LIE = FALSE,
                         silent = T)

  # Check output with expectations
  expect_equal(length(pl_iv_fit$coef), 5)
})#TEST_THAT

test_that("pl_iv computes with multiple ensemble procedures and
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
  # Estimate IV coefficient from partially linear model
  pl_iv_fit <- pl_iv(dat$y, dat$D, as(dat$Z, "sparseMatrix"),
                         as(dat$X, "sparseMatrix"),
                         models,
                         ens_type = c("stacking", "stacking_nn", "stacking_01",
                                      "cv", "average"),
                         cv_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(pl_iv_fit$coef), 5)
})#TEST_THAT

test_that("pl_iv computes with two sets of models", {
  # Simulate small dataset
  dat <- sim_dat(100)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  models <- list(list(fun = ols),
                 list(fun = ols),
                 list(fun = ols))
  models_DXZ <- list(list(fun = mdl_xgboost,
                      args = list(num_parallel_tree = 3)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 1)))
  # Estimate IV coefficient from partially linear model
  pl_iv_fit <- pl_iv(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         models_DXZ = models_DXZ,
                         ens_type = c("stacking", "stacking_nn", "stacking_01",
                                      "cv", "average"),
                         cv_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(pl_iv_fit$coef), 5)
})#TEST_THAT

test_that("pl_iv computes with constant X", {
  # Simulate small dataset
  dat <- sim_dat(100)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  models <- list(list(fun = ols),
                 list(fun = ols),
                 list(fun = ols))
  models_DXZ <- list(list(fun = mdl_xgboost,
                          args = list(num_parallel_tree = 3)),
                     list(fun = mdl_glmnet,
                          args = list(alpha = 0.5)),
                     list(fun = mdl_glmnet,
                          args = list(alpha = 1)))
  # Estimate IV coefficient from partially linear model
  pl_iv_fit <- pl_iv(dat$y, dat$D, dat$Z, matrix(1, 100, 1),
                     models = models,
                     models_DXZ = models_DXZ,
                     models_DX = models,
                     ens_type = c("stacking", "stacking_nn", "stacking_01",
                                  "cv", "average"),
                     cv_folds = 3,
                     silent = T)
  # Check output with expectations
  expect_equal(length(pl_iv_fit$coef), 5)
})#TEST_THAT
