library(ddml)
context("Testing ddml_tsls.")

sim_dat <- function(nobs) {
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- rnorm(nobs)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Organize and return output
  output <- list(y = y, D = D, Z = Z, X = X)
  return(output)
}#SIM_DAT

test_that("ddml_tsls computes with a single model", {
  # Simulate small dataset
  dat <- sim_dat(100)
  # Define arguments
  models <- list(what = mdl_glmnet,
                 args = list(alpha = 0.5))
  ddml_tsls_fit <- ddml_tsls(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         cv_folds = 3,
                         sample_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_tsls_fit$coef), 1)
})#TEST_THAT

test_that("ddml_tsls computes with an ensemble procedure", {
  # Simulate small dataset
  dat <- sim_dat(100)

  # Define arguments
  models <- list(list(fun = mdl_xgboost,
                      args = list(num_parallel_tree = 1,
                                  nrounds = 5)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_tsls_fit <- ddml_tsls(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         ens_type = c("stacking"),
                         cv_folds = 3,
                         sample_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_tsls_fit$coef), 1)
})#TEST_THAT


test_that("ddml_tsls computes with multiple ensemble procedures", {
  # Simulate small dataset
  dat <- sim_dat(100)

  # Define arguments
  models <- list(list(fun = mdl_xgboost,
                      args = list(num_parallel_tree = 3)),
                 list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_tsls_fit <- ddml_tsls(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         ens_type = c("stacking", "stacking_nn", "stacking_01",
                                      "cv", "average"),
                         cv_folds = 3,
                         sample_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_tsls_fit$coef), 5)
})#TEST_THAT


test_that("ddml_tsls computes with different sets of models", {
  # Simulate small dataset
  dat <- sim_dat(100)

  # Define arguments
  models <- list(list(fun = ols),
                 list(fun = ols),
                 list(fun = ols))
  models_ZX <- list(list(fun = mdl_xgboost,
                          args = list(num_parallel_tree = 3)),
                     list(fun = mdl_glmnet,
                          args = list(alpha = 0.5)),
                     list(fun = mdl_glmnet,
                          args = list(alpha = 1)))
  models_DX <- list(list(fun = ols),
                    list(fun = mdl_glmnet,
                         args = list(alpha = 0.5)),
                    list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_tsls_fit <- ddml_tsls(dat$y, dat$D, dat$Z, dat$X,
                         models,
                         models_ZX = models_ZX,
                         models_DX = models_DX,
                         ens_type = c("stacking", "stacking_nn", "stacking_01",
                                      "cv", "average"),
                         cv_folds = 3,
                         sample_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_tsls_fit$coef), 5)
})#TEST_THAT
