library(ddml)
context("Testing ensemble objects.")

test_that("ensemble_weights returns a weight matrix", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  Z <- matrix(rnorm(100*10), 100, 10)
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(100)
  # Define arguments
  models <- list(list(fun = rlasso,
                      args = list(include = NULL,
                                  iter_resid = 1, d = 5),
                      assign_X = c(1:ncol(X)),
                      assign_Z = c(1:ncol(Z))),
                 list(fun = rlasso,
                      args = list(include = c(35:45),
                                  iter_resid = 1, d = 5),
                      assign_X = c(1:ncol(X)),
                      assign_Z = c(1:ncol(Z))), # rlasso
                 list(fun = ols,
                      args = list(),
                      assign_X = c(1:ncol(X)),
                      assign_Z = c(1:ncol(Z)))) # ols w/ important features
  # Compute ensemble weights
  ens_w_res <- ensemble_weights(D, X, Z,
                                type = c("average", "stacking", "cv"),
                                models,
                                cv_folds = 3,
                                silent = T)
  # Compute via pass-through of cv_res
  ens_w_res_pt <- ensemble_weights(D, X, Z,
                                   type = c("average", "stacking", "cv"),
                                   models,
                                   cv_res = ens_w_res$cv_res,
                                   silent = T)
  # Check output with expectations
  expect_equal(colSums(ens_w_res$weights), rep(1, 3))
  expect_equal(colSums(ens_w_res_pt$weights), rep(1, 3))
})#TEST_THAT

test_that("ensemble returns a list of fitted models", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  Z <- matrix(rnorm(100*10), 100, 10)
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(100)
  # Define arguments
  models <- list(list(fun = rlasso,
                      args = list(include = NULL,
                                  iter_resid = 1, d = 5),
                      assign_X = c(1:ncol(X)),
                      assign_Z = c(1:ncol(Z))),
                 list(fun = rlasso,
                      args = list(include = c(35:45),
                                  iter_resid = 1, d = 5),
                      assign_X = c(1:ncol(X)),
                      assign_Z = c(1:ncol(Z))), # rlasso
                 list(fun = ols,
                      args = list(),
                      assign_X = c(1:ncol(X)),
                      assign_Z = c(1:ncol(Z)))) # ols w/ important features
  # Compute ensemble
  ens_fit <- ensemble(D, X, Z,
                      type = c("average", "stacking", "cv"),
                      models,
                      cv_folds = 3,
                      silent = T)
  # Check output with expectations
  expect_equal(length(ens_fit$mdl_fits), length(models))
  expect_is(ens_fit$mdl_w_iv, "integer")
})#TEST_THAT

test_that("prediction with ensemble models returns fitted values", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  Z <- matrix(rnorm(100*10), 100, 10)
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(100)
  # Define arguments
  models <- list(list(fun = rlasso,
                      args = list(include = NULL,
                                  iter_resid = 1, d = 5),
                      assign_X = c(1:ncol(X)),
                      assign_Z = c(1:ncol(Z))),
                 list(fun = rlasso,
                      args = list(include = c(35:45),
                                  iter_resid = 1, d = 5),
                      assign_X = c(1:ncol(X)),
                      assign_Z = c(1:ncol(Z))), # rlasso
                 list(fun = ols,
                      args = list(),
                      assign_X = c(1:ncol(X)),
                      assign_Z = c(1:ncol(Z)))) # ols w/ important features
  # Compute ensemble
  ens_fit <- ensemble(D, X, Z,
                      type = c("average", "stacking", "cv"),
                      models,
                      cv_folds = 3,
                      silent = T)
  ens_fitted <- predict(ens_fit)
  # Check output with expectations
  expect_equal(dim(ens_fitted), c(length(D), 3))
})#TEST_THAT
