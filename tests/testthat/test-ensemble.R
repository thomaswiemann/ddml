test_that("ensemble_weights returns a weight matrix", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  Z <- matrix(rnorm(100*10), 100, 10)
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(100)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet),
                   list(fun = ols))
  # Compute ensemble weights with and without passthrough
  ensemble_types = c("average", "cv", "stacking", "stacking_nn", "stacking_10")
  ens_w_res <- ensemble_weights(D, X, Z,
                                type = ensemble_types,
                                learners,
                                cv_folds = 3,
                                silent = T)
  ens_w_res_pt <- ensemble_weights(D, X, Z,
                                   type = ensemble_types,
                                   learners,
                                   cv_folds = 3,
                                   cv_results = ens_w_res$cv_results,
                                   silent = T)
  # Check output with expectations
  expect_equal(dim(ens_w_res$weights), c(length(learners), 5))
  expect_equal(dim(ens_w_res_pt$weights), c(length(learners), 5))
})#TEST_THAT

test_that("ensemble returns a list of fitted learners", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  Z <- matrix(rnorm(100*10), 100, 10)
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(100)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet),
                   list(fun = ols))
  # Compute ensemble
  ens_fit <- ensemble(D, X, Z,
                      type = c("average", "stacking", "cv"),
                      learners,
                      cv_folds = 3,
                      silent = T)
  # Check output with expectations
  expect_equal(length(ens_fit$mdl_fits), length(learners))
})#TEST_THAT

test_that("prediction with ensemble learners returns fitted values", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  Z <- matrix(rnorm(100*10), 100, 10)
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(100)
  # Define arguments
  learners <- list(list(fun = rlasso,
                      args = list(iter_resid = 1, d = 5)),
                 list(fun = rlasso,
                      args = list(include = c(35:45),
                                  iter_resid = 1, d = 5)), # rlasso
                 list(fun = ols)) # ols w/ important features
  # Compute ensemble
  ens_fit <- ensemble(D, X, Z,
                      type = c("average", "stacking", "cv"),
                      learners,
                      cv_folds = 3,
                      silent = T)
  ens_fitted <- predict(ens_fit)
  # Check output with expectations
  expect_equal(dim(ens_fitted), c(length(D), 3))
})#TEST_THAT
