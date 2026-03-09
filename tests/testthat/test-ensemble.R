test_that("ensemble_weights returns a weight matrix", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  y <- X %*% runif(40) + rnorm(100)
  # Define arguments
  learners <- list(list(what = mdl_glmnet),
                   list(what = ols))
  # Compute ensemble weights with and without passthrough
  ensemble_types = c("average", "singlebest", "ols", "nnls", "nnls1")
  ens_w_res <- ensemble_weights(y, X,
                                type = ensemble_types,
                                learners,
                                cv_folds = 3,
                                silent = TRUE)
  ens_w_res_pt <- ensemble_weights(y, X,
                                   type = ensemble_types,
                                   learners,
                                   cv_folds = 3,
                                   cv_results = ens_w_res$cv_results,
                                   silent = TRUE)
  # Check output with expectations
  expect_equal(dim(ens_w_res$weights), c(length(learners), 5))
  expect_equal(dim(ens_w_res_pt$weights), c(length(learners), 5))
})#TEST_THAT

test_that("prediction with ensemble learners returns fitted values", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  y <- X %*% runif(40) + rnorm(100)
  # Define arguments
  learners <- list(list(what = mdl_glmnet),
                   list(what = ols))
  # Compute ensemble
  ens_fit <- ensemble(y, X,
                      type = c("average", "ols", "singlebest"),
                      learners,
                      cv_folds = 3,
                      silent = TRUE)
  ens_fitted <- predict(ens_fit, newdata = X)
  # Check output with expectations
  expect_equal(dim(ens_fitted), c(length(y), 3))
})#TEST_THAT

test_that("ensemble_weights returns a weight matrix w/ custom weights", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  y <- X %*% runif(40) + rnorm(100)
  # Define arguments
  learners <- list(list(what = mdl_glmnet),
                   list(what = ols))
  # Define custom weights
  weights_DX <- diag(length(learners))
  # Compute ensemble weights
  ensemble_types = c("average", "singlebest", "ols", "nnls", "nnls1")
  ens_w_res <- ensemble_weights(y, X,
                                type = ensemble_types,
                                learners,
                                cv_folds = 3,
                                custom_weights = weights_DX,
                                silent = TRUE)
  # Check output with expectations
  expect_equal(dim(ens_w_res$weights), c(length(learners), 7))
})#TEST_THAT

test_that("prediction w/ ensembles returns fitted values w/ custom weights", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  y <- X %*% runif(40) + rnorm(100)
  # Define arguments
  learners <- list(list(what = mdl_glmnet),
                   list(what = ols))
  # Define custom weights
  weights_DX <- diag(length(learners))
  # Compute ensemble
  ens_fit <- ensemble(y, X,
                      type = c("average", "ols", "singlebest"),
                      learners,
                      cv_folds = 3,
                      custom_weights = weights_DX,
                      silent = TRUE)
  ens_fitted <- predict(ens_fit, newdata = X)
  # Check output with expectations
  expect_equal(dim(ens_fitted), c(length(y), 5))
})#TEST_THAT

test_that("ensemble returns mean_y for constant outcomes", {
  # Simulate dataset with constant y
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  y <- rep(42, 100)  # Constant outcome

  # Define arguments
  learners <- list(list(what = mdl_glmnet),
                   list(what = ols))

  # Compute ensemble and expect warning
  expect_warning({
    ens_fit <- ensemble(y, X,
                        type = c("average", "ols", "singlebest"),
                        learners = learners,
                        cv_folds = 3,
                        silent = TRUE)
  }, "Outcome variable y is constant")

  # Check predictions
  ens_fitted <- predict.ensemble(ens_fit, newdata = X)

  # Check output matches expectations
  expect_equal(dim(ens_fitted), c(nrow(X), length(learners)))
  expect_true(all(ens_fitted == 42))
  expect_true(ens_fit$constant_y)
  expect_null(ens_fit$mdl_fits)
})
