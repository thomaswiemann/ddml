test_that("ddml_fpliv computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(what = mdl_glmnet,
                   args = list(alpha = 0.5))
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               sample_folds = 2,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_fpliv computes with an ensemble procedure", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols"),
                               cv_folds = 2,
                               sample_folds = 2,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_fpliv computes with stacking w/o enforcing the LIE", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols"),
                               sample_folds = 2,
                               cv_folds = 2,
                               enforce_LIE = FALSE,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_fpliv computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "nnls1",
                                                 "singlebest", "average"),
                               cv_folds = 2,
                               sample_folds = 2,
                               silent = T)

  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 5)
})#TEST_THAT

test_that("ddml_fpliv computes with multiple ensembles w/o the LIE", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "nnls1",
                                                 "singlebest", "average"),
                               cv_folds = 2,
                               sample_folds = 2,
                               enforce_LIE = FALSE,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 5)
})#TEST_THAT

test_that("ddml_fpliv computes with multiple ensembles and sparse matrices", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D,
                               as(Z, "sparseMatrix"),
                               as(X, "sparseMatrix"),
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "nnls1",
                                                 "singlebest", "average"),
                               cv_folds = 2,
                               sample_folds = 2,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 5)
})#TEST_THAT

test_that("ddml_fpliv computes with different sets of learners", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols),
                   list(fun = ols))
  learners_DXZ <- list(list(fun = mdl_glmnet,
                            args = list(alpha = 0.5)),
                       list(fun = mdl_glmnet,
                            args = list(alpha = 1)))
  learners_DX <- list(list(fun = ols),
                      list(fun = mdl_glmnet),
                      list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               learners_DXZ = learners_DXZ,
                               learners_DX = learners_DX,
                               ensemble_type = c("ols", "nnls",
                                                 "nnls1",
                                                 "singlebest", "average"),
                               cv_folds = 2,
                               sample_folds = 2,
                               silent = T)

  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 5)
})#TEST_THAT
