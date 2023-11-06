test_that("ddml_fpliv computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols))
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               sample_folds = 3,
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
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = "ols",
                               cv_folds = 3,
                               sample_folds = 3,
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
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = "ols",
                               sample_folds = 3,
                               cv_folds = 3,
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
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)

  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 4)
})#TEST_THAT

test_that("ddml_fpliv computes with custom weights", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("average"),
                               cv_folds = 3,
                               custom_ensemble_weights = diag(1, 2),
                               sample_folds = 3,
                               silent = T)

  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 3)
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
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               enforce_LIE = FALSE,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 4)
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
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D,
                               Z = as(Z, "sparseMatrix"),
                               X = as(X, "sparseMatrix"),
                               learners = learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 4)
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
  learners_DXZ <- list(list(fun = ols),
                       list(fun = ols))
  learners_DX <- list(list(fun = ols),
                      list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               learners_DXZ = learners_DXZ,
                               learners_DX = learners_DX,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)

  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 4)
})#TEST_THAT

test_that("ddml_fpliv computes w/ ensembles & shortstack", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               shortstack = T,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)

  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 4)
})#TEST_THAT

test_that("ddml_fpliv computes w/ ensembles & shortstack but w/o the LIE ", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               enforce_LIE = FALSE,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 4)
})#TEST_THAT

test_that("summary.ddml_fpliv computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(what = ols)
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               sample_folds = 3,
                               silent = T)
  capture_output({inf_res <- summary(ddml_fpliv_fit, type = "HC1")})
  # Check output with expectations
  expect_equal(length(inf_res), 8)
})#TEST_THAT

test_that("ddml_fpliv computes with an ensemble procedure, multi D", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  cbind(X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1], rnorm(nobs))
  y <- rowSums(D) + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = "ols",
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 2)
})#TEST_THAT

test_that("ddml_fpliv computes with an ensemble procedure w/o LIE, multi D", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  cbind(X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1], rnorm(nobs))
  y <- rowSums(D) + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = "ols",
                               cv_folds = 3,
                               sample_folds = 3,
                               enforce_LIE = F,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 2)
})#TEST_THAT

test_that("ddml_fpliv computes with multiple ensemble procedures, multi D", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  cbind(X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1], rnorm(nobs))
  y <- rowSums(D) + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 8)
  })#TEST_THAT

test_that("ddml_fpliv computes with ensemble procedures w/o LIE, multi D", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  cbind(X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1], rnorm(nobs))
  y <- rowSums(D) + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               enforce_LIE = F,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_fpliv_fit$coef), 8)
})#TEST_THAT
