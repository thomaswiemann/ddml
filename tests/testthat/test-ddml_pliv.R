test_that("ddml_pliv computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(what = ols)
  # Compute DDML PLIV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_pliv computes with a single model and dependence", {
  # Simulate small dataset
  n_cluster <- 250
  nobs <- 500
  X <- cbind(1, matrix(rnorm(n_cluster*39), n_cluster, 39))
  Z_tld <-  X %*% runif(40) + rnorm(n_cluster)
  fun <- stepfun(quantile(Z_tld, probs = c(0.5)), c(0, 1))
  Z <- fun(Z_tld)
  cluster_variable <- sample(1:n_cluster, nobs, replace = TRUE)
  Z <- Z[cluster_variable, drop = F]
  X <- X[cluster_variable, , drop = F]
  eps <- rnorm(nobs)
  D <- Z + X %*% runif(40) + eps
  y <- D + X %*% runif(40) + 0.1 * eps + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  # Compute DDML PLIV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_pliv computes with an ensemble procedure", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             ensemble_type = "ols",
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_pliv computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "nnls1",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 5)
})#TEST_THAT


test_that("ddml_pliv computes with different sets of learners", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols),
                   list(fun = ols))
  learners_ZX <- list(list(fun = ols),
                      list(fun = ols))
  learners_DX <- list(list(fun = ols),
                      list(fun = ols),
                      list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             learners_ZX = learners_ZX,
                             learners_DX = learners_DX,
                             ensemble_type = c("ols", "nnls",
                                               "nnls1",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 5)
})#TEST_THAT

test_that("ddml_pliv computes with different sets of learners & shortstack", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols),
                   list(fun = ols))
  learners_ZX <- list(list(fun = ols),
                      list(fun = ols))
  learners_DX <- list(list(fun = ols),
                      list(fun = ols),
                      list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             learners_ZX = learners_ZX,
                             learners_DX = learners_DX,
                             ensemble_type = c("ols", "nnls",
                                               "nnls1",
                                               "singlebest", "average"),
                             shortstack = T,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 5)
})#TEST_THAT

test_that("summary.ddml_pliv computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(what = ols)
  # Compute DDML PLIV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             sample_folds = 3,
                             silent = T)
  inf_res <- summary(ddml_pliv_fit, type = "HC1")
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 1)
})#TEST_THAT


test_that("summary.ddml_pliv computes with custom ensemble weights", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute DDML PLIV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             sample_folds = 3,
                             custom_ensemble_weights = diag(1, 2),
                             silent = T)
  inf_res <- summary(ddml_pliv_fit, type = "HC1")
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 3)
})#TEST_THAT

test_that("ddml_pliv computes with a single model and multivariate D,Z", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  cbind(X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1], rnorm(nobs))
  Z <- cbind(Z, rnorm(nobs))
  y <- rowSums(D) + X %*% runif(40) + UV[, 2]

  # Define arguments
  learners <- list(what = ols)
  # Compute DDML PLIV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 2)
})#TEST_THAT

test_that("ddml_pliv computes with different ensembles and multivariate D,Z", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  cbind(X %*% runif(40) + Z %*% (1 + runif(1)) + UV[, 1], rnorm(nobs))
  Z <- cbind(Z, rnorm(nobs))
  y <- rowSums(D) + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "nnls1",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_pliv_fit$coef), 10)
})#TEST_THAT
