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
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_pliv_fit)), 2)
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
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_pliv_fit)), 2)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             ensemble_type = "ols",
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_pliv_fit)), 2)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "nnls1",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_pliv_fit)), 10)
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
  learners <- list(list(what = ols),
                   list(what = ols),
                   list(what = ols))
  learners_ZX <- list(list(what = ols),
                      list(what = ols))
  learners_DX <- list(list(what = ols),
                      list(what = ols),
                      list(what = ols))
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
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_pliv_fit)), 10)
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
  learners <- list(list(what = ols),
                   list(what = ols),
                   list(what = ols))
  learners_ZX <- list(list(what = ols),
                      list(what = ols))
  learners_DX <- list(list(what = ols),
                      list(what = ols),
                      list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             learners_ZX = learners_ZX,
                             learners_DX = learners_DX,
                             ensemble_type = c("ols", "nnls",
                                               "nnls1",
                                               "singlebest", "average"),
                             shortstack = TRUE,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_pliv_fit)), 10)
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
                             silent = TRUE)
  inf_res <- summary(ddml_pliv_fit, type = "HC1")
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_equal(length(coef(ddml_pliv_fit)), 2)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML PLIV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             sample_folds = 3,
                             custom_ensemble_weights = diag(1, 2),
                             silent = TRUE)
  inf_res <- summary(ddml_pliv_fit, type = "HC1")
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_equal(length(coef(ddml_pliv_fit)), 6)
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
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_pliv_fit)), 3)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_pliv_fit <- ddml_pliv(y, D, Z, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "nnls1",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_pliv_fit)), 15)
})#TEST_THAT

test_that("ddml_pliv HC0/HC1 SEs close to sandwich::vcovHC on iv_fit", {
  skip_if_not_installed("sandwich")
  skip_if_not_installed("AER")
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z <- matrix(rnorm(nobs), nobs, 1)
  UV <- matrix(rnorm(2 * nobs), nobs, 2) %*%
    chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <- X %*% c(1, 0.5, 0, 0, 0) +
    Z %*% 1.5 + UV[, 1]
  y <- 2 * D + X %*% c(0, 1, 0.5, 0.3, 0) + UV[, 2]

  fit <- ddml_pliv(y, D, Z, X,
                   learners = list(what = ols),
                   sample_folds = 5,
                   silent = TRUE)

  # Reconstruct the final partialing-out regression since iv_fit is removed:
  y_r <- as.vector(y - fit$fitted$y_X$cf_fitted_bylearner[, 1])
  D_r <- as.matrix(D - fit$fitted$D_X[[1]]$cf_fitted_bylearner[, 1])
  V_r <- as.matrix(Z - fit$fitted$Z_X[[1]]$cf_fitted_bylearner[, 1])
  colnames(D_r) <- colnames(D)
  colnames(V_r) <- colnames(Z)
  iv_fit <- AER::ivreg(y_r ~ D_r | V_r, x = TRUE)

  for (type in c("HC0", "HC1")) {
    V_ddml <- vcov(fit, type = type)
    V_sw <- sandwich::vcovHC(iv_fit, type = type)
    idx <- c(seq_len(nrow(V_sw))[-1], 1)
    V_sw_reord <- V_sw[idx, idx, drop = FALSE]
    expect_equal(as.numeric(V_ddml),
                 as.numeric(V_sw_reord),
                 tolerance = 1e-8,
                 info = paste("PLIV", type))
  }
  # HC3: regressor-based leverage matches sandwich exactly
  V_ddml <- vcov(fit, type = "HC3")
  V_sw <- sandwich::vcovHC(iv_fit, type = "HC3")
  idx <- c(seq_len(nrow(V_sw))[-1], 1)
  V_sw_reord <- V_sw[idx, idx, drop = FALSE]
  expect_equal(as.numeric(V_ddml),
               as.numeric(V_sw_reord),
               tolerance = 1e-8,
               info = "PLIV HC3")
})#TEST_THAT

test_that("ddml_pliv fitted pass-through works", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  Z <- matrix(rnorm(nobs * 2), nobs, 2)
  D <- X %*% c(1, 0.5, 0) + Z %*% c(0.5, 0.3) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_pliv(y, D, Z, X,
                   learners = learners,
                   ensemble_type = "average",
                   sample_folds = 2,
                   silent = TRUE)

  fit2 <- ddml_pliv(y, D, Z, X,
                    learners = learners,
                    ensemble_type = "average",
                    sample_folds = 2,
                    silent = TRUE,
                    fitted = fit$fitted,
                    splits = fit$splits)
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-6)

  expect_error(
    ddml_pliv(y, D, Z, X,
              learners = learners,
              ensemble_type = "average",
              sample_folds = 2,
              silent = TRUE,
              fitted = fit$fitted),
    "must be supplied when 'fitted' is supplied"
  )
})#TEST_THAT
