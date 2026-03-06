test_that("ddml_fpliv computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10)
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Define arguments
  learners <- list(list(what = ols))
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 2)
})#TEST_THAT

test_that("ddml_fpliv computes with a single model and dependence", {
  # Simulate small dataset
  n_cluster <- 250
  nobs <- 500
  X <- cbind(1, matrix(rnorm(n_cluster*39), n_cluster, 39))
  Z_tld <-  X %*% runif(40) + rnorm(n_cluster)
  fun <- stepfun(quantile(Z_tld, probs = 0.5), c(0, 1))
  Z <- fun(Z_tld)
  cluster_variable <- sample(1:n_cluster, nobs, replace = TRUE)
  Z <- as.matrix(Z[cluster_variable, drop = F])
  X <- X[cluster_variable, , drop = F]
  eps <- rnorm(nobs)
  D <- Z + X %*% runif(40) + eps
  y <- D + X %*% runif(40) + 0.1 * eps + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols))
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               cluster_variable = cluster_variable,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 2)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = "ols",
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 2)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = "ols",
                               sample_folds = 3,
                               cv_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 2)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)

  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 8)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("average"),
                               cv_folds = 3,
                               custom_ensemble_weights = diag(1, 2),
                               sample_folds = 3,
                               silent = TRUE)

  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 6)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 8)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D,
                               Z = as(Z, "sparseMatrix"),
                               X = as(X, "sparseMatrix"),
                               learners = learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 8)
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
  learners <- list(list(what = ols),
                   list(what = ols),
                   list(what = ols))
  learners_DXZ <- list(list(what = ols),
                       list(what = ols))
  learners_DX <- list(list(what = ols),
                      list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               learners_DXZ = learners_DXZ,
                               learners_DX = learners_DX,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)

  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 8)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               shortstack = TRUE,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)

  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 8)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 8)
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
                               silent = TRUE)
  inf_res <- summary(ddml_fpliv_fit)
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(2, 4, 1))
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = "ols",
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 3)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = "ols",
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 3)
})#TEST_THAT

test_that("ddml_fpliv computes with multiple ensemble procedures, multi D", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z <- matrix(rnorm(nobs * 3), nobs, 3)
  UV <- matrix(rnorm(2 * nobs), nobs, 2) %*%
    chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <- cbind(0.3 * X[, 1] + Z %*% c(0.5, 0.3, 0.1) + UV[, 1],
             0.2 * X[, 2] + Z %*% c(0.2, 0.4, 0.1) + rnorm(nobs))
  y <- rowSums(D) + 0.3 * X[, 1] + UV[, 2]
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                                 learners,
                                 ensemble_type = c("ols", "nnls",
                                                   "singlebest", "average"),
                                 cv_folds = 3,
                                 sample_folds = 5,
                                 silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 12)
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
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML IV estimator
  ddml_fpliv_fit <- ddml_fpliv(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_fpliv_fit)), 12)
})#TEST_THAT

test_that("ddml_fpliv HC0/HC1 SEs close to sandwich::vcovHC on iv_fit", {
  skip_if_not_installed("sandwich")
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z <- matrix(rnorm(nobs * 2), nobs, 2)
  UV <- matrix(rnorm(2 * nobs), nobs, 2) %*%
    chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <- X %*% c(1, 0.5, 0, 0, 0) +
    Z %*% c(0.5, 0.3) + UV[, 1]
  y <- 2 * D + X %*% c(0, 1, 0.5, 0.3, 0) + UV[, 2]

  fit <- ddml_fpliv(y, D, Z, X,
                    learners = list(what = ols),
                    sample_folds = 5,
                    silent = TRUE)

  for (type in c("HC0", "HC1")) {
    V_ddml <- vcov(fit, type = type)
    V_sw <- sandwich::vcovHC(fit$iv_fit[[1]], type = type)
    idx <- c(seq_len(nrow(V_sw))[-1], 1)
    V_sw_reord <- V_sw[idx, idx, drop = FALSE]
    expect_equal(as.numeric(V_ddml),
                 as.numeric(V_sw_reord),
                 tolerance = 1e-8,
                 info = paste("FPLIV", type))
  }
  # HC3: regressor-based leverage matches sandwich exactly
  V_ddml <- vcov(fit, type = "HC3")
  V_sw <- sandwich::vcovHC(fit$iv_fit[[1]], type = "HC3")
  idx <- c(seq_len(nrow(V_sw))[-1], 1)
  V_sw_reord <- V_sw[idx, idx, drop = FALSE]
  expect_equal(as.numeric(V_ddml),
               as.numeric(V_sw_reord),
               tolerance = 1e-8,
               info = "FPLIV HC3")
})#TEST_THAT

test_that("ddml_fpliv fitted pass-through works", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  Z <- matrix(rnorm(nobs * 2), nobs, 2)
  D <- X %*% c(1, 0.5, 0) + Z %*% c(0.5, 0.3) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_fpliv(y, D, Z, X,
                    learners = learners,
                    ensemble_type = "average",
                    sample_folds = 2,
                    silent = TRUE)
  fit2 <- ddml_fpliv(y, D, Z, X,
                     learners = learners,
                     ensemble_type = "average",
                     sample_folds = 2,
                     silent = TRUE,
                     fitted = fit$fitted,
                     splits = fit$splits)
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-6)

  # Error when fitted supplied without splits
  expect_error(
    ddml_fpliv(y, D, Z, X,
               learners = learners,
               sample_folds = 2,
               silent = TRUE,
               fitted = fit$fitted),
    "must be supplied when 'fitted' is supplied"
  )
})#TEST_THAT
