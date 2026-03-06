test_that("ddml_late computes with a single model", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z_tld <- 0.1 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.2 * Z + 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_late_fit)), 1)
})#TEST_THAT

test_that("ddml_late computes with stratify = FALSE", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z_tld <- 0.1 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.2 * Z + 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               stratify = FALSE,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_late_fit)), 1)
})#TEST_THAT

test_that("ddml_late computes with a single model and dependence", {
  # Simulate small dataset
  n_cluster <- 250
  nobs <- 500
  X <- matrix(rnorm(n_cluster * 5), n_cluster, 5)
  Z_tld <- 0.1 * X[, 1] + rnorm(n_cluster)
  Z <- 1 * (Z_tld > 0)
  cluster_variable <- sample(seq_len(n_cluster), nobs,
                             replace = TRUE)
  Z <- Z[cluster_variable]
  X <- X[cluster_variable, , drop = FALSE]
  eps <- rnorm(nobs)
  D <- Z + 0.1 * X[, 1] + eps
  y <- D + 0.1 * X[, 1] + 0.1 * eps + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               cluster_variable = cluster_variable,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_late_fit)), 1)
})#TEST_THAT

test_that("ddml_late computes with a single model & perfect non-compliance", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z_tld <- 0.1 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.5 * Z + 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  D[Z == 0] <- 0
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_late_fit)), 1)
})#TEST_THAT

test_that("ddml_late computes with an ensemble procedure", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z_tld <- 0.1 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.2 * Z + 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               ensemble_type = "ols",
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_late_fit)), 1)
})#TEST_THAT

test_that("ddml_late computes w/ multiple ensembles & custom weights", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z_tld <- 0.1 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.2 * Z + 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "nnls1",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               custom_ensemble_weights = diag(1, 2),
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_late_fit)), 7)
})#TEST_THAT

test_that("ddml_late computes with multiple ensemble procedures + perfect compliance", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z_tld <- 0.1 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.2 * Z + 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  D[Z == 1] <- 1 # perfect compliance
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "nnls1",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_late_fit)), 5)
})#TEST_THAT

test_that("ddml_late computes w/ mult ensembles, custom weights, & shortstack", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z_tld <- 0.1 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.2 * Z + 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "nnls1",
                                                 "singlebest", "average"),
                               shortstack = TRUE,
                               cv_folds = 3,
                               custom_ensemble_weights = diag(1, 2),
                               sample_folds = 3,
                               silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_late_fit)), 7)
})#TEST_THAT

test_that("summary.ddml_late computes with a single model", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  Z_tld <- 0.1 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.2 * Z + 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Compute inference results & test print
  inf_res <- summary(ddml_late_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(1, 4, 1))
})#TEST_THAT

test_that("summary.ddml_late computes with a single model and dependence", {
  # Simulate small dataset
  n_cluster <- 250
  nobs <- 500
  X <- matrix(rnorm(n_cluster * 5), n_cluster, 5)
  Z_tld <- 0.1 * X[, 1] + rnorm(n_cluster)
  Z <- 1 * (Z_tld > 0)
  cluster_variable <- sample(seq_len(n_cluster), nobs,
                             replace = TRUE)
  Z <- Z[cluster_variable]
  X <- X[cluster_variable, , drop = FALSE]
  eps <- rnorm(nobs)
  D <- Z + 0.1 * X[, 1] + eps
  y <- D + 0.1 * X[, 1] + 0.1 * eps + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               cluster_variable = cluster_variable,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = TRUE)
  # Compute inference results & test print
  inf_res <- summary(ddml_late_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(1, 4, 1))
})#TEST_THAT

test_that("ddml_late fitted pass-through works", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  Z_tld <- 0.1 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.2 * Z + 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_late(y, D, Z, X,
                   learners = learners,
                   ensemble_type = "average",
                   sample_folds = 2,
                   silent = TRUE)

  fit2 <- ddml_late(y, D, Z, X,
                    learners = learners,
                    ensemble_type = "average",
                    sample_folds = 2,
                    silent = TRUE,
                    fitted = fit$fitted,
                    splits = fit$splits)
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-6)

  expect_error(
    ddml_late(y, D, Z, X,
              learners = learners,
              ensemble_type = "average",
              sample_folds = 2,
              silent = TRUE,
              fitted = fit$fitted),
    "must be supplied when 'fitted' is supplied"
  )
})#TEST_THAT

test_that("ddml_late scores are mean-zero", {
  nobs <- 500
  set.seed(42)
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  Z_tld <- 0.2 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.3 * Z + 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  fit <- ddml_late(y, D, Z, X,
                   learners = list(what = ols),
                   sample_folds = 3, silent = TRUE)
  score_mean <- mean(fit$scores[[1]])
  expect_true(abs(score_mean) < 0.05,
              label = paste("score mean =", round(score_mean, 6)))
})

test_that("ddml_late has correct sign (matches Wald estimator)", {
  nobs <- 2000
  set.seed(123)
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  Z_tld <- 0.5 * X[, 1] + rnorm(nobs)
  Z <- 1 * (Z_tld > 0)
  D_tld <- 0.6 * Z + 0.2 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- 1.0 * D + 0.3 * X[, 1] + rnorm(nobs)
  fit <- ddml_late(y, D, Z, X,
                   learners = list(what = ols),
                   sample_folds = 5, silent = TRUE)
  # Wald estimator: cov(y,Z)/cov(D,Z)
  wald <- as.numeric(stats::cov(y, Z) / stats::cov(D, Z))
  late_hat <- as.numeric(coef(fit))
  # LATE should have the same sign as Wald
  expect_true(sign(late_hat) == sign(wald))
  # And be in a reasonable range
  expect_equal(late_hat, wald, tolerance = 0.5)
})
