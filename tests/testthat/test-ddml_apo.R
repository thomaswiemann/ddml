test_that("ddml_apo computes with a single model", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_apo_fit <- ddml_apo(y, D, X,
                             d = 1,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_apo_fit)), 1)
})#TEST_THAT

test_that("ddml_apo computes with weights", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  weights <- runif(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_apo_fit <- ddml_apo(y, D, X,
                             d = 0,
                             weights = weights,
                             learners = learners,
                             stratify = FALSE,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_apo_fit)), 1)
})#TEST_THAT

test_that("ddml_apo computes with a single model and dependence", {
  # Simulate small dataset
  n_cluster <- 200
  nobs <- 500
  X <- matrix(rnorm(n_cluster * 5), n_cluster, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(n_cluster)
  D <- 1 * (D_tld > 0)
  cluster_variable <- sample(seq_len(n_cluster), nobs,
                             replace = TRUE)
  D <- D[cluster_variable]
  X <- X[cluster_variable, , drop = FALSE]
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_apo_fit <- ddml_apo(y, D, X,
                             d = 1,
                             learners = learners,
                             cluster_variable = cluster_variable,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_apo_fit)), 1)
})#TEST_THAT

test_that("ddml_apo computes with an ensemble procedure", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_apo_fit <- ddml_apo(y, D, X,
                             d = 1,
                             learners = learners,
                             ensemble_type = "ols",
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_apo_fit)), 1)
})#TEST_THAT

test_that("ddml_apo computes w/ multiple ensembles + custom weights", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_apo_fit <- ddml_apo(y, D, X,
                             d = 0,
                             weights = rep(1, nobs),
                             learners = learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             custom_ensemble_weights = diag(1, 2),
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_apo_fit)), 6)
})#TEST_THAT

test_that("ddml_apo computes with multiple ensemble procedures & shortstack", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols))
  # Compute DDML PLM estimator
  ddml_apo_fit <- ddml_apo(y, D, X,
                             d = 1,
                             learners = learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             shortstack = TRUE,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_apo_fit)), 4)
})#TEST_THAT

test_that("summary.ddml_apo computes with a single model", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_apo_fit <- ddml_apo(y, D, X,
                             d = 1,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Compute inference results & test print
  inf_res <- summary(ddml_apo_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(1, 4, 1))
})#TEST_THAT

test_that("summary.ddml_apo computes with a single model and dependence", {
  # Simulate small dataset
  n_cluster <- 200
  nobs <- 500
  X <- matrix(rnorm(n_cluster * 5), n_cluster, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(n_cluster)
  D <- 1 * (D_tld > 0)
  cluster_variable <- sample(seq_len(n_cluster), nobs,
                             replace = TRUE)
  D <- D[cluster_variable]
  X <- X[cluster_variable, , drop = FALSE]
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_apo_fit <- ddml_apo(y, D, X,
                             d = 1,
                             learners = learners,
                             cluster_variable = cluster_variable,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Compute inference results & test print
  inf_res <- summary(ddml_apo_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(1, 4, 1))
})#TEST_THAT

test_that("summary.ddml_apo computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols))
  # Compute DDML PLM estimator
  ddml_apo_fit <- ddml_apo(y, D, X,
                             d = 1,
                             learners = learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Compute inference results & test print
  inf_res <- summary(ddml_apo_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(1, 4, 4))
})#TEST_THAT

test_that("ddml_apo fitted pass-through works", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_apo(y, D, X,
                  d = 1,
                  learners = learners,
                  ensemble_type = "average",
                  sample_folds = 2,
                  silent = TRUE)

  # Pass-through with average ensemble reproduces exactly
  fit2 <- ddml_apo(y, D, X,
                   d = 1,
                   learners = learners,
                   ensemble_type = "average",
                   sample_folds = 2,
                   silent = TRUE,
                   fitted = fit$fitted,
                   splits = fit$splits)
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-6)

  expect_error(
    ddml_apo(y, D, X,
             d = 1,
             learners = learners,
             ensemble_type = "average",
             sample_folds = 2,
             silent = TRUE,
             fitted = fit$fitted),
    "must be supplied when 'fitted' is supplied"
  )
})

