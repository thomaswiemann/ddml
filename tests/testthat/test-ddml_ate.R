test_that("ddml_ate computes with a single model", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_ate_fit)), 1)
})#TEST_THAT

test_that("ddml_ate computes with stratify = FALSE", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             stratify = FALSE,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_ate_fit)), 1)
})#TEST_THAT

test_that("ddml_ate computes with a single model and dependence", {
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
  ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             cluster_variable = cluster_variable,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_ate_fit)), 1)
})#TEST_THAT

test_that("ddml_ate computes with an ensemble procedure", {
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
  ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             ensemble_type = "ols",
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_ate_fit)), 1)
})#TEST_THAT

test_that("ddml_ate computes w/ multiple ensembles + custom weights", {
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
  ddml_ate_fit <- ddml_ate(y, D, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             custom_ensemble_weights = diag(1, 2),
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_ate_fit)), 6)
})#TEST_THAT

test_that("ddml_ate computes with multiple ensemble procedures & shortstack", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols))
  # Compute DDML PLM estimator
  ddml_ate_fit <- ddml_ate(y, D, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             shortstack = TRUE,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_ate_fit)), 4)
})#TEST_THAT

test_that("summary.ddml_ate computes with a single model", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Compute inference results & test print
  inf_res <- summary(ddml_ate_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(1, 4, 1))
})#TEST_THAT

test_that("summary.ddml_ate computes with a single model and dependence", {
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
  ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             cluster_variable = cluster_variable,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Compute inference results & test print
  inf_res <- summary(ddml_ate_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(1, 4, 1))
})#TEST_THAT

test_that("summary.ddml_ate computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols))
  # Compute DDML PLM estimator
  ddml_ate_fit <- ddml_ate(y, D, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = TRUE)
  # Compute inference results & test print
  inf_res <- summary(ddml_ate_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(1, 4, 4))
})#TEST_THAT

test_that("ddml_ate fitted pass-through works", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_ate(y, D, X,
                  learners = learners,
                  ensemble_type = "average",
                  sample_folds = 2,
                  silent = TRUE)

  # Pass-through with average ensemble reproduces exactly
  fit2 <- ddml_ate(y, D, X,
                   learners = learners,
                   ensemble_type = "average",
                   sample_folds = 2,
                   silent = TRUE,
                   fitted = fit$fitted,
                   splits = fit$splits)
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-6)

  expect_error(
    ddml_ate(y, D, X,
             learners = learners,
             ensemble_type = "average",
             sample_folds = 2,
             silent = TRUE,
             fitted = fit$fitted),
    "must be supplied when 'fitted' is supplied"
  )
})

test_that("ddml_ate legacy grouped split args warn and still work", {
  set.seed(202)
  nobs <- 600
  X <- matrix(rnorm(nobs * 4), nobs, 4)
  D_tld <- X %*% c(0.7, 0.2, 0.1, 0) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% c(0.2, 0.3, 0.1, 0.4) + rnorm(nobs)
  learners <- list(what = ols)
  splits <- get_sample_splits(
    cluster_variable = seq_len(nobs),
    sample_folds = 3,
    cv_folds = 3,
    D = D,
    stratify = TRUE
  )
  fit_splits <- ddml_ate(
    y, D, X,
    learners = learners,
    sample_folds = 3,
    cv_folds = 3,
    splits = list(
      subsamples = splits$subsamples,
      subsamples_byD = splits$subsamples_byD,
      cv_subsamples = splits$cv_subsamples,
      cv_subsamples_byD = splits$cv_subsamples_byD
    ),
    silent = TRUE
  )
  expect_warning(
    fit_legacy <- ddml_ate(
      y, D, X,
      learners = learners,
      sample_folds = 3,
      cv_folds = 3,
      subsamples = splits$subsamples,
      subsamples_byD = splits$subsamples_byD,
      cv_subsamples = splits$cv_subsamples,
      cv_subsamples_byD = splits$cv_subsamples_byD,
      silent = TRUE
    ),
    "Deprecated split arguments detected"
  )
  expect_equal(coef(fit_legacy), coef(fit_splits), tolerance = 1e-8)
})

test_that("ddml_ate scores are mean-zero", {
  nobs <- 500
  set.seed(42)
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- 0.2 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  fit <- ddml_ate(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 3, silent = TRUE)
  score_mean <- mean(fit$scores[, , 1])
  expect_true(abs(score_mean) < 0.01,
              label = paste("score mean =", round(score_mean, 6)))
})

test_that("ddml_ate point estimate is close to true ATE", {
  nobs <- 2000
  set.seed(123)
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- 0.5 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- 1.0 * D + 0.3 * X[, 1] + rnorm(nobs)
  fit <- ddml_ate(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 5, silent = TRUE)
  expect_equal(coef(fit), c(ATE = 1.0), tolerance = 0.3)
})

test_that("ddml_ate computes with sparse matrices", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  learners <- list(list(what = ols), list(what = ols))
  ddml_ate_fit <- ddml_ate(y, D, as(X, "sparseMatrix"),
                           learners = learners,
                           ensemble_type = "ols",
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = TRUE)
  expect_equal(length(coef(ddml_ate_fit)), 1)
})#TEST_THAT

test_that("ddml_ate computes with parallel", {
  skip_on_cran()
  skip_if_not_installed("parallel")
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  learners <- list(what = ols)
  # Sequential
  res_seq <- ddml_ate(y, D, X,
                      learners = learners,
                      sample_folds = 3,
                      stratify = FALSE,
                      silent = TRUE)
  # Parallel with same splits
  res_par <- ddml_ate(y, D, X,
                      learners = learners,
                      sample_folds = 3,
                      stratify = FALSE,
                      splits = res_seq$splits,
                      silent = TRUE,
                      parallel = list(cores = 2))
  expect_equal(coef(res_par), coef(res_seq))
})#TEST_THAT
