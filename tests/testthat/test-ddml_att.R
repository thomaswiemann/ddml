test_that("ddml_att computes with a single model", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_att_fit <- ddml_att(y, D, X,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_att_fit$att), 1)
})#TEST_THAT

test_that("ddml_att computes with stratify = FALSE", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_att_fit <- ddml_att(y, D, X,
                             learners = learners,
                             stratify = FALSE,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_att_fit$att), 1)
})#TEST_THAT

test_that("ddml_att computes with a single model and dependence", {
  # Simulate small dataset
  n_cluster <- 200
  nobs <- 500
  X <- cbind(1, matrix(rnorm(n_cluster*39), n_cluster, 39))
  D_tld <-  X %*% runif(40) + rnorm(n_cluster)
  fun <- stepfun(quantile(D_tld, probs = 0.5), c(0, 1))
  D <- fun(D_tld)
  cluster_variable <- sample(1:n_cluster, nobs, replace = TRUE)
  D <- D[cluster_variable, drop = F]
  X <- X[cluster_variable, , drop = F]
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_att_fit <- ddml_att(y, D, X,
                             learners = learners,
                             cluster_variable = cluster_variable,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_att_fit$att), 1)
})#TEST_THAT

test_that("ddml_att computes with an ensemble procedure", {
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
  ddml_att_fit <- ddml_att(y, D, X,
                             learners = learners,
                             ensemble_type = "ols",
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_att_fit$att), 1)
})#TEST_THAT

test_that("ddml_att computes w/ multiple ensembles + custom weights", {
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
  ddml_att_fit <- ddml_att(y, D, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             custom_ensemble_weights = diag(1, 2),
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_att_fit$att), 6)
})#TEST_THAT

test_that("ddml_att computes w/ multp ensembles, custom weights + shortstack", {
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
  ddml_att_fit <- ddml_att(y, D, X,
                             learners,
                             ensemble_type = c("ols", "average"),
                             shortstack = TRUE,
                             cv_folds = 3,
                             custom_ensemble_weights = diag(1, 2),
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_att_fit$att), 4)
})#TEST_THAT

test_that("summary.ddml_att computes with a single model", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_att_fit <- ddml_att(y, D, X,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Compute inference results & test print
  inf_res <- summary(ddml_att_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$inf_results), c(1, 4, 1))
})#TEST_THAT

test_that("summary.ddml_att computes with a single model and dependence", {
  # Simulate small dataset
  n_cluster <- 200
  nobs <- 500
  X <- cbind(1, matrix(rnorm(n_cluster*39), n_cluster, 39))
  D_tld <-  X %*% runif(40) + rnorm(n_cluster)
  fun <- stepfun(quantile(D_tld, probs = 0.5), c(0, 1))
  D <- fun(D_tld)
  cluster_variable <- sample(1:n_cluster, nobs, replace = TRUE)
  D <- D[cluster_variable, drop = F]
  X <- X[cluster_variable, , drop = F]
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_att_fit <- ddml_att(y, D, X,
                             learners = learners,
                             cluster_variable = cluster_variable,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Compute inference results & test print
  inf_res <- summary(ddml_att_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$inf_results), c(1, 4, 1))
})#TEST_THAT

test_that("summary.ddml_att computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols))
  # Compute DDML PLM estimator
  ddml_att_fit <- ddml_att(y, D, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Compute inference results & test print
  inf_res <- summary(ddml_att_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$inf_results), c(1, 4, 4))
})#TEST_THAT

test_that("ddml_att fitted pass-through works", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_att(y, D, X,
                  learners = learners,
                  ensemble_type = "average",
                  sample_folds = 2,
                  silent = TRUE)

  fit2 <- ddml_att(y, D, X,
                   learners = learners,
                   ensemble_type = "average",
                   sample_folds = 2,
                   silent = TRUE,
                   fitted = fit$fitted,
                   splits = fit$splits)
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-6)

  expect_error(
    ddml_att(y, D, X,
             learners = learners,
             ensemble_type = "average",
             sample_folds = 2,
             silent = TRUE,
             fitted = fit$fitted),
    "must be supplied when 'fitted' is supplied"
  )
})#TEST_THAT
