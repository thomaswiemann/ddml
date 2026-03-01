test_that("ddml_ate computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  suppressWarnings({
    ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  })
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 1)
})#TEST_THAT

test_that("ddml_ate computes with stratify = FALSE", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  suppressWarnings({
    ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             stratify = FALSE,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  })
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 1)
})#TEST_THAT

test_that("ddml_ate computes with a single model and dependence", {
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
  suppressWarnings({
    ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             cluster_variable = cluster_variable,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  })
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 1)
})#TEST_THAT

test_that("ddml_ate computes with an ensemble procedure", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML PLM estimator
  suppressWarnings({
    ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             ensemble_type = "ols",
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  })
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 1)
})#TEST_THAT

test_that("ddml_ate computes w/ multiple ensembles + custom weights", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols))
  # Compute DDML PLM estimator
  suppressWarnings({
    ddml_ate_fit <- ddml_ate(y, D, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             custom_ensemble_weights = diag(1, 2),
                             sample_folds = 3,
                             silent = T)
  })
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 6)
})#TEST_THAT

test_that("ddml_ate computes with multiple ensemble procedures & shortstack", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols))
  # Compute DDML PLM estimator
  suppressWarnings({
    ddml_ate_fit <- ddml_ate(y, D, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             shortstack = TRUE,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  })
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 4)
})#TEST_THAT

test_that("summary.ddml_ate computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  suppressWarnings({
    ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  })
  # Compute inference results & test print
  inf_res <- summary(ddml_ate_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$inf_results), c(1, 4, 1))
})#TEST_THAT

test_that("summary.ddml_ate computes with a single model and dependence", {
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
  suppressWarnings({
    ddml_ate_fit <- ddml_ate(y, D, X,
                             learners = learners,
                             cluster_variable = cluster_variable,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  })
  # Compute inference results & test print
  inf_res <- summary(ddml_ate_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$inf_results), c(1, 4, 1))
})#TEST_THAT

test_that("summary.ddml_ate computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols))
  # Compute DDML PLM estimator
  suppressWarnings({
    ddml_ate_fit <- ddml_ate(y, D, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  })
  # Compute inference results & test print
  inf_res <- summary(ddml_ate_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$inf_results), c(1, 4, 4))
})#TEST_THAT

test_that("ddml_ate fitted pass-through works", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- 0.3 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.3 * X[, 1] + rnorm(nobs)

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
