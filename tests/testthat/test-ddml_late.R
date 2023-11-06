test_that("ddml_late computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z_tld <-  X %*% runif(40) + rnorm(nobs)
  Z <- 1 * (Z_tld > mean(Z_tld))
  D_tld <-  0.5 * (1 - 2 * Z) + 0.2 * X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_late_fit <- ddml_late(y, D, Z, X,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_late_fit$late), 1)
})#TEST_THAT

test_that("ddml_late computes with a single model & perfect non-compliance", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z_tld <-  X %*% runif(40) + rnorm(nobs)
  Z <- 1 * (Z_tld > mean(Z_tld))
  D_tld <-  0.5 * (1 - 2 * Z) + 0.2 * X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  D[Z == 0] <- 0
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_late_fit <- ddml_late(y, D, Z, X,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_late_fit$late), 1)
})#TEST_THAT

test_that("ddml_late computes with an ensemble procedure", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z_tld <-  X %*% runif(40) + rnorm(nobs)
  Z <- 1 * (Z_tld > mean(Z_tld))
  D_tld <-  0.25 * (1 - 2 * Z) + 0.2 * X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_late_fit <- ddml_late(y, D, Z, X,
                             learners = learners,
                             ensemble_type = "ols",
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_late_fit$late), 1)
})#TEST_THAT

test_that("ddml_late computes w/ multiple ensembles & custom weights", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z_tld <-  X %*% runif(40) + rnorm(nobs)
  Z <- 1 * (Z_tld > mean(Z_tld))
  D_tld <-  0.25 * (1 - 2 * Z) + 0.2 * X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_late_fit <- ddml_late(y, D, Z, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "nnls1",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             custom_ensemble_weights = diag(1, 2),
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_late_fit$late), 7)
})#TEST_THAT

test_that("ddml_late computes with multiple ensemble procedures + perfect compliance", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z_tld <-  X %*% runif(40) + rnorm(nobs)
  Z <- 1 * (Z_tld > mean(Z_tld))
  D_tld <-  0.25 * (1 - 2 * Z) + 0.2 * X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  D[Z == 1] <- 1 # perfect compliance
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_late_fit <- ddml_late(y, D, Z, X,
                             learners,
                             ensemble_type = c("ols", "nnls",
                                               "nnls1",
                                               "singlebest", "average"),
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_late_fit$late), 5)
})#TEST_THAT

test_that("ddml_late computes w/ mult ensembles, custom weights, & shortstack", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z_tld <-  X %*% runif(40) + rnorm(nobs)
  Z <- 1 * (Z_tld > mean(Z_tld))
  D_tld <-  0.25 * (1 - 2 * Z) + 0.2 * X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
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
                             silent = T)
  # Check output with expectations
  expect_equal(length(ddml_late_fit$late), 7)
})#TEST_THAT

test_that("summary.ddml_late computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z_tld <-  X %*% runif(40) + rnorm(nobs)
  Z <- 1 * (Z_tld > mean(Z_tld))
  D_tld <-  0.5 * (1 - 2 * Z) + 0.2 * X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  ddml_late_fit <- ddml_late(y, D, Z, X,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  capture_output({inf_res <- summary(ddml_late_fit)})
  # Check output with expectations
  expect_equal(length(inf_res), 4)
})#TEST_THAT
