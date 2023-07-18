test_that("ddml_ate computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = mdl_glmnet,
                   args = list(alpha = 0.5))
  ddml_ate_fit <- ddml_ate(y, D, X,
                           learners = learners,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 1)
})#TEST_THAT

test_that("ddml_plm computes with an ensemble procedure", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_ate_fit <- ddml_ate(y, D, X,
                           learners = learners,
                           ensemble_type = c("ols"),
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 1)
})#TEST_THAT

test_that("ddml_plm computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_ate_fit <- ddml_ate(y, D, X,
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 5)
})#TEST_THAT

test_that("ddml_plm computes with multiple ensemble procedures & shortstack", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_ate_fit <- ddml_ate(y, D, X,
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = TRUE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 5)
})#TEST_THAT