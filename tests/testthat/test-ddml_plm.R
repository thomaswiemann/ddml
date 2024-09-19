test_that("ddml_plm computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = mdl_glmnet,
                   args = list(alpha = 0.5))
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 1)
})#TEST_THAT

test_that("ddml_plm computes with clustered observations", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  cluster_variable <- sample(1:100, nobs, replace = T)
  # Define arguments
  learners <- list(what = mdl_glmnet,
                   args = list(alpha = 0.5))
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           cv_folds = 3,
                           sample_folds = 3,
                           cluster_variable = cluster_variable,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 1)
})#TEST_THAT

test_that("ddml_plm computes with an ensemble procedure", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           ensemble_type = "ols",
                           shortstack = FALSE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 1)
})#TEST_THAT

test_that("ddml_plm computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = FALSE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 5)
})#TEST_THAT

test_that("ddml_plm computes with multiple ensemble procedures & sparse mats", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, as(X, "sparseMatrix"),
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = FALSE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 5)
})#TEST_THAT
test_that("ddml_plm computes w/ an ensemble procedure & shortstacking", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           ensemble_type = "ols",
                           shortstack = TRUE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 1)
})#TEST_THAT

test_that("ddml_plm computes w/ multiple ensemble procedures & shortstacking", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = TRUE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 5)
})#TEST_THAT

test_that("ddml_plm computes w/ ensemble procedures & custom weights", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = TRUE,
                           cv_folds = 3,
                           custom_ensemble_weights = diag(1, length(learners)),
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 7)
})#TEST_THAT

test_that("summary.ddml_plm computes with a single model", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(what = mdl_glmnet,
                   args = list(alpha = 0.5))
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  inf_res <- summary(ddml_plm_fit, type = "HC1")
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_equal(length(inf_res), 8)
})#TEST_THAT

test_that("summary.ddml_plm computes with a single model and dependence", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  cluster_variable <- sample(1:100, nobs, replace = T)
  # Define arguments
  learners <- list(what = mdl_glmnet,
                   args = list(alpha = 0.5))
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           cv_folds = 3,
                           sample_folds = 3,
                           cluster_variable = cluster_variable,
                           silent = T)
  inf_res <- summary(ddml_plm_fit)
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_equal(length(inf_res), 8)
})#TEST_THAT

test_that("summary.ddml_plm computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = FALSE,
                           cv_folds = 3,
                           custom_ensemble_weights = diag(1, length(learners)),
                           sample_folds = 3,
                           silent = T)
  inf_res <- summary(ddml_plm_fit, type = "HC1")
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_equal(length(inf_res), 8 * 7)
})#TEST_THAT

test_that("ddml_plm computes with an ensemble procedure and multivariate D", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  cbind(X %*% runif(40) + rnorm(nobs), rnorm(nobs))
  y <- rowSums(D) + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           ensemble_type = "ols",
                           shortstack = FALSE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 2)
})#TEST_THAT

test_that("ddml_plm computes with multiple ensemble types and multivariate D", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  cbind(X %*% runif(40) + rnorm(nobs), rnorm(nobs))
  y <- rowSums(D) + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = FALSE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 10)
})#TEST_THAT
