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
                         ensemble_type = c("stacking"),
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
                           ensemble_type = c("stacking", "stacking_nn",
                                             "stacking_01",
                                        "stacking_best", "average"),
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 5)
})#TEST_THAT
