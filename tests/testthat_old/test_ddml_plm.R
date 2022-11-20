library(ddml)
context("Testing ddml_plm.")

sim_dat <- function(nobs) {
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Organize and return output
  output <- list(y = y, D = D, X = X)
  return(output)
}#SIM_DAT

test_that("ddml_plm computes with a single model", {
  # Simulate small dataset
  dat <- sim_dat(100)
  # Define arguments
  learners <- list(what = mdl_glmnet,
                 args = list(alpha = 0.5))
  ddml_plm_fit <- ddml_plm(dat$y, dat$D, dat$X,
                         learners = learners,
                         cv_folds = 3,
                         sample_folds = 3,
                         silent = T)
  # Check output with expectations
  expect_equal(length(ddml_plm_fit$coef), 1)
})#TEST_THAT

test_that("ddml_plm computes with an ensemble procedure", {
  # Simulate small dataset
  dat <- sim_dat(100)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(dat$y, dat$D, dat$X,
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
  dat <- sim_dat(100)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                      args = list(alpha = 0.5)),
                 list(fun = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(dat$y, dat$D, dat$X,
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
