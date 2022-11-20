library(ddml)
context("Testing ddml_epliv.")

sim_dat <- function(nobs) {
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10) # overidentified
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Organize and return output
  output <- list(y = y, D = D, Z = Z, X = X)
  return(output)
}#SIM_DAT

test_that("ddml_epliv computes with a single model", {
  # Simulate small dataset
  dat <- sim_dat(1000)
  # Define arguments
  learners <- list(what = mdl_glmnet,
                   args = list(alpha = 0.5))
  ddml_epliv_fit <- ddml_epliv(dat$y, dat$D, dat$Z, dat$X,
                               learners,
                               sample_folds = 2,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_epliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_epliv computes with an ensemble procedure", {
  # Simulate small dataset
  dat <- sim_dat(500)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_epliv_fit <- ddml_epliv(dat$y, dat$D, dat$Z, dat$X,
                               learners,
                               ensemble_type = c("stacking"),
                               cv_folds = 2,
                               sample_folds = 2,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_epliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_epliv computes with stacking w/o enforcing the LIE", {
  # Simulate small dataset
  dat <- sim_dat(500)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_epliv_fit <- ddml_epliv(dat$y, dat$D, dat$Z, dat$X,
                               learners,
                               ensemble_type = c("stacking"),
                               sample_folds = 2,
                               cv_folds = 2,
                               enforce_LIE = FALSE,
                               silent = F)
  # Check output with expectations
  expect_equal(length(ddml_epliv_fit$coef), 1)
})#TEST_THAT

test_that("ddml_epliv computes with multiple ensemble procedures", {
  # Simulate small dataset
  dat <- sim_dat(500)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_epliv_fit <- ddml_epliv(dat$y, dat$D, dat$Z, dat$X,
                               learners,
                               ensemble_type = c("stacking", "stacking_nn",
                                                 "stacking_01",
                                                 "stacking_best", "average"),
                               cv_folds = 2,
                               sample_folds = 2,
                               silent = T)

  # Check output with expectations
  expect_equal(length(ddml_epliv_fit$coef), 5)
})#TEST_THAT

test_that("ddml_epliv computes with multiple ensembles w/o the LIE", {
  # Simulate small dataset
  dat <- sim_dat(500)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_epliv_fit <- ddml_epliv(dat$y, dat$D, dat$Z, dat$X,
                               learners,
                               ensemble_type = c("stacking", "stacking_nn",
                                                 "stacking_01",
                                                 "stacking_best", "average"),
                               cv_folds = 2,
                               sample_folds = 2,
                               enforce_LIE = FALSE,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_epliv_fit$coef), 5)
})#TEST_THAT

test_that("ddml_epliv computes with multiple ensembles and sparse matrices", {
  # Simulate small dataset
  dat <- sim_dat(500)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  learners <- list(list(fun = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_epliv_fit <- ddml_epliv(dat$y, dat$D,
                               as(dat$Z, "sparseMatrix"),
                               as(dat$X, "sparseMatrix"),
                               learners,
                               ensemble_type = c("stacking", "stacking_nn",
                                                 "stacking_01",
                                                 "stacking_best", "average"),
                               cv_folds = 2,
                               sample_folds = 2,
                               silent = T)
  # Check output with expectations
  expect_equal(length(ddml_epliv_fit$coef), 5)
})#TEST_THAT

test_that("ddml_epliv computes with different sets of learners", {
  # Simulate small dataset
  dat <- sim_dat(500)
  ncol_X <- ncol(dat$X)
  ncol_Z <- ncol(dat$Z)
  # Define arguments
  learners <- list(list(fun = ols),
                   list(fun = ols),
                   list(fun = ols))
  learners_DXZ <- list(list(fun = mdl_glmnet,
                            args = list(alpha = 0.5)),
                       list(fun = mdl_glmnet,
                            args = list(alpha = 1)))
  learners_DX <- list(list(fun = ols),
                      list(fun = mdl_glmnet,
                           args = list(alpha = 0.5)),
                      list(fun = ols))
  # Compute LIE-conform DDML IV estimator
  ddml_epliv_fit <- ddml_epliv(dat$y, dat$D, dat$Z, dat$X,
                               learners,
                               learners_DXZ = learners_DXZ,
                               learners_DX = learners_DX,
                               ensemble_type = c("stacking", "stacking_nn",
                                                 "stacking_01",
                                                 "stacking_best", "average"),
                               cv_folds = 2,
                               sample_folds = 2,
                               silent = T)

  # Check output with expectations
  expect_equal(length(ddml_epliv_fit$coef), 5)
})#TEST_THAT
