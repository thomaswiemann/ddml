test_that("crossval_compute returns residuals (w/o instruments)", {
  # Simulate small dataset
  X <- matrix(rnorm(100*100), 100, 100) # Simulate features
  y <- 1 + X %*% (10*runif(100) * (runif(100) < 0.05)) + rnorm(100)
  # Define arguments
  test_sample <- sample(1:length(y), 33)
  learner <- list(fun = ols)
  # Compute cross-validation instance
  oos_resid <- crossval_compute(test_sample, learner,
                                y, X, Z = NULL)
  # Check output with expectations
  expect_equal(length(oos_resid), 33)
})#TEST_THAT

test_that("crossval returns residuals by learner (w/o instruments)", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*99), 100, 99)) # Simulate features
  nonzero_X <- (runif(100) < 0.05)
  y <- X %*% (10*runif(100) * nonzero_X) + rnorm(100)
  # Define arguments
  learners <- list(list(fun = ols),
                 list(fun = ols),
                 list(fun = ols,
                      assign_X = which(nonzero_X)))
  # Compute cross-validation instance
  cv_res <- crossval(y, X, Z = NULL,
                     learners,
                     cv_folds = 3,
                     silent = T)
  # Check output with expectations
  expect_equal(dim(cv_res$oos_resid), c(length(y), length(learners)))
})#TEST_THAT

test_that("crossval returns residuals by learner in correct order", {
  # Simulate small dataset
  n <- 100
  X <- rnorm(n)
  y <- X + rnorm(n)
  # split data to two folds and compute residuals manually
  subsample_list <- generate_subsamples(n, 2)
  oos_res_cv <- numeric(n)
  for(i in seq_along(subsample_list)) {
    idx_i <- subsample_list[[i]]
    ols_fit <- ols(y[-idx_i], X[-idx_i])
    oos_res_cv[idx_i] <- y[idx_i] - X[idx_i] %*% ols_fit$coef
  }#FOR
  # Compute cross-validation with crossval using the same subsamples
  cv_res <- crossval(y, as.matrix(X),
                     learners = list(list(fun = ols),
                                     list(fun = ols)),
                     cv_subsamples = subsample_list,
                     silent = T)
  # Check output with expectations
  expect_equal(cv_res$oos_resid[, 1], oos_res_cv)
  expect_equal(cv_res$oos_resid[, 2], oos_res_cv)
})#TEST_THAT

test_that("crossval returns residuals by learner (w/ instruments)", {
  # Simulate small dataset
  X <- cbind(1, matrix(rnorm(100*39), 100, 39))
  Z <- matrix(rnorm(100*10), 100, 10)
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + rnorm(100)
  # Define arguments
  learners <- list(list(fun = ols),
                 list(fun = ols),
                 list(fun = ols))
  # Compute cross-validation instance
  cv_res <- crossval(D, X, Z,
                     learners,
                     cv_folds = 3,
                     silent = T)
  # Check output with expectations
  expect_equal(all(cv_res$oos_resid[,2] == cv_res$oos_resid[,3]), TRUE)
})#TEST_THAT
