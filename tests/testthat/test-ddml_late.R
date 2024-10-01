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
  expect_warning({
    ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)
  })
  # Check output with expectations
  expect_equal(length(ddml_late_fit$late), 1)
})#TEST_THAT

test_that("ddml_late computes with a single model and dependence", {
  # Simulate small dataset
  n_cluster <- 250
  nobs <- 500
  X <- cbind(1, matrix(rnorm(n_cluster*39), n_cluster, 39))
  Z_tld <-  X %*% runif(40) + rnorm(n_cluster)
  fun <- stepfun(quantile(Z_tld, probs = 0.5), c(0, 1))
  Z <- fun(Z_tld)
  cluster_variable <- sample(1:n_cluster, nobs, replace = TRUE)
  Z <- Z[cluster_variable, drop = F]
  X <- X[cluster_variable, , drop = F]
  eps <- rnorm(nobs)
  D <- Z + X %*% runif(40) + eps
  y <- D + X %*% runif(40) + 0.1 * eps + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  expect_warning({
    ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               cluster_variable = cluster_variable,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)
  })
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
  expect_warning({
    ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)
  })
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
  expect_warning({
    ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               ensemble_type = "ols",
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)
  })
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
  expect_warning({
    ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "nnls1",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               custom_ensemble_weights = diag(1, 2),
                               sample_folds = 3,
                               silent = T)
  })
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
  expect_warning({
    ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners,
                               ensemble_type = c("ols", "nnls",
                                                 "nnls1",
                                                 "singlebest", "average"),
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)
  })
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
  expect_warning({
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
  })
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
  expect_warning({
    ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)
  })
  # Compute inference results & test print
  inf_res <- summary(ddml_late_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_equal(length(inf_res), 4)
})#TEST_THAT

test_that("summary.ddml_late computes with a single model and dependence", {
  # Simulate small dataset
  n_cluster <- 250
  nobs <- 500
  X <- cbind(1, matrix(rnorm(n_cluster*39), n_cluster, 39))
  Z_tld <-  X %*% runif(40) + rnorm(n_cluster)
  fun <- stepfun(quantile(Z_tld, probs = 0.5), c(0, 1))
  Z <- fun(Z_tld)
  cluster_variable <- sample(1:n_cluster, nobs, replace = TRUE)
  Z <- Z[cluster_variable, drop = F]
  X <- X[cluster_variable, , drop = F]
  eps <- rnorm(nobs)
  D <- Z + X %*% runif(40) + eps
  y <- D + X %*% runif(40) + 0.1 * eps + rnorm(nobs)
  # Define arguments
  learners <- list(what = ols)
  expect_warning({
    ddml_late_fit <- ddml_late(y, D, Z, X,
                               learners = learners,
                               cluster_variable = cluster_variable,
                               cv_folds = 3,
                               sample_folds = 3,
                               silent = T)
  })
  # Compute inference results & test print
  inf_res <- summary(ddml_late_fit)
  capture_output({print(inf_res)}, print = FALSE)
  # Check output with expectations
  expect_equal(length(inf_res), 4)
})#TEST_THAT
