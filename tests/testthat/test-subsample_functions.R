test_that("subsample construction works with multi-valued D", {
  # Simulate small dataset
  nobs <- 500
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  fun <- stepfun(quantile(D_tld, probs = c(0.25, 0.5, 0.75)), c(1, 2, 3, 4))
  D <- fun(D_tld)
  y <- D + X %*% runif(40) + rnorm(nobs)

  # Compute crossfit indices w/o splitting by D
  crossfit_indices_1 <- get_crossfit_indices(cluster_variable = 1:500,
                                           sample_folds = 5, cv_folds = 3,
                                           D = NULL)

  # Compute crossfit indices w/ splitting by D
  crossfit_indices_2 <- get_crossfit_indices(cluster_variable = 1:500,
                                           sample_folds = 5, cv_folds = 3,
                                           D = D)

  # Check that subsamples and cv subsamples are of matching size
  expect_equal(nobs - length(crossfit_indices_1$subsamples[[1]]),
               length(unlist(crossfit_indices_1$cv_subsamples_list[[1]])))
  expect_equal(nobs - length(crossfit_indices_2$subsamples[[1]]),
               length(unlist(crossfit_indices_2$cv_subsamples_list[[1]])))
  expect_equal(sum(D==1) - length(crossfit_indices_2$subsamples_byD[[1]][[1]]),
               length(unlist(crossfit_indices_2$cv_subsamples_byD[[1]][[1]])))
})#TEST_THAT

test_that("subsample construction works with multi-valued D and dependence", {
  # Simulate small dataset
  n_cluster <- 200
  nobs <- 500
  X <- cbind(1, matrix(rnorm(n_cluster*39), n_cluster, 39))
  D_tld <-  X %*% runif(40) + rnorm(n_cluster)
  fun <- stepfun(quantile(D_tld, probs = c(0.25, 0.5, 0.75)), c(1, 2, 3, 4))
  D <- fun(D_tld)
  cluster_variable <- sample(1:n_cluster, nobs, replace = TRUE)
  D <- D[cluster_variable]
  X <- X[cluster_variable, ]
  y <- D + X %*% runif(40) + rnorm(nobs)

  # Compute crossfit indices w/o splitting by D
  crossfit_indices_1 <- get_crossfit_indices(cluster_variable =
                                               cluster_variable,
                                             sample_folds = 5, cv_folds = 3,
                                             D = NULL)
  expect_identical(
    setdiff(unique(cluster_variable[crossfit_indices_1$subsamples[[1]]]),
          unique(cluster_variable[crossfit_indices_1$subsamples[[2]]])),
    unique(cluster_variable[crossfit_indices_1$subsamples[[1]]]))

  # Compute crossfit indices w/ splitting by D
  crossfit_indices_2 <- get_crossfit_indices(cluster_variable =
                                               cluster_variable,
                                             sample_folds = 5, cv_folds = 3,
                                             D = D)

  # Check that cluster variables are unique across folds
  expect_identical(
  setdiff(cluster_variable[D==1][crossfit_indices_2$subsamples_byD[[1]][[1]]],
          cluster_variable[D==1][crossfit_indices_2$subsamples_byD[[1]][[2]]]),
  unique(cluster_variable[D==1][crossfit_indices_2$subsamples_byD[[1]][[1]]]))
  expect_identical(
    setdiff(cluster_variable[D==2][crossfit_indices_2$subsamples_byD[[2]][[1]]],
            cluster_variable[D==2][crossfit_indices_2$subsamples_byD[[2]][[2]]]),
    unique(cluster_variable[D==2][crossfit_indices_2$subsamples_byD[[2]][[1]]]))
  expect_identical(
    setdiff(unique(cluster_variable[crossfit_indices_2$subsamples[[1]]]),
            unique(cluster_variable[crossfit_indices_2$subsamples[[2]]])),
    unique(cluster_variable[crossfit_indices_2$subsamples[[1]]]))

  # Check that subsamples and cv subsamples are of matching size
  expect_equal(nobs - length(crossfit_indices_1$subsamples[[1]]),
               length(unlist(crossfit_indices_1$cv_subsamples_list[[1]])))
  expect_equal(nobs - length(crossfit_indices_2$subsamples[[1]]),
               length(unlist(crossfit_indices_2$cv_subsamples_list[[1]])))
  expect_equal(sum(D==1) - length(crossfit_indices_2$subsamples_byD[[1]][[1]]),
               length(unlist(crossfit_indices_2$cv_subsamples_byD[[1]][[1]])))
})#TEST_THAT

test_that("auxiliary_X construction works with multi-valued D", {
  # Simulate small dataset
  nobs <- 500
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  fun <- stepfun(quantile(D_tld, probs = c(0.25, 0.5, 0.75)), c(1, 2, 3, 4))
  D <- fun(D_tld)
  y <- D + X %*% runif(40) + rnorm(nobs)

  # Compute crossfit indices w/ splitting by D
  crossfit_indices <- get_crossfit_indices(cluster_variable = 1:500,
                                             sample_folds = 5, cv_folds = 3,
                                             D = D)

  # Get auxiliary_indx
  auxiliary_indx <- get_auxiliary_indx(crossfit_indices$subsamples_byD, D)

  # Get auxiliary_X for first sample_folds
  auxiliary_X_d <- get_auxiliary_X(auxiliary_indx[[1]], X)

  # Check that auxiliary X is of correct size
  expect_equal(length(unlist(auxiliary_indx[[1]])),
               nobs - length(crossfit_indices$subsamples[[1]]) -
                 length(crossfit_indices$subsamples_byD[[1]][[1]]))
  expect_equal(dim(auxiliary_X_d[[1]])[1], length(auxiliary_indx[[1]][[1]]))
})#TEST_THAT


