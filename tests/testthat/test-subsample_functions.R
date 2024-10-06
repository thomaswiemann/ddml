test_that("get_cf_indices_stratified computes ", {

  set.seed(20241004)

  # Simulate imbalanced panel with varying treatment status
  nclusters <- 100
  obs_per_cluster <- sample(1:5, nclusters, replace = TRUE)
  nobs <- sum(obs_per_cluster)
  cluster_variable <- rep(1:nclusters, times = obs_per_cluster)
  D <- sample(x = c(1, 2, 3), size = nobs, replace = TRUE,
              prob = c(0.1, 0.2, 0.7))

  # Compute indices balanced by observations and by clusters
  cf_obs <- get_cf_indices_stratified(cluster_variable = cluster_variable,
                                      sample_folds = 10,  D = D,
                                      balance_on = "observations")
  cf_cl <-  get_cf_indices_stratified(cluster_variable = cluster_variable,
                                      sample_folds = 10, D = D,
                                      balance_on = "clusters")

  # Check whether subsample dimensions are correct
  expect_equal(length(unlist(cf_obs$subsamples)), length(D))
  expect_equal(length(unlist(cf_obs$subsamples_byD[[1]])), sum(D==1))
  expect_equal(length(unlist(cf_obs$subsamples_byD[[2]])), sum(D==2))
  expect_equal(length(unlist(cf_obs$subsamples_byD[[3]])), sum(D==3))
  expect_equal(length(unlist(cf_cl$subsamples)), length(D))
  expect_equal(length(unlist(cf_cl$subsamples_byD[[1]])), sum(D==1))
  expect_equal(length(unlist(cf_cl$subsamples_byD[[2]])), sum(D==2))
  expect_equal(length(unlist(cf_cl$subsamples_byD[[3]])), sum(D==3))

  # Compute indices w/o panel
  cf_alt <- get_cf_indices_stratified(cluster_variable = seq_along(D),
                                      sample_folds = 10,  D = D)

  # Check whether subsample dimensions are correct
  expect_equal(length(unlist(cf_alt$subsamples)), length(D))
  expect_equal(length(unlist(cf_alt$subsamples_byD[[1]])), sum(D==1))
  expect_equal(length(unlist(cf_alt$subsamples_byD[[2]])), sum(D==2))
  expect_equal(length(unlist(cf_alt$subsamples_byD[[3]])), sum(D==3))
})#TEST_THAT

test_that("get_cf_indices_simple computes ", {

  set.seed(20241006)

  # Simulate imbalanced panel with varying treatment status
  nclusters <- 100
  obs_per_cluster <- sample(1:5, nclusters, replace = TRUE)
  nobs <- sum(obs_per_cluster)
  cluster_variable <- rep(1:nclusters, times = obs_per_cluster)
  D <- sample(x = c(1, 2, 3), size = nobs, replace = TRUE,
              prob = c(0.1, 0.2, 0.7))

  # Compute indices balanced by observations and by clusters
  cf_cl <- get_cf_indices_simple(cluster_variable = cluster_variable,
                                 sample_folds = 10,  D = D, by_D = TRUE)
  cf_cross <- get_cf_indices_simple(cluster_variable = seq_along(D),
                                    sample_folds = 10,  D = D, by_D = TRUE)

  # Check whether subsample dimensions are correct
  expect_equal(length(unlist(cf_cl$subsamples)), length(D))
  expect_equal(length(unlist(cf_cl$subsamples_byD[[1]])), sum(D==1))
  expect_equal(length(unlist(cf_cl$subsamples_byD[[2]])), sum(D==2))
  expect_equal(length(unlist(cf_cl$subsamples_byD[[3]])), sum(D==3))
  expect_equal(length(unlist(cf_cross$subsamples)), length(D))
  expect_equal(length(unlist(cf_cross$subsamples_byD[[1]])), sum(D==1))
  expect_equal(length(unlist(cf_cross$subsamples_byD[[2]])), sum(D==2))
  expect_equal(length(unlist(cf_cross$subsamples_byD[[3]])), sum(D==3))
})#TEST_THAT

test_that("get_crossfit_indices computes", {

  set.seed(20241006)

  # Simulate imbalanced panel with varying treatment status
  nclusters <- 50
  obs_per_cluster <- sample(1:5, nclusters, replace = TRUE)
  nobs <- sum(obs_per_cluster)
  cluster_variable <- rep(1:nclusters, times = obs_per_cluster)
  D <- sample(x = c(1, 2, 3), size = nobs, replace = TRUE,
              prob = c(0.1, 0.2, 0.7))

  # Compute indices, check for warnings
  cf_0 <- get_crossfit_indices(cluster_variable, sample_folds = 2,
                               D = NULL, by_D = FALSE, stratify = FALSE)
  cf_1 <- get_crossfit_indices(cluster_variable, sample_folds = 6,
                               D = D, by_D = TRUE, stratify = FALSE)
  cf_2 <- get_crossfit_indices(cluster_variable, sample_folds = 7,
                               D = D, by_D = TRUE, stratify = TRUE)

  # Check that warnings are also given when subsamples are user-provided
  cf_00 <- get_crossfit_indices(cluster_variable, sample_folds = 2,
                                D = NULL, by_D = FALSE, stratify = FALSE,
                                subsamples = cf_0$subsamples)
  expect_identical(cf_0$subsamples, cf_00$subsamples)
  cf_11 <- get_crossfit_indices(cluster_variable, sample_folds = 6,
                                D = D, by_D = TRUE, stratify = FALSE,
                                subsamples = cf_1$subsamples,
                                subsamples_byD = cf_1$subsamples_byD)
  expect_identical(cf_1$subsamples, cf_11$subsamples)
  expect_identical(cf_1$subsamples_byD, cf_11$subsamples_byD)
})#TEST_THAT


test_that("get_crossval_indices computes", {

  set.seed(20241006)

  # Simulate imbalanced panel with varying treatment status
  nclusters <- 100
  obs_per_cluster <- sample(1:5, nclusters, replace = TRUE)
  nobs <- sum(obs_per_cluster)
  cluster_variable <- rep(1:nclusters, times = obs_per_cluster)
  D <- sample(x = c(1, 2, 3), size = nobs, replace = TRUE,
              prob = c(0.1, 0.2, 0.7))

  # Compute indices, check for warnings

  cf_indx <- get_crossfit_indices(cluster_variable, sample_folds = 5,
                                  D = D, by_D = TRUE, stratify = TRUE)
  cv_indx <- get_crossval_indices(cf_indx$subsamples,
                                  cluster_variable,
                                  cv_folds = 3,
                                  D = D,
                                  by_D = TRUE,
                                  stratify = TRUE,
                                  cv_subsamples = NULL,
                                  cv_subsamples_byD = NULL)

  cv_indx_2 <- get_crossval_indices(cf_indx$subsamples,
                                    cluster_variable,
                                    cv_folds = 3,
                                    D = D,
                                    by_D = TRUE,
                                    stratify = TRUE,
                                    cv_subsamples = cv_indx$cv_subsamples,
                                    cv_subsamples_byD =
                                      cv_indx$cv_subsamples_byD)

  expect_equal(length(cv_indx$cv_subsamples_byD), length(unique(D)))
  expect_equal(length(cv_indx_2$cv_subsamples_byD), length(unique(D)))





  # Check that is also works w/o D
  cf_indx2 <- get_crossfit_indices(cluster_variable, sample_folds = 3,
                                   by_D = FALSE)
  cv_indx2 <- get_crossval_indices(cf_indx2$subsamples,
                                   cluster_variable,
                                   cv_folds = 5,
                                   cv_subsamples = NULL,
                                   cv_subsamples_byD = NULL)
  expect_equal(class(cv_indx2$cv_subsamples_byD[[1]][[1]]), "NULL")
})#TEST_THAT

test_that("get_all_indices computes", {

  set.seed(20241006)

  # Simulate imbalanced panel with varying treatment status
  nclusters <- 100
  obs_per_cluster <- sample(1:5, nclusters, replace = TRUE)
  nobs <- sum(obs_per_cluster)
  cluster_variable <- rep(1:nclusters, times = obs_per_cluster)
  D <- sample(x = c(1, 2, 3), size = nobs, replace = TRUE,
              prob = c(0.1, 0.2, 0.7))

  # Compute indices, check for warnings

  all_indx <- get_all_indx(cluster_variable=cluster_variable,
                           sample_folds = 10,
                           cv_folds = 4,
                           D = D,
                           by_D = TRUE,
                           stratify = TRUE,
                           balance_on = "clusters",
                           subsamples = NULL,
                           subsamples_byD = NULL,
                           cv_subsamples = NULL,
                           cv_subsamples_byD = NULL,
                           compute_cv_indices = TRUE,
                           compute_aux_X_indices = TRUE)

  all_indx2 <- get_all_indx(cluster_variable=cluster_variable,
                            sample_folds = 10,
                            cv_folds = 4,
                            D = D,
                            by_D = TRUE,
                            stratify = TRUE,
                            balance_on = "clusters",
                            subsamples = all_indx$subsamples,
                            subsamples_byD = all_indx$subsamples_byD,
                            cv_subsamples = all_indx$cv_subsamples,
                            cv_subsamples_byD = all_indx$cv_subsamples_byD,
                            compute_cv_indices = TRUE,
                            compute_aux_X_indices = TRUE)

  # Check whether pass-through worked
  expect_identical(all_indx$subsamples, all_indx2$subsamples)
  expect_identical(all_indx$subsamples_byD, all_indx2$subsamples_byD)
  expect_identical(all_indx$cv_subsamples, all_indx2$cv_subsamples)
  expect_identical(all_indx$cv_subsamples_byD, all_indx2$cv_subsamples_byD)

  # Check that is also works w/o D
  all_indx <- get_all_indx(cluster_variable=cluster_variable,
                           sample_folds = 10,
                           cv_folds = 10,
                           compute_cv_indices = TRUE)
  expect_equal(length(unlist(all_indx$cv_subsamples[[1]])),
               nobs - length(all_indx$subsamples[[1]]))
})#TEST_THAT


test_that("auxiliary_X construction works with multi-valued D", {
  # Simulate small dataset
  nobs <- 294
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  fun <- stepfun(quantile(D_tld, probs = c(0.25, 0.5, 0.75)), c(1, 2, 3, 4))
  D <- fun(D_tld)

  # Compute indices, check for warnings

  all_indx <- get_all_indx(cluster_variable=seq_along(D),
                           sample_folds = 10,
                           cv_folds = 10,
                           D = D,
                           by_D = TRUE,
                           stratify = TRUE,
                           balance_on = "clusters",
                           subsamples = NULL,
                           subsamples_byD = NULL,
                           cv_subsamples = NULL,
                           cv_subsamples_byD = NULL,
                           compute_cv_indices = FALSE,
                           compute_aux_X_indices = TRUE)

  # Get auxiliary_indx
  auxiliary_indx <- get_auxiliary_indx(all_indx$subsamples_byD, D)

  # Get auxiliary_X for first sample fold
  auxiliary_X_d <- get_auxiliary_X(auxiliary_indx[[1]], X)

  # Check that auxiliary X is of correct size
  expect_equal(length(auxiliary_indx[[1]][[1]]),
               sum(D[all_indx$subsamples[[1]]]!=1))
  expect_equal(length(auxiliary_indx[[2]][[1]]),
               sum(D[all_indx$subsamples[[1]]]!=2))
  expect_equal(length(auxiliary_indx[[3]][[1]]),
               sum(D[all_indx$subsamples[[1]]]!=3))
  expect_equal(dim(auxiliary_X_d[[1]])[1], length(auxiliary_indx[[1]][[1]]))
  expect_identical(all_indx$aux_indx, auxiliary_indx)
})#TEST_THAT
