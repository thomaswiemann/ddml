# --- get_sample_splits: basic (no D) ---

test_that("get_sample_splits returns correct structure without D", {
  nobs <- 200
  cl <- seq_len(nobs)

  # Without CV
  res <- get_sample_splits(cl, sample_folds = 5)
  expect_equal(length(res$subsamples), 5)
  expect_equal(length(unlist(res$subsamples)), nobs)
  expect_null(res$cv_subsamples)
  expect_null(res$subsamples_byD)
  expect_null(res$cv_subsamples_byD)
  expect_null(res$aux_indx)

  # With CV
  res_cv <- get_sample_splits(cl, sample_folds = 5, cv_folds = 3)
  expect_equal(length(res_cv$subsamples), 5)
  expect_equal(length(res_cv$cv_subsamples), 5)
  expect_equal(length(res_cv$cv_subsamples[[1]]), 3)
  # CV fold sizes should sum to training set size
  expect_equal(
    length(unlist(res_cv$cv_subsamples[[1]])),
    nobs - length(res_cv$subsamples[[1]]))
})#TEST_THAT

# --- get_sample_splits: with D ---

test_that("get_sample_splits returns correct structure with D", {
  nobs <- 200
  cl <- seq_len(nobs)
  D <- rep(c(0, 1), each = nobs / 2)

  res <- get_sample_splits(cl, sample_folds = 5, cv_folds = 3, D = D)

  # Basic structure
  expect_equal(length(res$subsamples), 5)
  expect_equal(length(res$subsamples_byD), 2)
  expect_equal(length(res$subsamples_byD[[1]]), 5)
  expect_equal(length(res$cv_subsamples), 5)
  expect_equal(length(res$cv_subsamples_byD), 2)
  expect_false(is.null(res$aux_indx))
  expect_equal(length(res$aux_indx), 2)

  # CV fold sizes match training set for D=0 subsample
  n_D0 <- sum(D == 0)
  n_D0_fold1 <- length(res$subsamples_byD[[1]][[1]])
  expect_equal(
    length(unlist(res$cv_subsamples_byD[[1]][[1]])),
    n_D0 - n_D0_fold1)
})#TEST_THAT

# --- get_sample_splits: multi-valued D ---

test_that("get_sample_splits works with multi-valued D", {
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- X %*% runif(5) + rnorm(nobs)
  fun <- stepfun(quantile(D_tld, probs = c(0.25, 0.5, 0.75)),
                 c(1, 2, 3, 4))
  D <- fun(D_tld)

  res <- get_sample_splits(seq_len(nobs), sample_folds = 5,
                           cv_folds = 3, D = D)

  # 4 treatment levels
  expect_equal(length(res$subsamples_byD), 4)
  expect_equal(length(res$cv_subsamples_byD), 4)
  expect_equal(length(res$aux_indx), 4)

  # Subsamples partition 1:nobs
  all_indx <- sort(unlist(res$subsamples))
  expect_equal(all_indx, seq_len(nobs))
})#TEST_THAT

# --- get_sample_splits: clustered ---

test_that("get_sample_splits respects clusters", {
  nobs <- 500
  n_cluster <- 100
  cluster_variable <- sample(seq_len(n_cluster), nobs, replace = TRUE)
  D <- ifelse(rnorm(nobs) > 0, 1, 0)

  res <- get_sample_splits(cluster_variable, sample_folds = 5,
                           cv_folds = 3, D = D)

  # Clusters should be unique across folds
  cl_fold1 <- unique(cluster_variable[res$subsamples[[1]]])
  cl_fold2 <- unique(cluster_variable[res$subsamples[[2]]])
  expect_equal(length(intersect(cl_fold1, cl_fold2)), 0)
})#TEST_THAT

# --- get_sample_splits: pre-specified subsamples pass through ---

test_that("pre-specified subsamples pass through unchanged", {
  nobs <- 100
  cl <- seq_len(nobs)
  my_subsamples <- list(1:25, 26:50, 51:75, 76:100)

  res <- get_sample_splits(cl, subsamples = my_subsamples)

  for (k in seq_along(my_subsamples)) {
    expect_equal(res$subsamples[[k]], my_subsamples[[k]])
  }#FOR
})#TEST_THAT

# --- get_sample_splits: auto-merge subsamples_byD ---

test_that("subsamples_byD without subsamples auto-merges", {
  nobs <- 100
  D <- rep(c(0, 1), each = 50)
  cl <- seq_len(nobs)

  # Create by-D subsamples manually
  byD <- list(
    list(1:25, 26:50),   # D=0 folds (indices within D=0 subsample)
    list(1:25, 26:50)    # D=1 folds
  )

  res <- get_sample_splits(cl, D = D, subsamples_byD = byD)

  # subsamples should have been auto-merged
  expect_equal(length(res$subsamples), 2)
  all_indx <- sort(unlist(res$subsamples))
  expect_equal(all_indx, seq_len(nobs))
})#TEST_THAT

# --- get_cf_indices_stratified: balance ---

test_that("stratified folds have balanced treatment counts", {
  nobs <- 500
  D <- c(rep(0, 400), rep(1, 100))
  cl <- seq_len(nobs)

  res <- get_sample_splits(cl, sample_folds = 5, D = D,
                           stratify = TRUE)

  # Count D=1 observations in each fold
  d1_per_fold <- vapply(res$subsamples, function(idx) {
    sum(D[idx] == 1)
  }, FUN.VALUE = integer(1))

  # Folds should be roughly balanced (max - min <= 1 for integer rounding)
  expect_lte(max(d1_per_fold) - min(d1_per_fold), 1)
})#TEST_THAT

# --- check_subsamples: warning for small training sets ---

test_that("check_subsamples warns for small training sets", {
  # 50 obs split into 5 folds = 10 obs per fold = 40 training
  small_subsamples <- split(seq_len(50), rep(1:5, each = 10))
  names(small_subsamples) <- NULL

  expect_warning(
    check_subsamples(small_subsamples, NULL, stratify = FALSE),
    "only uses")
})#TEST_THAT

test_that("check_subsamples does not warn for large training sets", {
  large_subsamples <- split(seq_len(500), rep(1:5, each = 100))
  names(large_subsamples) <- NULL

  expect_silent(
    check_subsamples(large_subsamples, NULL, stratify = FALSE))
})#TEST_THAT

# --- auxiliary_X construction ---

test_that("auxiliary_X construction works with multi-valued D", {
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- X %*% runif(5) + rnorm(nobs)
  fun <- stepfun(quantile(D_tld, probs = c(0.25, 0.5, 0.75)),
                 c(1, 2, 3, 4))
  D <- fun(D_tld)

  res <- get_sample_splits(seq_len(nobs), sample_folds = 5, D = D)

  # Get auxiliary_X for first treatment level
  auxiliary_X_d <- get_auxiliary_X(res$aux_indx[[1]], X)

  # Check sizes
  expect_equal(
    length(unlist(res$aux_indx[[1]])),
    nobs - length(res$subsamples[[1]]) -
      length(res$subsamples_byD[[1]][[1]]))
  expect_equal(dim(auxiliary_X_d[[1]])[1],
               length(res$aux_indx[[1]][[1]]))
})#TEST_THAT
