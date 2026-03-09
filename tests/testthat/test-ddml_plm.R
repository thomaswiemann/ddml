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
                           silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_plm_fit)), 2)
})#TEST_THAT

test_that("ddml_plm computes with clustered observations", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  cluster_variable <- sample(1:100, nobs, replace = TRUE)
  # Define arguments
  learners <- list(what = mdl_glmnet,
                   args = list(alpha = 0.5))
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           cv_folds = 3,
                           sample_folds = 3,
                           cluster_variable = cluster_variable,
                           silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_plm_fit)), 2)
})#TEST_THAT

test_that("ddml_plm computes with an ensemble procedure", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           ensemble_type = "ols",
                           shortstack = FALSE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_plm_fit)), 2)
})#TEST_THAT

test_that("ddml_plm computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = FALSE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_plm_fit)), 10)
})#TEST_THAT

test_that("ddml_plm computes with multiple ensemble procedures & sparse mats", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, as(X, "sparseMatrix"),
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = FALSE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_plm_fit)), 10)
})#TEST_THAT
test_that("ddml_plm computes w/ an ensemble procedure & shortstacking", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           ensemble_type = "ols",
                           shortstack = TRUE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_plm_fit)), 2)
})#TEST_THAT

test_that("ddml_plm computes w/ multiple ensemble procedures & shortstacking", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = TRUE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_plm_fit)), 10)
})#TEST_THAT

test_that("ddml_plm computes w/ ensemble procedures & custom weights", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = ols),
                   list(what = ols))
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
                           silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_plm_fit)), 14)
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
                           silent = TRUE)
  inf_res <- summary(ddml_plm_fit)
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(2, 4, 1))
})#TEST_THAT

test_that("summary.ddml_plm computes with a single model and dependence", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  cluster_variable <- sample(1:100, nobs, replace = TRUE)
  # Define arguments
  learners <- list(what = mdl_glmnet,
                   args = list(alpha = 0.5))
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           cv_folds = 3,
                           sample_folds = 3,
                           cluster_variable = cluster_variable,
                           silent = TRUE)
  inf_res <- summary(ddml_plm_fit)
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(2, 4, 1))
})#TEST_THAT

test_that("summary.ddml_plm computes with multiple ensemble procedures", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(what = ols))
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
                           silent = TRUE)
  inf_res <- summary(ddml_plm_fit)
  capture_output(print(inf_res), print = FALSE)
  # Check output with expectations
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(2, 4, 7))
})#TEST_THAT

test_that("ddml_plm computes with an ensemble procedure and multivariate D", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  cbind(X %*% runif(40) + rnorm(nobs), rnorm(nobs))
  y <- rowSums(D) + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners = learners,
                           ensemble_type = "ols",
                           shortstack = FALSE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_plm_fit)), 3)
})#TEST_THAT

test_that("ddml_plm computes with multiple ensemble types and multivariate D", {
  # Simulate small dataset
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D <-  cbind(X %*% runif(40) + rnorm(nobs), rnorm(nobs))
  y <- rowSums(D) + X %*% runif(40) + rnorm(nobs)
  # Define arguments
  learners <- list(list(what = mdl_glmnet,
                        args = list(alpha = 0.5)),
                   list(what = ols))
  # Compute DDML PLM estimator
  ddml_plm_fit <- ddml_plm(y, D, X,
                           learners,
                           ensemble_type = c("ols", "nnls",
                                             "nnls1",
                                             "singlebest", "average"),
                           shortstack = FALSE,
                           cv_folds = 3,
                           sample_folds = 3,
                           silent = TRUE)
  # Check output with expectations
  expect_equal(length(coef(ddml_plm_fit)), 15)
})#TEST_THAT

test_that("ddml_plm backward-compat cv_subsamples_list works with message", {
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs * 5), nobs, 5))
  D <- rnorm(nobs)
  y <- D + X %*% runif(6) + rnorm(nobs)
  # Pre-generate splits
  splits <- get_sample_splits(seq_len(nobs),
                              sample_folds = 3, cv_folds = 3)
  learners <- list(list(what = ols), list(what = ols))
  # Call with deprecated cv_subsamples_list — should emit message + warning
  expect_warning(
    expect_message(
      ddml_plm(y, D, X,
               learners = learners,
               ensemble_type = "ols",
               sample_folds = 3,
               splits = list(subsamples = splits$subsamples),
               cv_subsamples_list = splits$cv_subsamples,
               silent = TRUE),
      "cv_subsamples_list has been renamed"),
    "Deprecated split arguments detected")
})#TEST_THAT

test_that("ddml_plm computes with parallel", {
  skip_on_cran()
  skip_if_not_installed("parallel")
  set.seed(42)
  nobs <- 200
  X <- cbind(1, matrix(rnorm(nobs * 39), nobs, 39))
  D <- X %*% runif(40) + rnorm(nobs)
  y <- D + X %*% runif(40) + rnorm(nobs)
  learners <- list(what = ols)
  splits <- get_sample_splits(seq_len(nobs), sample_folds = 3)
  # Sequential
  res_seq <- ddml_plm(y, D, X,
                      learners = learners,
                      sample_folds = 3,
                      splits = list(subsamples = splits$subsamples),
                      silent = TRUE)
  # Parallel
  res_par <- ddml_plm(y, D, X,
                      learners = learners,
                      sample_folds = 3,
                      splits = list(subsamples = splits$subsamples),
                      silent = TRUE,
                      parallel = list(cores = 2))
  expect_equal(coef(res_par), coef(res_seq))
})#TEST_THAT

test_that("ddml_plm fitted pass-through works", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  presplits <- get_sample_splits(seq_len(nobs),
                                 sample_folds = 2,
                                 cv_folds = 2)
  fit <- ddml_plm(y, D, X,
                  learners = learners,
                  ensemble_type = "average",
                  sample_folds = 2,
                  splits = list(
                    subsamples = presplits$subsamples,
                    cv_subsamples = presplits$cv_subsamples),
                  silent = TRUE)

  # fitted should be stored with per-equation crossfit data
  expect_true(!is.null(fit$fitted))
  expect_true(!is.null(fit$fitted$y_X$cf_fitted_bylearner))
  expect_true(!is.null(fit$fitted$y_X$cf_resid_bylearner))

  # Pass-through with average ensemble reproduces exactly
  fit2 <- ddml_plm(y, D, X,
                   learners = learners,
                   ensemble_type = "average",
                   sample_folds = 2,
                   silent = TRUE,
                   fitted = fit$fitted,
                   splits = fit$splits)
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-6)

  # Pass-through with different ensemble gives valid result
  fit_nnls <- ddml_plm(y, D, X,
                       learners = learners,
                       ensemble_type = "nnls1",
                       sample_folds = 2,
                       silent = TRUE,
                       fitted = fit$fitted,
                       splits = fit$splits)
  expect_s3_class(fit_nnls, "ddml")
  expect_true(is.numeric(coef(fit_nnls)))

  expect_error(
    ddml_plm(y, D, X,
             learners = learners,
             ensemble_type = "average",
             sample_folds = 2,
             silent = TRUE,
             fitted = fit$fitted),
    "must be supplied when 'fitted' is supplied"
  )
})

test_that("ddml_plm pass-through with nnls/nnls1 reproduces exactly via cv_resid_byfold", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  presplits <- get_sample_splits(seq_len(nobs),
                                 sample_folds = 2,
                                 cv_folds = 2)
  splits <- list(
    subsamples = presplits$subsamples,
    cv_subsamples = presplits$cv_subsamples)

  for (ens in c("nnls", "nnls1")) {
    fit_fresh <- ddml_plm(y, D, X,
                          learners = learners,
                          ensemble_type = ens,
                          sample_folds = 2,
                          splits = splits,
                          silent = TRUE)

    # cv_resid_byfold must be stored for exact reproduction
    expect_true(
      !is.null(fit_fresh$fitted$y_X$cv_resid_byfold),
      info = paste(ens, "y_X cv_resid_byfold"))

    # Pass-through with same ensemble reproduces exactly
    fit_pt <- ddml_plm(y, D, X,
                       learners = learners,
                       ensemble_type = ens,
                       sample_folds = 2,
                       splits = splits,
                       silent = TRUE,
                       fitted = fit_fresh$fitted)
    expect_equal(coef(fit_pt), coef(fit_fresh),
                 tolerance = 1e-10,
                 info = paste(ens, "exact reproduction"))
  }
})

test_that("ddml_plm pass-through with save_crossval = FALSE uses approximate path", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  presplits <- get_sample_splits(seq_len(nobs),
                                 sample_folds = 2,
                                 cv_folds = 2)
  splits <- list(
    subsamples = presplits$subsamples,
    cv_subsamples = presplits$cv_subsamples)

  # Fit with save_crossval = FALSE strips cv_resid_byfold
  fit_no_cv <- ddml_plm(y, D, X,
                        learners = learners,
                        ensemble_type = "nnls1",
                        sample_folds = 2,
                        splits = splits,
                        save_crossval = FALSE,
                        silent = TRUE)
  expect_null(fit_no_cv$fitted$y_X$cv_resid_byfold)

  # average ensemble: still exact without cv_resid_byfold
  fit_avg <- ddml_plm(y, D, X,
                      learners = learners,
                      ensemble_type = "average",
                      sample_folds = 2,
                      splits = splits,
                      silent = TRUE)
  fit_avg_pt <- ddml_plm(y, D, X,
                         learners = learners,
                         ensemble_type = "average",
                         sample_folds = 2,
                         splits = splits,
                         silent = TRUE,
                         fitted = fit_no_cv$fitted)
  expect_equal(coef(fit_avg_pt), coef(fit_avg),
               tolerance = 1e-10)

  # nnls1 ensemble: approximate (close but not identical)
  fit_nnls_fresh <- ddml_plm(y, D, X,
                             learners = learners,
                             ensemble_type = "nnls1",
                             sample_folds = 2,
                             splits = splits,
                             silent = TRUE)
  fit_nnls_pt <- ddml_plm(y, D, X,
                          learners = learners,
                          ensemble_type = "nnls1",
                          sample_folds = 2,
                          splits = splits,
                          silent = TRUE,
                          fitted = fit_no_cv$fitted)
  expect_true(is.numeric(coef(fit_nnls_pt)))
  expect_equal(coef(fit_nnls_pt), coef(fit_nnls_fresh),
               tolerance = 0.5)
})

test_that("ddml_plm HC0/HC1/HC3 SEs close to sandwich::vcovHC on ols_fit", {
  skip_if_not_installed("sandwich")
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D <- X %*% c(1, 0.5, 0, 0, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5, 0.3, 0) + rnorm(nobs)

  fit <- ddml_plm(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 5,
                  silent = TRUE)

  # Reconstruct the final partialing-out regression since ols_fit is removed:
  y_r <- as.vector(y - fit$fitted$y_X$cf_fitted_bylearner[, 1])
  D_r <- as.matrix(D - fit$fitted$D_X[[1]]$cf_fitted_bylearner[, 1])
  colnames(D_r) <- colnames(D)
  ols_fit <- stats::lm(y_r ~ D_r)

  for (type in c("HC0", "HC1")) {
    V_ddml <- vcov(fit, type = type)
    V_sw <- sandwich::vcovHC(ols_fit, type = type)
    # Reorder sandwich: intercept first -> last
    idx <- c(seq_len(nrow(V_sw))[-1], 1)
    V_sw_reord <- V_sw[idx, idx, drop = FALSE]
    expect_equal(as.numeric(V_ddml),
                 as.numeric(V_sw_reord),
                 tolerance = 1e-8,
                 info = paste("PLM", type))
  }
  # HC3: regressor-based leverage matches sandwich exactly
  V_ddml <- vcov(fit, type = "HC3")
  V_sw <- sandwich::vcovHC(ols_fit, type = "HC3")
  idx <- c(seq_len(nrow(V_sw))[-1], 1)
  V_sw_reord <- V_sw[idx, idx, drop = FALSE]
  expect_equal(as.numeric(V_ddml),
               as.numeric(V_sw_reord),
               tolerance = 1e-8,
               info = "PLM HC3")
})#TEST_THAT

test_that("ddml_plm legacy split args warn and still work", {
  set.seed(101)
  nobs <- 200
  X <- matrix(rnorm(nobs * 4), nobs, 4)
  D <- X %*% c(1, 0.5, 0.2, 0) + rnorm(nobs)
  y <- D + X %*% c(0.3, 0.1, 0.4, 0.2) + rnorm(nobs)
  learners <- list(what = ols)
  splits <- get_sample_splits(seq_len(nobs),
                              sample_folds = 3,
                              cv_folds = 3)
  fit_splits <- ddml_plm(
    y, D, X,
    learners = learners,
    sample_folds = 3,
    cv_folds = 3,
    splits = list(
      subsamples = splits$subsamples,
      cv_subsamples = splits$cv_subsamples
    ),
    silent = TRUE
  )
  expect_warning(
    fit_legacy <- ddml_plm(
      y, D, X,
      learners = learners,
      sample_folds = 3,
      cv_folds = 3,
      subsamples = splits$subsamples,
      cv_subsamples = splits$cv_subsamples,
      silent = TRUE
    ),
    "Deprecated split arguments detected"
  )
  expect_equal(coef(fit_legacy), coef(fit_splits), tolerance = 1e-8)
})
