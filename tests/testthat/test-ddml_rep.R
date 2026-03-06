test_that("ddml_rep validates inputs correctly", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X[, 1] + rnorm(nobs)
  y <- 2 * D + X[, 2] + rnorm(nobs)

  fits <- lapply(1:3, function(r) {
    ddml_plm(y, D, X,
             learners = list(what = ols),
             sample_folds = 2,
             silent = TRUE)
  })

  # Valid construction works
  reps <- ddml_rep(fits)
  expect_s3_class(reps, "ddml_rep")
  expect_false(inherits(reps, "ddml"))
  expect_equal(reps$nresamples, 3)
  expect_equal(reps$model_type, "ddml_plm")
  expect_equal(reps$nobs, nobs)

  # List of length 1 fails
  expect_error(ddml_rep(fits[1]))

  # Non-ddml object fails
  expect_error(ddml_rep(list(fits[[1]], "not_ddml")))

  # Mismatched nobs fails
  fits_bad <- fits
  fits_bad[[2]]$nobs <- 999
  expect_error(ddml_rep(fits_bad))

  # Mismatched coef_names fails
  fits_bad2 <- fits
  fits_bad2[[2]]$coef_names <- "wrong_name"
  expect_error(ddml_rep(fits_bad2))
})

test_that("ddml_replicate produces valid ddml_rep object", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X[, 1] + rnorm(nobs)
  y <- 2 * D + X[, 2] + rnorm(nobs)

  reps <- ddml_replicate(ddml_plm,
                         y = y, D = D, X = X,
                         learners = list(what = ols),
                         sample_folds = 2,
                         resamples = 3,
                         silent = TRUE)

  expect_s3_class(reps, "ddml_rep")
  expect_equal(reps$nresamples, 3)
  expect_equal(reps$nobs, nobs)

  # Individual fits are accessible
  expect_s3_class(reps[[1]], "ddml")
  expect_s3_class(reps[[2]], "ddml_plm")
  expect_equal(length(reps), 3)

  # Different resamples produce different coefficients
  expect_false(identical(coef(reps[[1]]), coef(reps[[2]])))
})

test_that("ddml_replicate is reproducible with set.seed", {
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X[, 1] + rnorm(nobs)
  y <- 2 * D + X[, 2] + rnorm(nobs)

  set.seed(123)
  reps1 <- ddml_replicate(ddml_plm,
                          y = y, D = D, X = X,
                          learners = list(what = ols),
                          sample_folds = 2,
                          resamples = 3,
                          silent = TRUE)

  set.seed(123)
  reps2 <- ddml_replicate(ddml_plm,
                          y = y, D = D, X = X,
                          learners = list(what = ols),
                          sample_folds = 2,
                          resamples = 3,
                          silent = TRUE)

  expect_equal(coef(reps1), coef(reps2))
  expect_equal(coef(reps1[[1]]), coef(reps2[[1]]))
})

test_that("coef, vcov, confint work on ddml_rep", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X[, 1] + rnorm(nobs)
  y <- 2 * D + X[, 2] + rnorm(nobs)

  reps <- ddml_replicate(ddml_plm,
                         y = y, D = D, X = X,
                         learners = list(what = ols),
                         sample_folds = 2,
                         resamples = 3,
                         silent = TRUE)

  # coef (D + intercept = 2 elements)
  cf <- coef(reps)
  expect_true(is.numeric(cf))
  expect_length(cf, 2)

  # coef with mean aggregation differs from median
  cf_mean <- coef(reps, aggregation = "mean")
  expect_true(is.numeric(cf_mean))

  # Aggregated D coef is within range of per-resample D coefs
  per_rep_coefs <- sapply(seq_len(3),
    function(i) coef(reps[[i]])[1])
  expect_true(cf[1] >= min(per_rep_coefs) - 0.01)
  expect_true(cf[1] <= max(per_rep_coefs) + 0.01)

  # vcov
  V <- vcov(reps)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(2, 2))
  expect_true(V[1, 1] > 0)

  # vcov names
  expect_equal(rownames(V), reps$coef_names)
  expect_equal(colnames(V), reps$coef_names)

  # confint
  ci <- confint(reps)
  expect_true(is.matrix(ci))
  expect_equal(dim(ci), c(2, 2))
  expect_identical(colnames(ci), c(" 2.5 %", "97.5 %"))
  expect_true(ci[1, 1] < ci[1, 2])

  # confint at different level
  ci90 <- confint(reps, level = 0.90)
  expect_true((ci90[1, 2] - ci90[1, 1]) <
              (ci[1, 2] - ci[1, 1]))
})

test_that("summary and print work on ddml_rep", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  reps <- ddml_replicate(ddml_ate,
                         y = y, D = D, X = X,
                         learners = list(what = ols),
                         sample_folds = 2,
                         resamples = 3,
                         silent = TRUE)

  # summary
  s <- summary(reps)
  expect_s3_class(s, "summary.ddml_rep")
  expect_false(inherits(s, "summary.ddml"))
  expect_equal(s$nobs, nobs)
  expect_equal(s$nresamples, 3)
  expect_equal(s$aggregation, "median")
  expect_true(is.array(s$coefficients))
  expect_equal(dim(s$coefficients)[2], 4)

  # summary with mean aggregation
  s_mean <- summary(reps, aggregation = "mean")
  expect_equal(s_mean$aggregation, "mean")

  # summary with HC3
  s_hc3 <- summary(reps, type = "HC3")
  expect_equal(s_hc3$type, "HC3")

  # print.summary
  out <- capture_output(print(s))
  expect_true(grepl("DDML estimation", out))
  expect_true(grepl("Average Treatment Effect", out))
  expect_true(grepl("Resamples: 3", out))
  expect_true(grepl("median", out, ignore.case = TRUE))

  # print.ddml_rep
  out2 <- capture_output(print(reps))
  expect_true(grepl("replicated fits", out2,
                     ignore.case = TRUE))
  expect_true(grepl("Resamples: 3", out2))
})

test_that("tidy and glance work on ddml_rep", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  reps <- ddml_replicate(ddml_ate,
                         y = y, D = D, X = X,
                         learners = learners,
                         ensemble_type = c("ols", "nnls"),
                         cv_folds = 2,
                         sample_folds = 2,
                         resamples = 3,
                         silent = TRUE)

  # tidy single ensemble
  td <- tidy(reps, ensemble_idx = 1)
  expect_s3_class(td, "data.frame")
  expect_true(all(c("term", "estimate", "std.error",
    "statistic", "p.value", "ensemble_type",
    "aggregation") %in% colnames(td)))
  expect_equal(nrow(td), 1)
  expect_equal(td$aggregation, "median")

  # tidy all ensembles
  td_all <- tidy(reps, ensemble_idx = NULL)
  expect_equal(nrow(td_all), 2)

  # tidy with confidence intervals
  td_ci <- tidy(reps, conf.int = TRUE)
  expect_true(all(c("conf.low", "conf.high")
    %in% colnames(td_ci)))
  expect_true(td_ci$conf.low < td_ci$conf.high)

  # tidy with mean aggregation
  td_mean <- tidy(reps, aggregation = "mean")
  expect_s3_class(td_mean, "data.frame")
  expect_equal(td_mean$aggregation, "mean")

  # glance
  gl <- glance(reps)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_true("nresamples" %in% colnames(gl))
  expect_equal(gl$nresamples, 3)
  expect_equal(gl$model_type, "ddml_ate")
})

test_that("median aggregation formula is correct", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X[, 1] + rnorm(nobs)
  y <- 2 * D + X[, 2] + rnorm(nobs)

  R <- 5
  reps <- ddml_replicate(ddml_plm,
                         y = y, D = D, X = X,
                         learners = list(what = ols),
                         sample_folds = 2,
                         resamples = R,
                         silent = TRUE)

  coefs <- sapply(seq_len(R),
                   function(i) coef(reps[[i]]))
  vcovs <- lapply(seq_len(R),
                   function(i) vcov(reps[[i]]))
  p <- nrow(coefs)

  # theta_tilde = coordinate-wise median
  expected_coef <- apply(coefs, 1, median)

  # V_r = Sigma_r + (theta_r - theta_tilde)(...)^T
  V_list <- lapply(seq_len(R), function(r) {
    bdiff <- coefs[, r] - expected_coef
    vcovs[[r]] + tcrossprod(bdiff)
  })

  # tilde_Sigma_ij = median_r(V_r_ij)
  V_arr <- array(unlist(V_list), dim = c(p, p, R))
  expected_vcov <- apply(V_arr, c(1, 2), median)
  expected_se <- sqrt(diag(expected_vcov))

  actual_coef <- coef(reps, aggregation = "median")
  actual_vcov <- vcov(reps, aggregation = "median")
  s <- summary(reps, aggregation = "median")
  actual_se <- s$coefficients[, 2, 1]

  expect_equal(unname(actual_coef),
               unname(expected_coef),
               tolerance = 1e-10)
  expect_equal(unname(actual_vcov),
               unname(expected_vcov),
               tolerance = 1e-10)
  expect_equal(unname(actual_se),
               unname(expected_se),
               tolerance = 1e-10)
})

test_that("mean aggregation formula is correct", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X[, 1] + rnorm(nobs)
  y <- 2 * D + X[, 2] + rnorm(nobs)

  R <- 5
  reps <- ddml_replicate(ddml_plm,
                         y = y, D = D, X = X,
                         learners = list(what = ols),
                         sample_folds = 2,
                         resamples = R,
                         silent = TRUE)

  coefs <- sapply(seq_len(R),
                   function(i) coef(reps[[i]]))
  vcovs <- lapply(seq_len(R),
                   function(i) vcov(reps[[i]]))
  p <- nrow(coefs)

  # theta_tilde = arithmetic mean
  expected_coef <- rowMeans(coefs)

  # V_r = Sigma_r + (theta_r - theta_tilde)(...)^T
  V_list <- lapply(seq_len(R), function(r) {
    bdiff <- coefs[, r] - expected_coef
    vcovs[[r]] + tcrossprod(bdiff)
  })

  # tilde_Sigma = (1/R) sum_r V_r
  expected_vcov <- Reduce(`+`, V_list) / R
  expected_se <- sqrt(diag(expected_vcov))

  actual_coef <- coef(reps, aggregation = "mean")
  actual_vcov <- vcov(reps, aggregation = "mean")
  s <- summary(reps, aggregation = "mean")
  actual_se <- s$coefficients[, 2, 1]

  expect_equal(unname(actual_coef),
               unname(expected_coef),
               tolerance = 1e-10)
  expect_equal(unname(actual_vcov),
               unname(expected_vcov),
               tolerance = 1e-10)
  expect_equal(unname(actual_se),
               unname(expected_se),
               tolerance = 1e-10)
})

test_that("spectral aggregation works and equals median for p=1", {
  skip_if_not_installed("CVXR")

  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  reps <- ddml_replicate(ddml_ate,
                         y = y, D = D, X = X,
                         learners = list(what = ols),
                         sample_folds = 2,
                         resamples = 3,
                         silent = TRUE)

  # For p = 1, spectral reduces to scalar median
  cf_spec <- coef(reps, aggregation = "spectral")
  cf_med <- coef(reps, aggregation = "median")
  expect_equal(cf_spec, cf_med, tolerance = 1e-10)

  V_spec <- vcov(reps, aggregation = "spectral")
  V_med <- vcov(reps, aggregation = "median")
  expect_equal(V_spec, V_med, tolerance = 1e-10)

  s <- summary(reps, aggregation = "spectral")
  expect_equal(s$aggregation, "spectral")
})

test_that("spectral aggregation gives PSD matrix for p>1", {
  skip_if_not_installed("CVXR")

  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X[, 1] + rnorm(nobs)
  y <- 2 * D + X[, 2] + rnorm(nobs)

  reps <- ddml_replicate(ddml_plm,
                         y = y, D = D, X = X,
                         learners = list(what = ols),
                         sample_folds = 2,
                         resamples = 5,
                         silent = TRUE)

  V <- vcov(reps, aggregation = "spectral")
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(2, 2))

  # PSD: all eigenvalues >= 0
  eigs <- eigen(V, symmetric = TRUE, only.values = TRUE)
  expect_true(all(eigs$values >= -1e-10))

  # Symmetric
  expect_equal(V, t(V), tolerance = 1e-10)

  # Positive diagonal
  expect_true(all(diag(V) > 0))
})

test_that("ddml_rep works with ddml_ate", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  reps <- ddml_replicate(ddml_ate,
                         y = y, D = D, X = X,
                         learners = list(what = ols),
                         sample_folds = 2,
                         resamples = 3,
                         silent = TRUE)

  expect_s3_class(reps, "ddml_rep")
  expect_equal(reps$model_type, "ddml_ate")

  cf <- coef(reps)
  expect_length(cf, 1)

  s <- summary(reps)
  expect_s3_class(s, "summary.ddml_rep")
})

test_that("type argument threads through ddml_rep methods", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X[, 1] + rnorm(nobs)
  y <- 2 * D + X[, 2] + rnorm(nobs)

  reps <- ddml_replicate(ddml_plm,
                         y = y, D = D, X = X,
                         learners = list(what = ols),
                         sample_folds = 2,
                         resamples = 3,
                         silent = TRUE)

  V_hc0 <- vcov(reps, type = "HC0")
  V_hc1 <- vcov(reps, type = "HC1")
  V_hc3 <- vcov(reps, type = "HC3")

  # All positive
  expect_true(V_hc0[1, 1] > 0)
  expect_true(V_hc1[1, 1] > 0)
  expect_true(V_hc3[1, 1] > 0)

  # Coefficients identical across HC types
  expect_equal(coef(reps), coef(reps))

  # Summary stores type
  s3 <- summary(reps, type = "HC3")
  expect_equal(s3$type, "HC3")

  # HC3 shows in print
  out <- capture_output(print(s3))
  expect_true(grepl("HC3", out))
})
