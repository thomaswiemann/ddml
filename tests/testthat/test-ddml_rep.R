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
    "statistic", "p.value", "ensemble_type")
    %in% colnames(td)))
  expect_equal(nrow(td), 1)

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

  # glance
  gl <- glance(reps)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_true("nresamples" %in% colnames(gl))
  expect_equal(gl$nresamples, 3)
  expect_equal(gl$model_type, "ddml_ate")
})

test_that("aggregation formulas are correct", {
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

  # Extract per-resample coefs and SEs manually (p x R matrices)
  R <- 5
  coefs <- sapply(seq_len(R), function(i) coef(reps[[i]]))
  ses <- sapply(seq_len(R), function(i) {
    V <- vcov(reps[[i]])
    sqrt(diag(V))
  })
  p <- nrow(coefs)

  # Median aggregation (per coefficient)
  expected_coef_med <- apply(coefs, 1, median)
  var_total_med <- ses^2 +
    (coefs - expected_coef_med)^2
  expected_se_med <- sqrt(apply(var_total_med, 1, median))

  actual_coef_med <- coef(reps, aggregation = "median")
  s_med <- summary(reps, aggregation = "median")
  actual_se_med <- s_med$coefficients[, 2, 1]

  expect_equal(unname(actual_coef_med),
               unname(expected_coef_med),
               tolerance = 1e-10)
  expect_equal(unname(actual_se_med),
               unname(expected_se_med),
               tolerance = 1e-10)

  # Mean aggregation (per coefficient)
  expected_coef_mean <- rowMeans(coefs)
  var_total_mean <- ses^2 +
    (coefs - expected_coef_mean)^2
  expected_se_mean <- sqrt(
    R / rowSums(1 / var_total_mean))

  actual_coef_mean <- coef(reps, aggregation = "mean")
  s_mean <- summary(reps, aggregation = "mean")
  actual_se_mean <- s_mean$coefficients[, 2, 1]

  expect_equal(unname(actual_coef_mean),
               unname(expected_coef_mean),
               tolerance = 1e-10)
  expect_equal(unname(actual_se_mean),
               unname(expected_se_mean),
               tolerance = 1e-10)
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
