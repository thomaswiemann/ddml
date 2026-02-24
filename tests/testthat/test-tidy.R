test_that("tidy and glance broom formatters work correctly", {
  # Generate a minimal fitted object to test the plumbing
  set.seed(42)
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs * 4), nobs, 4))
  D_tld <- X %*% runif(5) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(5) + rnorm(nobs)
  
  learners <- list(list(fun = ols), list(fun = ols))
  suppressWarnings({
    fit <- ddml_ate(y, D, X,
                    learners = learners,
                    ensemble_type = c("ols", "nnls"),
                    cv_folds = 2,
                    sample_folds = 2,
                    silent = TRUE)
  })
  
  # tidy for a single ensemble
  td1 <- tidy(fit, ensemble_idx = 1)
  expect_s3_class(td1, "data.frame")
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value", "ensemble_type") %in% colnames(td1)))
  expect_equal(nrow(td1), 1)
  expect_equal(td1$ensemble_type, "ols")
  
  # tidy for all ensembles
  td_all <- tidy(fit, ensemble_idx = NULL)
  expect_s3_class(td_all, "data.frame")
  expect_equal(nrow(td_all), 2)
  
  # tidy with confidence intervals
  td_ci <- tidy(fit, ensemble_idx = 2, conf.int = TRUE)
  expect_s3_class(td_ci, "data.frame")
  expect_true(all(c("conf.low", "conf.high") %in% colnames(td_ci)))
  expect_true(td_ci$conf.low < td_ci$conf.high)
  expect_equal(td_ci$ensemble_type, "nnls")
  
  # glance
  gl <- glance(fit)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_true(all(c("nobs", "sample_folds", "shortstack", "ensemble_type", "model_type") %in% colnames(gl)))
})

test_that("tidy respects type argument", {
  set.seed(42)
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs * 4), nobs, 4))
  D_tld <- X %*% runif(5) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(5) + rnorm(nobs)

  learners <- list(list(fun = ols), list(fun = ols))
  suppressWarnings({
    fit <- ddml_ate(y, D, X,
                    learners = learners,
                    ensemble_type = c("ols"),
                    cv_folds = 2,
                    sample_folds = 2,
                    silent = TRUE)
  })

  td_hc0 <- tidy(fit, type = "HC0")
  td_hc1 <- tidy(fit, type = "HC1")
  td_hc3 <- tidy(fit, type = "HC3")

  # Estimates are identical across HC types
  expect_equal(td_hc0$estimate, td_hc1$estimate)
  expect_equal(td_hc0$estimate, td_hc3$estimate)

  # HC1 SE > HC0 SE (dof correction)
  expect_true(td_hc1$std.error > td_hc0$std.error)

  # Confidence intervals widen with HC3
  ci_hc1 <- tidy(fit, type = "HC1", conf.int = TRUE)
  ci_hc3 <- tidy(fit, type = "HC3", conf.int = TRUE)
  width_hc1 <- ci_hc1$conf.high - ci_hc1$conf.low
  width_hc3 <- ci_hc3$conf.high - ci_hc3$conf.low
  expect_true(width_hc3 >= width_hc1)
})

test_that("tidy works with PLM (multi-covariate D)", {
  set.seed(42)
  nobs <- 100
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  suppressWarnings({
    fit <- ddml_plm(y, D, X,
                    learners = list(what = ols),
                    sample_folds = 2,
                    silent = TRUE)
  })

  td <- tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 1)
  expect_true("term" %in% colnames(td))

  gl <- glance(fit)
  expect_equal(gl$model_type, "ddml_plm")
  expect_equal(gl$nobs, nobs)
})
