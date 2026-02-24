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
