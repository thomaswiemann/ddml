test_that("standard S3 generic methods work correctly", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_ate(y, D, X,
                    learners = learners,
                    ensemble_type = c("ols", "nnls"),
                    cv_folds = 2,
                    sample_folds = 2,
                    silent = TRUE)

  # coef
  cf <- coef(fit)
  expect_true(is.numeric(cf))
  expect_length(cf, 2)

  # vcov
  V <- vcov(fit, ensemble_idx = 2)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(1, 1))

  # confint
  ci <- confint(fit, ensemble_idx = 1)
  expect_true(is.matrix(ci))
  expect_equal(dim(ci), c(1, 2))
  expect_identical(colnames(ci), c(" 2.5 %", "97.5 %"))

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.ddml")
  expect_equal(s$nobs, nobs)

  # print.summary
  out <- capture_output(print(s))
  expect_true(grepl("DDML estimation", out))
  expect_true(grepl("Average Treatment Effect", out))
})

test_that("type argument threads through S3 methods", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_ate(y, D, X,
                    learners = learners,
                    ensemble_type = c("ols", "nnls"),
                    cv_folds = 2,
                    sample_folds = 2,
                    silent = TRUE)

  # vcov with all three types
  V_hc0 <- vcov(fit, type = "HC0")
  V_hc1 <- vcov(fit, type = "HC1")
  V_hc3 <- vcov(fit, type = "HC3")
  expect_true(V_hc1[1, 1] > V_hc0[1, 1])

  # summary stores type
  s0 <- summary(fit, type = "HC0")
  s1 <- summary(fit, type = "HC1")
  s3 <- summary(fit, type = "HC3")
  expect_equal(s0$type, "HC0")
  expect_equal(s1$type, "HC1")
  expect_equal(s3$type, "HC3")

  # HC0 SEs < HC1 SEs (dof correction)
  se_hc0 <- s0$inf_results[1, 2, 1]
  se_hc1 <- s1$inf_results[1, 2, 1]
  expect_true(se_hc1 > se_hc0)

  # Coefficients are identical across vcov types
  expect_equal(s0$inf_results[, 1, ], s1$inf_results[, 1, ])
  expect_equal(s0$inf_results[, 1, ], s3$inf_results[, 1, ])

  # confint with HC3
  ci_hc1 <- confint(fit, type = "HC1")
  ci_hc3 <- confint(fit, type = "HC3")
  expect_true(is.matrix(ci_hc3))
  expect_equal(dim(ci_hc3), c(1, 2))

  # HC3 print shows SE type
  out <- capture_output(print(s3))
  expect_true(grepl("HC3", out))

  # HC1 print does not show SE type (default)
  out1 <- capture_output(print(s1))
  expect_false(grepl("HC1", out1))
})

test_that("type works with PLM estimator", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  fit <- ddml_plm(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 2,
                  silent = TRUE)

  V_hc0 <- vcov(fit, type = "HC0")
  V_hc1 <- vcov(fit, type = "HC1")
  V_hc3 <- vcov(fit, type = "HC3")

  expect_true(is.matrix(V_hc3))
  expect_equal(dim(V_hc3), c(1, 1))
  expect_true(V_hc1[1, 1] > V_hc0[1, 1])
})
