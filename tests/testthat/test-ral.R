# Tests for ral base class =====================================================

test_that("ral() constructs a valid object", {
  n <- 50
  p <- 2

  coefficients <- matrix(c(1.5, -0.3), p, 1,
                          dimnames = list(c("a", "b"), "fit1"))
  inf_func <- array(rnorm(n * p), dim = c(n, p, 1))
  dinf_dtheta <- array(1, dim = c(n, p, p, 1))

  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             dinf_dtheta = dinf_dtheta,
             nobs = n,
             coef_names = c("a", "b"),
             estimator_name = "Test RAL")

  expect_s3_class(obj, "ral")
  expect_equal(class(obj), "ral")
  expect_equal(obj$nobs, n)
  expect_equal(obj$nfit, 1)
  expect_equal(obj$coef_names, c("a", "b"))
  expect_equal(obj$fit_labels, "fit1")
  expect_equal(obj$estimator_name, "Test RAL")
})#TEST_THAT

test_that("ral() with subclass", {
  n <- 30
  coefficients <- matrix(1.5, 1, 1,
                          dimnames = list("x", "e1"))
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))

  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n,
             coef_names = "x",
             subclass = "my_ral")

  expect_s3_class(obj, "my_ral")
  expect_s3_class(obj, "ral")
  expect_equal(class(obj), c("my_ral", "ral"))
})#TEST_THAT

test_that("ral() stores extra ... args", {
  n <- 30
  coefficients <- matrix(1, 1, 1,
                          dimnames = list("x", "e1"))
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))

  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = "x",
             custom_field = "hello",
             fixed_R = TRUE)

  expect_equal(obj$custom_field, "hello")
  expect_true(obj$fixed_R)
})#TEST_THAT

test_that("ral() validates inf_func dimensions", {
  n <- 30
  coefficients <- matrix(1, 2, 1)
  inf_func_bad <- array(rnorm(n * 3), dim = c(n, 3, 1))

  expect_error(ral(coefficients = coefficients,
                    inf_func = inf_func_bad,
                    nobs = n, coef_names = c("a", "b")),
               "inf_func")
})#TEST_THAT

test_that("ral() validates dinf_dtheta dimensions", {
  n <- 30
  coefficients <- matrix(1, 1, 1)
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))
  dinf_bad <- array(1, dim = c(n, 2, 2, 1))

  expect_error(ral(coefficients = coefficients,
                    inf_func = inf_func,
                    dinf_dtheta = dinf_bad,
                    nobs = n, coef_names = "x"),
               "dinf_dtheta")
})#TEST_THAT

# coef.ral =====================================================================

test_that("coef.ral returns named vector for single fit", {
  n <- 50
  coefficients <- matrix(c(1.5, -0.3), 2, 1,
                          dimnames = list(c("a", "b"), "e1"))
  inf_func <- array(rnorm(n * 2), dim = c(n, 2, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = c("a", "b"))

  cf <- coef(obj)
  expect_true(is.numeric(cf))
  expect_length(cf, 2)
  expect_equal(names(cf), c("a", "b"))
  expect_equal(cf[[1]], 1.5, tolerance = 1e-12)
})#TEST_THAT

test_that("coef.ral returns matrix for multiple fits", {
  n <- 50
  coefficients <- matrix(c(1, 2, 3, 4), 2, 2,
                          dimnames = list(c("a", "b"),
                                          c("e1", "e2")))
  inf_func <- array(rnorm(n * 2 * 2),
                     dim = c(n, 2, 2))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = c("a", "b"))

  cf <- coef(obj)
  expect_true(is.matrix(cf))
  expect_equal(dim(cf), c(2, 2))
})#TEST_THAT

# nobs.ral =====================================================================

test_that("nobs.ral returns correct count", {
  n <- 42
  coefficients <- matrix(1, 1, 1)
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = "x")

  expect_equal(nobs(obj), 42)
})#TEST_THAT

# hatvalues.ral ================================================================

test_that("hatvalues.ral returns leverage when dinf_dtheta present", {
  n <- 50
  coefficients <- matrix(1, 1, 1)
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))
  dinf_dtheta <- array(runif(n, 0.5, 1.5),
                        dim = c(n, 1, 1, 1))

  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             dinf_dtheta = dinf_dtheta,
             nobs = n, coef_names = "x")

  h <- hatvalues(obj)
  expect_length(h, n)
  expect_true(all(is.finite(h)))
})#TEST_THAT

test_that("hatvalues.ral warns when dinf_dtheta is NULL", {
  n <- 50
  coefficients <- matrix(1, 1, 1)
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = "x")

  expect_warning(hatvalues(obj), "dinf_dtheta")
})#TEST_THAT

# vcov.ral =====================================================================

test_that("vcov.ral produces valid covariance matrix", {
  n <- 200
  p <- 2
  coefficients <- matrix(c(1, 2), p, 1,
                          dimnames = list(c("a", "b"), "e1"))
  inf_func <- array(rnorm(n * p), dim = c(n, p, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = c("a", "b"))

  V <- vcov(obj)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(p, p))
  expect_true(all(diag(V) > 0))
  expect_equal(rownames(V), c("a", "b"))
  expect_equal(colnames(V), c("a", "b"))

  # PSD
  evals <- eigen(V, only.values = TRUE)$values
  expect_true(all(evals >= -1e-10))
})#TEST_THAT

test_that("vcov.ral HC0 vs HC1 differ by dof correction", {
  n <- 100
  coefficients <- matrix(1, 1, 1)
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = "x")

  V0 <- vcov(obj, type = "HC0")
  V1 <- vcov(obj, type = "HC1")
  expect_equal(V1[1, 1], V0[1, 1] * n / (n - 1),
               tolerance = 1e-10)
})#TEST_THAT

test_that("vcov.ral HC3 uses leverage correction", {
  n <- 100
  coefficients <- matrix(1, 1, 1)
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))
  dinf_dtheta <- array(1, dim = c(n, 1, 1, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             dinf_dtheta = dinf_dtheta,
             nobs = n, coef_names = "x")

  V3 <- vcov(obj, type = "HC3")
  V0 <- vcov(obj, type = "HC0")
  # HC3 >= HC0 in general
  expect_true(V3[1, 1] >= V0[1, 1] - 1e-12)
})#TEST_THAT

# validate_fit_idx =============================================================

test_that("fit_idx out of range errors", {
  n <- 50
  coefficients <- matrix(1, 1, 1)
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = "x")

  expect_error(vcov(obj, fit_idx = 0), "fit_idx")
  expect_error(vcov(obj, fit_idx = 2), "fit_idx")
  expect_error(vcov(obj, fit_idx = 99), "fit_idx")
  expect_error(confint(obj, fit_idx = 99), "fit_idx")
})#TEST_THAT

# confint.ral ==================================================================

test_that("confint.ral returns valid intervals", {
  n <- 200
  coefficients <- matrix(c(1, 2), 2, 1,
                          dimnames = list(c("a", "b"), "e1"))
  inf_func <- array(rnorm(n * 2), dim = c(n, 2, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = c("a", "b"))

  ci <- confint(obj)
  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)
  expect_true(all(ci[, 1] < ci[, 2]))
  expect_identical(colnames(ci), c(" 2.5 %", "97.5 %"))
  expect_equal(attr(ci, "crit_val"), qnorm(0.975),
               tolerance = 1e-6)
})#TEST_THAT

test_that("confint.ral parm subsetting works", {
  n <- 200
  coefficients <- matrix(c(1, 2, 3), 3, 1,
                          dimnames = list(c("a", "b", "c"),
                                          "e1"))
  inf_func <- array(rnorm(n * 3), dim = c(n, 3, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = c("a", "b", "c"))

  # By name
  ci_ab <- confint(obj, parm = c("a", "b"))
  expect_equal(nrow(ci_ab), 2)
  expect_equal(rownames(ci_ab), c("a", "b"))

  # By index
  ci_2 <- confint(obj, parm = 2)
  expect_equal(nrow(ci_2), 1)
  expect_equal(rownames(ci_2), "b")
})#TEST_THAT

test_that("confint.ral uniform wider than pointwise", {
  set.seed(42)
  n <- 200
  p <- 3
  coefficients <- matrix(rep(0, p), p, 1,
                          dimnames = list(letters[1:p], "e1"))
  inf_func <- array(rnorm(n * p), dim = c(n, p, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n,
             coef_names = letters[1:p])

  ci_pw <- confint(obj)
  set.seed(1)
  ci_uf <- confint(obj, uniform = TRUE, bootstraps = 499)

  width_pw <- ci_pw[, 2] - ci_pw[, 1]
  width_uf <- ci_uf[, 2] - ci_uf[, 1]
  expect_true(all(width_uf >= width_pw))
  expect_true(attr(ci_uf, "crit_val") > qnorm(0.975))
})#TEST_THAT

# summary.ral =================================================================

test_that("summary.ral produces valid output", {
  n <- 200
  coefficients <- matrix(c(1, 2), 2, 1,
                          dimnames = list(c("a", "b"), "e1"))
  inf_func <- array(rnorm(n * 2), dim = c(n, 2, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = c("a", "b"))

  s <- summary(obj)
  expect_s3_class(s, "summary.ral")
  expect_true(is.array(s$coefficients))
  expect_equal(dim(s$coefficients), c(2, 4, 1))
  expect_equal(s$nobs, n)
  expect_equal(s$type, "HC1")

  # Print does not error
  out <- capture_output(print(s))
  expect_true(grepl("RAL estimation", out))
  expect_true(grepl("Obs:", out))
})#TEST_THAT

# tidy.ral =====================================================================

test_that("tidy.ral returns valid data.frame", {
  n <- 200
  coefficients <- matrix(c(1, 2), 2, 1,
                          dimnames = list(c("a", "b"), "e1"))
  inf_func <- array(rnorm(n * 2), dim = c(n, 2, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = c("a", "b"))

  td <- tidy(obj)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 2)
  expect_true(all(c("term", "estimate", "std.error",
                     "statistic", "p.value",
                     "fit_label") %in% names(td)))
  expect_equal(td$term, c("a", "b"))
})#TEST_THAT

test_that("tidy.ral with conf.int", {
  n <- 200
  coefficients <- matrix(c(1, 2), 2, 1,
                          dimnames = list(c("a", "b"), "e1"))
  inf_func <- array(rnorm(n * 2), dim = c(n, 2, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = c("a", "b"))

  td <- tidy(obj, conf.int = TRUE)
  expect_true(all(c("conf.low", "conf.high") %in%
                    names(td)))
  expect_true(all(td$conf.low < td$conf.high))
})#TEST_THAT

test_that("tidy.ral fit_idx=NULL returns all fits", {
  n <- 50
  coefficients <- matrix(c(1, 2, 3, 4), 2, 2,
                          dimnames = list(c("a", "b"),
                                          c("e1", "e2")))
  inf_func <- array(rnorm(n * 2 * 2),
                     dim = c(n, 2, 2))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = c("a", "b"))

  td_all <- tidy(obj, fit_idx = NULL)
  expect_equal(nrow(td_all), 4)  # 2 params x 2 fits
})#TEST_THAT

# glance.ral ==================================================================

test_that("glance.ral returns one-row data.frame", {
  n <- 50
  coefficients <- matrix(1, 1, 1)
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))
  obj <- ral(coefficients = coefficients,
             inf_func = inf_func,
             nobs = n, coef_names = "x",
             estimator_name = "Test")

  gl <- glance(obj)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_equal(gl$nobs, n)
  expect_equal(gl$estimator_name, "Test")
})#TEST_THAT

# plot.ral =====================================================================

test_that("plot.ral produces a plot without error", {
  n <- 100
  p <- 3
  inf <- array(stats::rnorm(n * p), c(n, p, 1))
  theta <- matrix(c(0.5, -0.3, 0.1), p, 1)
  obj <- ral(theta, inf, nobs = n,
             coef_names = c("b1", "b2", "b3"))

  # Plot should run without error
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  res <- plot(obj)

  expect_type(res, "list")
  expect_named(res, c("coefficients", "ci", "labels"))
  expect_equal(length(res$coefficients), p)
  expect_equal(nrow(res$ci), p)
})#TEST_THAT

test_that("plot.ral with uniform = TRUE produces wider bands", {
  set.seed(123)
  n <- 200
  p <- 3
  inf <- array(stats::rnorm(n * p), c(n, p, 1))
  theta <- matrix(c(0.5, -0.3, 0.1), p, 1)
  obj <- ral(theta, inf, nobs = n,
             coef_names = c("b1", "b2", "b3"))

  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)

  res_pw <- plot(obj, uniform = FALSE)
  res_uf <- plot(obj, uniform = TRUE)

  # Uniform bands should be at least as wide as pointwise
  width_pw <- res_pw$ci[, 2] - res_pw$ci[, 1]
  width_uf <- res_uf$ci[, 2] - res_uf$ci[, 1]
  expect_true(all(width_uf >= width_pw - 1e-10))
})#TEST_THAT

test_that("plot.ral with parm selects subset", {
  n <- 100
  p <- 3
  inf <- array(stats::rnorm(n * p), c(n, p, 1))
  theta <- matrix(c(0.5, -0.3, 0.1), p, 1)
  obj <- ral(theta, inf, nobs = n,
             coef_names = c("b1", "b2", "b3"))

  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  res <- plot(obj, parm = c("b1", "b3"))

  expect_equal(length(res$coefficients), 2)
  expect_equal(res$labels, c("b1", "b3"))
})#TEST_THAT
