# Tests for ral_rep base class =================================================

# Helper: construct a minimal ral_rep from synthetic data.
make_ral_rep <- function(n = 200, p = 2, R = 3) {
  fits <- lapply(seq_len(R), function(r) {
    set.seed(r * 100)
    coefficients <- matrix(rnorm(p), p, 1,
                            dimnames = list(paste0("x", 1:p),
                                            "fit1"))
    inf_func <- array(rnorm(n * p), dim = c(n, p, 1))
    dinf_dtheta <- array(1, dim = c(n, p, p, 1))
    ral(coefficients = coefficients,
        inf_func = inf_func,
        dinf_dtheta = dinf_dtheta,
        nobs = n,
        coef_names = paste0("x", 1:p),
        estimator_name = "Test RAL")
  })
  ral_rep(fits)
}#MAKE_RAL_REP

# Constructor ==================================================================

test_that("ral_rep() constructs valid object", {
  rr <- make_ral_rep()

  expect_s3_class(rr, "ral_rep")
  expect_equal(rr$nresamples, 3)
  expect_equal(rr$nobs, 200)
  expect_equal(rr$nfit, 1)
  expect_equal(rr$coef_names, c("x1", "x2"))
})#TEST_THAT

test_that("ral_rep() with subclass", {
  fits <- lapply(1:2, function(r) {
    set.seed(r)
    ral(coefficients = matrix(1, 1, 1),
        inf_func = array(rnorm(50), dim = c(50, 1, 1)),
        nobs = 50, coef_names = "x")
  })

  rr <- ral_rep(fits, subclass = "my_rep",
                custom_field = "hello")

  expect_equal(class(rr), c("my_rep", "ral_rep"))
  expect_s3_class(rr, "my_rep")
  expect_s3_class(rr, "ral_rep")
  expect_equal(rr$custom_field, "hello")
})#TEST_THAT

test_that("ral_rep() validates inputs", {
  n <- 50
  fit1 <- ral(coefficients = matrix(1, 1, 1),
              inf_func = array(rnorm(n), dim = c(n, 1, 1)),
              nobs = n, coef_names = "x")
  fit2 <- ral(coefficients = matrix(1, 1, 1),
              inf_func = array(rnorm(n), dim = c(n, 1, 1)),
              nobs = n, coef_names = "x")

  # Valid
  expect_no_error(ral_rep(list(fit1, fit2)))

  # < 2 fits
  expect_error(ral_rep(list(fit1)), "at least 2")

  # Non-ral element
  expect_error(ral_rep(list(fit1, "bad")),
               "does not inherit")

  # Mismatched nobs
  fit3 <- ral(coefficients = matrix(1, 1, 1),
              inf_func = array(rnorm(100), dim = c(100, 1, 1)),
              nobs = 100, coef_names = "x")
  expect_error(ral_rep(list(fit1, fit3)),
               "different 'nobs'")

  # Mismatched coef_names
  fit4 <- ral(coefficients = matrix(1, 1, 1),
              inf_func = array(rnorm(n), dim = c(n, 1, 1)),
              nobs = n, coef_names = "y")
  expect_error(ral_rep(list(fit1, fit4)),
               "different 'coef_names'")
})#TEST_THAT

# Accessors ====================================================================

test_that("[[ and length work", {
  rr <- make_ral_rep()

  expect_equal(length(rr), 3)
  expect_s3_class(rr[[1]], "ral")
  expect_s3_class(rr[[2]], "ral")
})#TEST_THAT

test_that("nobs.ral_rep works", {
  rr <- make_ral_rep()
  expect_equal(nobs(rr), 200)
})#TEST_THAT

# coef.ral_rep =================================================================

test_that("coef.ral_rep with median aggregation", {
  rr <- make_ral_rep()
  cf <- coef(rr)
  expect_true(is.numeric(cf))
  expect_length(cf, 2)

  # Median of per-rep coefs
  per_rep <- sapply(rr$fits, function(f)
    f$coefficients[, 1])
  expected <- apply(per_rep, 1, median)
  expect_equal(unname(cf), unname(expected),
               tolerance = 1e-10)
})#TEST_THAT

test_that("coef.ral_rep with mean aggregation", {
  rr <- make_ral_rep()
  cf <- coef(rr, aggregation = "mean")
  per_rep <- sapply(rr$fits, function(f)
    f$coefficients[, 1])
  expected <- apply(per_rep, 1, mean)
  expect_equal(unname(cf), unname(expected),
               tolerance = 1e-10)
})#TEST_THAT

# vcov.ral_rep =================================================================

test_that("vcov.ral_rep produces valid matrix", {
  rr <- make_ral_rep()

  V <- vcov(rr)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(2, 2))
  expect_true(all(diag(V) > 0))
  expect_equal(rownames(V), c("x1", "x2"))
})#TEST_THAT

test_that("vcov.ral_rep median matches manual computation", {
  rr <- make_ral_rep()
  R <- rr$nresamples
  p <- 2

  agg_coef <- coef(rr)
  V_arr <- array(NA, dim = c(p, p, R))
  for (r in seq_len(R)) {
    Sigma_r <- vcov(rr$fits[[r]])
    bdiff <- coef(rr$fits[[r]]) - agg_coef
    V_arr[, , r] <- Sigma_r + tcrossprod(bdiff)
  }#FOR
  V_manual <- apply(V_arr, c(1, 2), median)
  dimnames(V_manual) <- dimnames(vcov(rr))

  expect_equal(vcov(rr), V_manual, tolerance = 1e-10)
})#TEST_THAT

# confint.ral_rep ==============================================================

test_that("confint.ral_rep returns valid intervals", {
  rr <- make_ral_rep()

  ci <- confint(rr)
  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)
  expect_true(all(ci[, 1] < ci[, 2]))
  expect_identical(colnames(ci), c(" 2.5 %", "97.5 %"))
})#TEST_THAT

test_that("confint.ral_rep uniform wider than pointwise", {
  rr <- make_ral_rep()

  ci_pw <- confint(rr)
  set.seed(1)
  ci_uf <- confint(rr, uniform = TRUE, bootstraps = 499)

  width_pw <- ci_pw[, 2] - ci_pw[, 1]
  width_uf <- ci_uf[, 2] - ci_uf[, 1]
  expect_true(all(width_uf >= width_pw))
  expect_false(is.null(attr(ci_uf, "crit_val")))
})#TEST_THAT

# summary.ral_rep =============================================================

test_that("summary.ral_rep works", {
  rr <- make_ral_rep()

  s <- summary(rr)
  expect_s3_class(s, "summary.ral_rep")
  expect_true(is.array(s$coefficients))
  expect_equal(dim(s$coefficients), c(2, 4, 1))
  expect_equal(s$nobs, 200)
  expect_equal(s$nresamples, 3)
  expect_equal(s$aggregation, "median")

  # Print does not error
  out <- capture_output(print(s))
  expect_true(grepl("RAL estimation", out))
  expect_true(grepl("Resamples:", out))
  expect_true(grepl("median", out, ignore.case = TRUE))
})#TEST_THAT

test_that("summary.ral_rep with mean aggregation", {
  rr <- make_ral_rep()
  s <- summary(rr, aggregation = "mean")
  expect_equal(s$aggregation, "mean")
})#TEST_THAT

# tidy.ral_rep =================================================================

test_that("tidy.ral_rep returns valid data.frame", {
  rr <- make_ral_rep()

  td <- tidy(rr)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 2)
  expect_true(all(c("term", "estimate", "std.error",
                     "statistic", "p.value",
                     "fit_label", "aggregation")
    %in% names(td)))
  expect_equal(td$aggregation[1], "median")
})#TEST_THAT

test_that("tidy.ral_rep with conf.int", {
  rr <- make_ral_rep()
  td <- tidy(rr, conf.int = TRUE)

  expect_true(all(c("conf.low", "conf.high")
    %in% names(td)))
  expect_true(all(td$conf.low < td$conf.high))
})#TEST_THAT

# glance.ral_rep ===============================================================

test_that("glance.ral_rep returns one-row df", {
  rr <- make_ral_rep()
  gl <- glance(rr)

  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_equal(gl$nobs, 200)
  expect_equal(gl$nresamples, 3)
})#TEST_THAT

# print.ral_rep ================================================================

test_that("print.ral_rep works", {
  rr <- make_ral_rep()
  out <- capture_output(print(rr))
  expect_true(grepl("RAL replicated fits", out))
  expect_true(grepl("Resamples:", out))
})#TEST_THAT
