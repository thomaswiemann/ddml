test_that("ddml() constructs a valid ddml object", {
  # Minimal valid inputs
  n <- 50
  y <- rnorm(n)
  theta <- mean(y)

  coef <- matrix(theta, 1, 1, dimnames = list("mean", "custom"))
  scores <- array(y - theta, dim = c(n, 1, 1))
  J <- array(-1, dim = c(1, 1, 1))
  inf_func <- array(y - theta, dim = c(n, 1, 1))
  dinf_dtheta <- array(1, dim = c(n, 1, 1, 1))

  fit <- ddml(
    coefficients = coef, scores = scores, J = J,
    inf_func = inf_func, dinf_dtheta = dinf_dtheta,
    nobs = n, coef_names = "mean",
    estimator_name = "Sample Mean",
    sample_folds = 5,
    ensemble_weights = list(y = matrix(1, 1, 1,
      dimnames = list("learner1", "custom"))))

  # Class

  expect_s3_class(fit, "ddml")

  # S3 methods work
  expect_no_error(summary(fit))
  expect_no_error(confint(fit))
  expect_no_error(coef(fit))
  expect_no_error(vcov(fit))
  expect_no_error(nobs(fit))
  expect_no_error(tidy(fit))
  expect_no_error(glance(fit))
})

test_that("ddml() with subclass", {
  n <- 30
  coef <- matrix(1.5, 1, 1, dimnames = list("param", "ens"))
  scores <- array(rnorm(n), dim = c(n, 1, 1))
  J <- array(-1, dim = c(1, 1, 1))
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))
  dinf_dtheta <- array(1, dim = c(n, 1, 1, 1))

  fit <- ddml(
    coefficients = coef, scores = scores, J = J,
    inf_func = inf_func, dinf_dtheta = dinf_dtheta,
    nobs = n, coef_names = "param",
    estimator_name = "My Custom",
    subclass = "my_custom")

  expect_s3_class(fit, "my_custom")
  expect_s3_class(fit, "ddml")
  expect_equal(class(fit), c("my_custom", "ddml", "ral"))
})

test_that("ddml() passes extra args via ...", {
  n <- 30
  coef <- matrix(1, 1, 1, dimnames = list("x", "ens"))
  scores <- array(rnorm(n), dim = c(n, 1, 1))
  J <- array(-1, dim = c(1, 1, 1))
  inf_func <- array(rnorm(n), dim = c(n, 1, 1))
  dinf_dtheta <- array(1, dim = c(n, 1, 1, 1))

  fit <- ddml(
    coefficients = coef, scores = scores, J = J,
    inf_func = inf_func, dinf_dtheta = dinf_dtheta,
    nobs = n, coef_names = "x",
    estimator_name = "Test",
    my_extra = "hello",
    learners = list(what = ols))

  expect_equal(fit$my_extra, "hello")
  expect_equal(fit$learners$what, ols)
})

test_that("ddml() rejects bad inputs", {
  n <- 20

  # coefficients not a matrix
  expect_error(
    ddml(coefficients = 1:3, scores = array(1, c(n, 3, 1)),
         J = array(1, c(3, 3, 1)),
         inf_func = array(1, c(n, 3, 1)),
         nobs = n, coef_names = "x",
         estimator_name = "Bad"),
    "matrix")

  # scores not 3D
  expect_error(
    ddml(coefficients = matrix(1, 1, 1),
         scores = matrix(1, n, 1),
         J = array(1, c(1, 1, 1)),
         inf_func = array(1, c(n, 1, 1)),
         nobs = n, coef_names = "x",
         estimator_name = "Bad"),
    "3D array")

  # J not 3D
  expect_error(
    ddml(coefficients = matrix(1, 1, 1),
         scores = array(1, c(n, 1, 1)),
         J = matrix(1, 1, 1),
         inf_func = array(1, c(n, 1, 1)),
         nobs = n, coef_names = "x",
         estimator_name = "Bad"),
    "3D array")

  # inf_func not 3D
  expect_error(
    ddml(coefficients = matrix(1, 1, 1),
         scores = array(1, c(n, 1, 1)),
         J = array(1, c(1, 1, 1)),
         inf_func = matrix(1, n, 1),
         nobs = n, coef_names = "x",
         estimator_name = "Bad"),
    "3D numeric array")

  # dinf_dtheta not 4D array
  expect_error(
    ddml(coefficients = matrix(1, 1, 1),
         scores = array(1, c(n, 1, 1)),
         J = array(1, c(1, 1, 1)),
         inf_func = array(1, c(n, 1, 1)),
         dinf_dtheta = 123,
         nobs = n, coef_names = "x",
         estimator_name = "Bad"),
    "dinf_dtheta.*4D array")

  # Dimension mismatch: scores wrong nobs
  expect_error(
    ddml(coefficients = matrix(1, 1, 1),
         scores = array(1, c(n + 5, 1, 1)),
         J = array(1, c(1, 1, 1)),
         inf_func = array(1, c(n, 1, 1)),
         nobs = n, coef_names = "x",
         estimator_name = "Bad"),
    "nobs x p x nensb")

  # Dimension mismatch: J wrong shape
  expect_error(
    ddml(coefficients = matrix(1, 2, 1),
         scores = array(1, c(n, 2, 1)),
         J = array(1, c(1, 1, 1)),
         inf_func = array(1, c(n, 2, 1)),
         nobs = n, coef_names = c("a", "b"),
         estimator_name = "Bad"),
    "p x p x nensb")

  # Dimension mismatch: inf_func wrong wrong length/shape
  expect_error(
    ddml(coefficients = matrix(1, 1, 2),
         scores = array(1, c(n, 1, 2)),
         J = array(1, c(1, 1, 2)),
         inf_func = array(1, c(n, 1, 1)),
         nobs = n, coef_names = "x",
         estimator_name = "Bad"),
    "nobs x p x nensb")
  
  # Dimension mismatch: dinf_dtheta wrong dims
  expect_error(
    ddml(coefficients = matrix(1, 1, 2),
         scores = array(1, c(n, 1, 2)),
         J = array(1, c(1, 1, 2)),
         inf_func = array(1, c(n, 1, 2)),
         dinf_dtheta = array(1, dim = c(n, 1, 1, 3)),
         nobs = n, coef_names = "x",
         estimator_name = "Bad"),
    "'dinf_dtheta' dimensions must be")
})

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
  V <- vcov(fit, fit_idx = 2)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(1, 1))

  # confint
  ci <- confint(fit, fit_idx = 1)
  expect_true(is.matrix(ci))
  expect_equal(dim(ci), c(1, 2))
  expect_identical(colnames(ci), c(" 2.5 %", "97.5 %"))

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.ddml")
  expect_equal(s$nobs, nobs)

  # nobs
  expect_equal(nobs(fit), nobs)

  # call
  expect_false(is.null(fit$call))

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
  se_hc0 <- s0$coefficients[1, 2, 1]
  se_hc1 <- s1$coefficients[1, 2, 1]
  expect_true(se_hc1 > se_hc0)

  # Coefficients are identical across vcov types
  expect_equal(s0$coefficients[, 1, ], s1$coefficients[, 1, ])
  expect_equal(s0$coefficients[, 1, ], s3$coefficients[, 1, ])

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
  expect_equal(dim(V_hc3), c(2, 2))
  expect_true(V_hc1[1, 1] > V_hc0[1, 1])
})

test_that("vcov rejects invalid type", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- 1 * (rnorm(nobs) > 0)
  y <- D + rnorm(nobs)

  fit <- ddml_ate(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 2, silent = TRUE)

  expect_error(vcov(fit, type = "HC2"))
  expect_error(vcov(fit, type = "hc1"))
})

test_that("confint parm subsetting works", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- cbind(D1 = rnorm(nobs), D2 = rnorm(nobs))
  y <- D[, 1] + D[, 2] + rnorm(nobs)

  fit <- ddml_plm(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 2,
                  silent = TRUE)

  ci_all <- confint(fit)
  expect_equal(nrow(ci_all), 3)

  ci_d1 <- confint(fit, parm = "D1")
  expect_equal(nrow(ci_d1), 1)
  expect_equal(rownames(ci_d1), "D1")

  ci_num <- confint(fit, parm = 1:2)
  expect_equal(nrow(ci_num), 2)

  expect_error(confint(fit, parm = "nonexistent"))
})

test_that("hatvalues returns correct length", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- 1 * (rnorm(nobs) > 0)
  y <- D + rnorm(nobs)

  fit <- ddml_ate(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 2, silent = TRUE)

  h <- hatvalues(fit)
  expect_length(h, nobs)
  expect_true(is.numeric(h))
  expect_true(all(is.finite(h)))
})

test_that("hatvalues works for multi-parameter estimator", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + rnorm(nobs)

  fit <- ddml_plm(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 2,
                  silent = TRUE)

  h <- hatvalues(fit)
  expect_length(h, nobs)
  expect_true(all(is.finite(h)))
})

test_that("fit_idx out of range errors", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- 1 * (rnorm(nobs) > 0)
  y <- D + rnorm(nobs)

  fit <- ddml_ate(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 2, silent = TRUE)

  expect_error(vcov(fit, fit_idx = 99))
  expect_error(hatvalues(fit, fit_idx = 0))
})

test_that("tidy and glance broom formatters work correctly", {
  # Generate a minimal fitted object to test the plumbing
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
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_ate(y, D, X,
                    learners = learners,
                    ensemble_type = c("ols"),
                    cv_folds = 2,
                    sample_folds = 2,
                    silent = TRUE)

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
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  fit <- ddml_plm(y, D, X,
                  learners = list(what = ols),
                  sample_folds = 2,
                  silent = TRUE)

  td <- tidy(fit)
  expect_s3_class(td, "data.frame")
  expect_equal(nrow(td), 2)
  expect_true("term" %in% colnames(td))

  gl <- glance(fit)
  expect_equal(gl$model_type, "ddml_plm")
  expect_equal(gl$nobs, nobs)
})

test_that("confint uniform produces wider bands", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- cbind(D1 = rnorm(nobs), D2 = rnorm(nobs))
  y <- D[, 1] + D[, 2] + rnorm(nobs)

  fit <- ddml_plm(y, D, X,
                   learners = list(what = ols),
                   sample_folds = 2,
                   silent = TRUE)

  ci_pw <- confint(fit)

  set.seed(1)
  ci_uf <- confint(fit, uniform = TRUE,
                    bootstraps = 499)

  # Uniform bands are wider (p > 1)
  width_pw <- ci_pw[, 2] - ci_pw[, 1]
  width_uf <- ci_uf[, 2] - ci_uf[, 1]
  expect_true(all(width_uf >= width_pw))

  # Critical value attribute
  expect_false(is.null(attr(ci_uf, "crit_val")))
  expect_true(attr(ci_uf, "crit_val") > qnorm(0.975))

  # Pointwise has Gaussian quantile
  expect_equal(attr(ci_pw, "crit_val"), qnorm(0.975), tolerance = 1e-6)
})

test_that("confint uniform with scalar estimator", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + rnorm(nobs)

  fit <- ddml_plm(y, D, X,
                   learners = list(what = ols),
                   sample_folds = 2,
                   silent = TRUE)

  ci_pw <- confint(fit)
  set.seed(1)
  ci_uf <- confint(fit, uniform = TRUE,
                    bootstraps = 499)

  # For p=1, uniform ~ pointwise (within tolerance)
  width_pw <- as.numeric(ci_pw[, 2] - ci_pw[, 1])
  width_uf <- as.numeric(ci_uf[, 2] - ci_uf[, 1])
  expect_true(all(abs(width_uf - width_pw) / width_pw
    < 0.25))
})

test_that("tidy uniform produces wider conf.int", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- cbind(D1 = rnorm(nobs), D2 = rnorm(nobs))
  y <- D[, 1] + D[, 2] + rnorm(nobs)

  fit <- ddml_plm(y, D, X,
                   learners = list(what = ols),
                   sample_folds = 2,
                   silent = TRUE)

  td_pw <- tidy(fit, conf.int = TRUE)

  set.seed(1)
  td_uf <- tidy(fit, conf.int = TRUE,
                 uniform = TRUE, bootstraps = 499)

  # Uniform conf intervals are wider
  width_pw <- td_pw$conf.high - td_pw$conf.low
  width_uf <- td_uf$conf.high - td_uf$conf.low
  expect_true(all(width_uf >= width_pw))
})

test_that("confint uniform produces wider bands for ddml_rep", {
  set.seed(42)
  nobs <- 300
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- cbind(D1 = rnorm(nobs), D2 = rnorm(nobs))
  y <- D[, 1] + D[, 2] + rnorm(nobs)

  reps <- ddml_replicate(ddml_plm, y = y, D = D, X = X,
                          learners = list(what = ols),
                          sample_folds = 2,
                          resamples = 3, silent = TRUE)

  ci_pw <- confint(reps)

  set.seed(1)
  ci_uf <- confint(reps, uniform = TRUE, bootstraps = 499)

  # Uniform bands are wider (p > 1)
  width_pw <- ci_pw[, 2] - ci_pw[, 1]
  width_uf <- ci_uf[, 2] - ci_uf[, 1]
  expect_true(all(width_uf >= width_pw))

  # Critical value attribute
  expect_false(is.null(attr(ci_uf, "crit_val")))
  expect_true(attr(ci_uf, "crit_val") > qnorm(0.975))

  # Pointwise has Gaussian quantile
  expect_equal(attr(ci_pw, "crit_val"), qnorm(0.975), tolerance = 1e-6)
})

test_that("hatvalues warns and returns NA if dinf_dtheta is missing", {
  fit <- ddml(
    coefficients = matrix(0, 1, 1, dimnames = list("ATE", "nnls")),
    scores = array(0, dim = c(10, 1, 1)),
    J = array(1, dim = c(1, 1, 1)),
    inf_func = array(runif(10), dim = c(10, 1, 1)),
    dinf_dtheta = NULL,
    nobs = 10,
    coef_names = "ATE",
    estimator_name = "test",
    ensemble_type = "nnls",
    cluster_variable = seq_len(10)
  )

  expect_warning(
    expect_equal(hatvalues(fit), rep(NA_real_, 10)),
    "dinf_dtheta not available",
    ignore.case = TRUE
  )
})

# as.list.ddml ================================================================

test_that("as.list.ddml splits by ensemble", {
  set.seed(42)
  y <- AE98[1:500, "worked"]
  D <- AE98[1:500, "morekids"]
  X <- AE98[1:500, c("age", "agefst", "black")]

  fit <- ddml_plm(y, D, X,
                  learners = list(
                    list(what = ols),
                    list(what = ols, args = list(const = FALSE))),
                  ensemble_type = c("nnls", "singlebest"),
                  sample_folds = 2, silent = TRUE)

  L <- as.list(fit)
  expect_length(L, 2)
  expect_equal(names(L), c("nnls", "singlebest"))
  for (j in 1:2) {
    expect_s3_class(L[[j]], "ddml")
    expect_equal(ncol(L[[j]]$coefficients), 1L)
    # tidy and glance work
    expect_true(is.data.frame(tidy(L[[j]])))
    expect_true(is.data.frame(glance(L[[j]])))
  }
  # coef values match parent
  expect_equal(L[[1]]$coefficients[, 1],
               fit$coefficients[, 1], tolerance = 1e-10)
  expect_equal(L[[2]]$coefficients[, 1],
               fit$coefficients[, 2], tolerance = 1e-10)
})#TEST_THAT
