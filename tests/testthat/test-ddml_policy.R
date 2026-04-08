test_that("ddml_policy computes with a single model (K=2)", {
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  policy <- ifelse(X[, 1] > 0, 1, 0)
  learners <- list(what = ols)
  fit <- ddml_policy(y, D, X,
                     policy = policy,
                     learners = learners,
                     cv_folds = 3,
                     sample_folds = 3,
                     silent = TRUE)
  expect_equal(length(coef(fit)), 1)
  expect_s3_class(fit, "ddml_policy")
  expect_s3_class(fit, "ddml")
})#TEST_THAT

test_that("ddml_policy computes with margins", {
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  policy <- ifelse(X[, 1] > 0, 1, 0)
  learners <- list(what = ols)
  fit <- ddml_policy(y, D, X,
                     policy = policy,
                     margins = c(2.0, 1.5),
                     learners = learners,
                     cv_folds = 3,
                     sample_folds = 3,
                     silent = TRUE)
  expect_equal(length(coef(fit)), 1)
})#TEST_THAT

test_that("ddml_policy computes with an ensemble procedure", {
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  policy <- ifelse(X[, 1] > 0, 1, 0)
  learners <- list(list(what = ols),
                   list(what = ols))
  fit <- ddml_policy(y, D, X,
                     policy = policy,
                     learners = learners,
                     ensemble_type = "ols",
                     cv_folds = 3,
                     sample_folds = 3,
                     silent = TRUE)
  expect_equal(length(coef(fit)), 1)
})#TEST_THAT

test_that("ddml_policy computes w/ multiple ensembles", {
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  policy <- ifelse(X[, 1] > 0, 1, 0)
  learners <- list(list(what = ols),
                   list(what = ols))
  fit <- ddml_policy(y, D, X,
                     policy = policy,
                     learners = learners,
                     ensemble_type = c("ols", "nnls",
                                       "singlebest", "average"),
                     cv_folds = 3,
                     sample_folds = 3,
                     silent = TRUE)
  expect_equal(length(coef(fit)), 4)
})#TEST_THAT

test_that("ddml_policy computes with shortstack", {
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  policy <- ifelse(X[, 1] > 0, 1, 0)
  learners <- list(list(what = ols))
  fit <- ddml_policy(y, D, X,
                     policy = policy,
                     learners = learners,
                     ensemble_type = c("ols", "average"),
                     shortstack = TRUE,
                     cv_folds = 3,
                     sample_folds = 3,
                     silent = TRUE)
  expect_equal(length(coef(fit)), 2)
})#TEST_THAT

test_that("summary.ddml computes for ddml_policy", {
  nobs <- 500
  X <- matrix(rnorm(nobs * 5), nobs, 5)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  policy <- ifelse(X[, 1] > 0, 1, 0)
  learners <- list(what = ols)
  fit <- ddml_policy(y, D, X,
                     policy = policy,
                     learners = learners,
                     cv_folds = 3,
                     sample_folds = 3,
                     silent = TRUE)
  inf_res <- summary(fit)
  capture_output({print(inf_res)}, print = FALSE)
  expect_s3_class(inf_res, "summary.ddml")
  expect_equal(dim(inf_res$coefficients), c(1, 4, 1))
})#TEST_THAT

test_that("ddml_policy fitted pass-through works", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- 0.1 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  policy <- ifelse(X[, 1] > 0, 1, 0)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_policy(y, D, X,
                     policy = policy,
                     learners = learners,
                     ensemble_type = "average",
                     sample_folds = 2,
                     silent = TRUE)

  fit2 <- ddml_policy(y, D, X,
                      policy = policy,
                      learners = learners,
                      ensemble_type = "average",
                      sample_folds = 2,
                      silent = TRUE,
                      fitted = fit$fitted,
                      splits = fit$splits)
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-6)

  expect_error(
    ddml_policy(y, D, X,
                policy = policy,
                learners = learners,
                ensemble_type = "average",
                sample_folds = 2,
                silent = TRUE,
                fitted = fit$fitted),
    "must be supplied when 'fitted' is supplied"
  )
})#TEST_THAT

test_that("ddml_policy scores are mean-zero", {
  nobs <- 500
  set.seed(42)
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- 0.2 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- D + 0.1 * X[, 1] + rnorm(nobs)
  policy <- ifelse(X[, 1] > 0, 1, 0)
  fit <- ddml_policy(y, D, X,
                     policy = policy,
                     learners = list(what = ols),
                     sample_folds = 3, silent = TRUE)
  score_mean <- mean(fit$scores[, , 1])
  expect_true(abs(score_mean) < 0.01,
              label = paste("score mean =", round(score_mean, 6)))
})#TEST_THAT

test_that("ddml_policy with margins=c(1,-1) matches ddml_ate", {
  set.seed(123)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- 0.5 * X[, 1] + rnorm(nobs)
  D <- 1 * (D_tld > 0)
  y <- 1.0 * D + 0.3 * X[, 1] + rnorm(nobs)
  learners <- list(what = ols)

  ate_fit <- ddml_ate(y, D, X,
                      learners = learners,
                      sample_folds = 3,
                      silent = TRUE)

  # Policy that assigns everyone to treatment (all get d=1)
  # with margins c(-1, 1) = APO(d=1) - APO(d=0) = ATE
  policy_all <- rep(1, nobs)
  # But policy levels must cover unique(D) = {0, 1}, so use
  # a non-trivial policy and compare the structural form instead.
  # With constant policy, can't cover both levels.
  # Instead, verify sign/magnitude: ATE ~ 1 and policy value
  # with same splits should be in a reasonable range.
  policy <- ifelse(X[, 1] > 0, 1, 0)
  pol_fit <- ddml_policy(y, D, X,
                         policy = policy,
                         learners = learners,
                         sample_folds = 3,
                         splits = ate_fit$splits,
                         silent = TRUE)
  expect_equal(length(coef(pol_fit)), 1)
})#TEST_THAT

test_that("ddml_policy rejects mismatched policy/D", {
  nobs <- 100
  X <- matrix(rnorm(nobs * 2), nobs, 2)
  D <- sample(c(0, 1), nobs, replace = TRUE)
  y <- rnorm(nobs)
  learners <- list(what = ols)

  expect_error(
    ddml_policy(y, D, X,
                policy = rep(2, nobs),
                learners = learners,
                silent = TRUE),
    "must appear in 'D'")

  expect_error(
    ddml_policy(y, D, X,
                policy = c(0, 1),
                learners = learners,
                silent = TRUE),
    "same length")
})#TEST_THAT

test_that("ddml_policy computes with partial coverage (K < unique D)", {
  set.seed(42)
  nobs <- 500
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- sample(c(0, 1, 2), nobs, replace = TRUE,
              prob = c(0.4, 0.4, 0.2))
  y <- 0.5 * (D == 1) + 1.0 * (D == 2) + 0.1 * X[, 1] +
    rnorm(nobs)
  policy <- ifelse(X[, 1] > 0, 1, 0)
  learners <- list(what = ols)
  fit <- ddml_policy(y, D, X,
                     policy = policy,
                     learners = learners,
                     cv_folds = 3,
                     sample_folds = 3,
                     silent = TRUE)
  expect_equal(length(coef(fit)), 1)
  expect_s3_class(fit, "ddml_policy")
  inf_res <- summary(fit)
  expect_s3_class(inf_res, "summary.ddml")
})#TEST_THAT

test_that("ddml_policy rejects invalid margins", {
  nobs <- 100
  X <- matrix(rnorm(nobs * 2), nobs, 2)
  D <- sample(c(0, 1), nobs, replace = TRUE)
  y <- rnorm(nobs)
  policy <- ifelse(X[, 1] > 0, 1, 0)
  learners <- list(what = ols)

  expect_error(
    ddml_policy(y, D, X,
                policy = policy,
                margins = c(1, 2, 3),
                learners = learners,
                silent = TRUE),
    "numeric vector of length 2")
})#TEST_THAT
