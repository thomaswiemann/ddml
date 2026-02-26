test_that("csolve returns a (generalized) inverse", {
  # Simulate small matrices, one of them singular
  X <- matrix(rnorm(30*30), 30, 30)
  X_sing <- rbind(cbind(X[1:15, 1:15], X[1:15, 1:15]),
                  cbind(X[1:15, 1:15], X[1:15, 1:15]))
  # Check combuted inverses
  expect_equal(csolve(X) %*% X, diag(1, 30))
  expect_equal((csolve(X_sing) %*% X_sing) %*% X_sing, X_sing)
})#TEST_THAT

test_that("is_single_learner detects correctly", {
  expect_true(is_single_learner(list(what = ols)))
  expect_true(is_single_learner(list(what = ols,
                                     args = list())))
  expect_false(is_single_learner(
    list(list(what = ols), list(what = ols))))
  expect_false(is_single_learner(
    list(list(fun = ols))))
})

test_that("normalize_learners resolves what and fun", {
  # what works
  l <- normalize_learners(list(what = ols))
  expect_identical(l$what, ols)

  # fun triggers deprecation; reset the flag first
  options(ddml.fun_deprecated_warned = NULL)
  expect_message(
    normalize_learners(list(list(fun = ols))),
    "deprecated")

  # what takes precedence over fun
  options(ddml.fun_deprecated_warned = NULL)
  l <- suppressMessages(
    normalize_learners(list(list(what = ols, fun = mdl_glmnet))))
  expect_identical(l[[1]]$what, ols)

  # neither → error
  expect_error(
    normalize_learners(list(list(args = list()))),
    "what.*fun")

  # reset flag
  options(ddml.fun_deprecated_warned = NULL)
})

test_that("ddml_plm works with list(fun=) for backward compat", {
  set.seed(42)
  nobs <- 100
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  # Reset deprecation flag
  options(ddml.fun_deprecated_warned = NULL)
  suppressWarnings(suppressMessages({
    fit <- ddml_plm(y, D, X,
                    learners = list(list(fun = ols),
                                    list(fun = ols)),
                    ensemble_type = "nnls1",
                    sample_folds = 2,
                    silent = TRUE)
  }))
  expect_s3_class(fit, "ddml")
  expect_true(!is.null(coef(fit)))
  options(ddml.fun_deprecated_warned = NULL)
})
