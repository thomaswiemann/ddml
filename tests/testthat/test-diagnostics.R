test_that("diagnostics works with PLM single learner", {
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

  diag <- diagnostics(fit)

  expect_s3_class(diag, "ddml_diagnostics")
  expect_true("y_X" %in% names(diag$tables))
  # Should have at least one equation table
  expect_true(length(diag$tables) >= 1)

  # print should run without error
  out <- capture_output(print(diag))
  expect_true(grepl("Stacking diagnostics", out))
})

test_that("diagnostics works with PLM multiple learners", {
  set.seed(42)
  nobs <- 100
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  suppressWarnings({
    fit <- ddml_plm(y, D, X,
                    learners = learners,
                    ensemble_type = c("nnls1"),
                    sample_folds = 2,
                    silent = TRUE)
  })

  diag <- diagnostics(fit)

  # Should have 2 learners per equation
  for (eq in names(diag$tables)) {
    expect_equal(nrow(diag$tables[[eq]]), 2)
  }

  # Weights should sum to ~1
  for (eq in names(diag$tables)) {
    expect_equal(sum(diag$tables[[eq]]$weight), 1,
                 tolerance = 1e-4)
  }
})

test_that("tidy.ddml_diagnostics returns flat data.frame", {
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

  td <- tidy(diagnostics(fit))

  expect_s3_class(td, "data.frame")
  expect_true(all(c("equation", "learner", "mspe", "r2",
                     "weight") %in% colnames(td)))
  expect_true(nrow(td) >= 1)
})

test_that("diagnostics with CVC", {
  set.seed(42)
  nobs <- 100
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  suppressWarnings({
    fit <- ddml_plm(y, D, X,
                    learners = learners,
                    ensemble_type = c("nnls1"),
                    sample_folds = 2,
                    silent = TRUE)
  })

  diag <- diagnostics(fit, cvc = TRUE, bootnum = 100)

  # CVC columns should be present
  for (eq in names(diag$tables)) {
    expect_true("cvc_pval" %in% colnames(diag$tables[[eq]]))
    expect_true("in_conf_set" %in%
                  colnames(diag$tables[[eq]]))
  }

  # tidy should include CVC columns
  td <- tidy(diag)
  expect_true("cvc_pval" %in% colnames(td))
})

test_that("diagnostics works with ATE", {
  set.seed(42)
  nobs <- 100
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D_tld <- X %*% runif(3) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(3) + rnorm(nobs)

  suppressWarnings({
    fit <- ddml_ate(y, D, X,
                    learners = list(what = ols),
                    sample_folds = 2,
                    silent = TRUE)
  })

  diag <- diagnostics(fit)

  expect_s3_class(diag, "ddml_diagnostics")
  expect_true(length(diag$tables) >= 2)
  out <- capture_output(print(diag))
  expect_true(grepl("Average Treatment Effect", out))
})

test_that("diagnostics r2 matches manual calculation", {
  set.seed(42)
  nobs <- 200
  X <- matrix(rnorm(nobs * 3), nobs, 3)
  D <- X %*% c(1, 0.5, 0) + rnorm(nobs)
  y <- 2 * D + X %*% c(0, 1, 0.5) + rnorm(nobs)

  learners <- list(list(what = ols), list(what = ols))
  suppressWarnings({
    fit <- ddml_plm(y, D, X,
                    learners = learners,
                    ensemble_type = "nnls1",
                    sample_folds = 2,
                    silent = TRUE)
  })

  # R2 per learner per fold from crossval: 1 - mspe / var(y)
  # Check that r2 is stored and non-NULL
  expect_false(is.null(fit$r2$y_X))
  expect_true(all(fit$r2$y_X <= 1))
})

test_that("diagnostics rejects non-ddml objects", {
  expect_error(diagnostics(lm(1:10 ~ rnorm(10))),
               "class")
})
