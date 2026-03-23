# Tests for ddml_attgt -------------------------------------------------------
# Follows the exact pattern from test-ddml_att.R

test_that("ddml_attgt computes with a single learner", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)
  for (i in seq_len(n)) {
    if (is.finite(G[i])) {
      for (tt in seq_len(T_)) {
        if (tt >= G[i]) y[i, tt] <- y[i, tt] + 1
      }
    }
  }

  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                  learners = list(what = ols),
                  sample_folds = 3, cv_folds = 3,
                  silent = TRUE)

  # Class and structure
  expect_s3_class(fit, "ddml_attgt")
  expect_s3_class(fit, "ddml")

  # Coefficient count = number of (g,t) cells
  expect_true(length(coef(fit)) > 0)

  # All cell info populated
  expect_equal(nrow(fit$cell_info), length(coef(fit)))
  expect_true(all(fit$cell_info$n_treated > 0))
  expect_true(all(fit$cell_info$n_control > 0))

  # Post-treatment cells: true effect is 1.0
  post <- fit$cell_info$time >= fit$cell_info$group
  expect_true(all(abs(coef(fit)[post] - 1) < 0.5))
  # Pre-treatment cells: placebo, true effect is 0
  if (any(!post)) {
    expect_true(all(abs(coef(fit)[!post]) < 0.5))
  }#IF
})#TEST_THAT

test_that("ddml_attgt computes with multiple ensemble types", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)

  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                  learners = list(list(what = ols),
                                  list(what = ols)),
                  ensemble_type = c("ols", "average"),
                  sample_folds = 2, cv_folds = 3,
                  silent = TRUE)

  # Multiple ensembles: coefficients is a matrix
  expect_true(is.matrix(fit$coefficients))
  expect_equal(ncol(fit$coefficients), 2)
})#TEST_THAT

test_that("summary.ddml works for ddml_attgt", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, Inf), n, replace = TRUE,
              prob = c(0.4, 0.6))
  y <- matrix(rnorm(n * T_), n, T_)

  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                  learners = list(what = ols),
                  sample_folds = 2, cv_folds = 3,
                  silent = TRUE)

  inf_res <- summary(fit)
  capture_output({print(inf_res)}, print = FALSE)
  expect_s3_class(inf_res, "summary.ddml")
  # Rows = number of cells, cols = 4 (est, se, z, p)
  expect_equal(dim(inf_res$coefficients)[2], 4)
})#TEST_THAT

test_that("vcov.ddml returns correct dimensions for ddml_attgt", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)

  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                  learners = list(what = ols),
                  sample_folds = 2, cv_folds = 3,
                  silent = TRUE)

  V <- vcov(fit)
  C_ <- length(coef(fit))
  expect_equal(dim(V), c(C_, C_))
  # Diagonal must be positive
  expect_true(all(diag(V) > 0))
})#TEST_THAT

test_that("ddml_attgt control_group = 'nevertreated' works", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)

  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                  learners = list(what = ols),
                  control_group = "nevertreated",
                  sample_folds = 2, cv_folds = 3,
                  silent = TRUE)

  expect_s3_class(fit, "ddml_attgt")
  # All control counts must equal number of never-treated
  n_never <- sum(is.infinite(G))
  for (idx in seq_len(nrow(fit$cell_info))) {
    expect_equal(fit$cell_info$n_control[idx], n_never)
  }
})#TEST_THAT

test_that("ddml_attgt cross-fitting consistency across cells", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)

  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                  learners = list(what = ols),
                  sample_folds = 3, cv_folds = 3,
                  silent = TRUE)

  # Structure: cell_info has right dimensions
  expect_equal(nrow(fit$cell_info), length(coef(fit)))
  expect_true(all(fit$cell_info$n_treated > 0))
})#TEST_THAT

test_that("ddml_attgt with X = NULL works", {
  set.seed(42)
  n <- 800; T_ <- 4
  G <- sample(c(3, Inf), n, replace = TRUE,
              prob = c(0.4, 0.6))
  y <- matrix(rnorm(n * T_), n, T_)

  fit <- ddml_attgt(y, X = NULL, t = 1:T_, G = G,
                  learners = list(what = ols),
                  sample_folds = 2, cv_folds = 3,
                  silent = TRUE)

  expect_s3_class(fit, "ddml_attgt")
})#TEST_THAT

test_that("ddml_attgt diagnostics works", {
  set.seed(42)
  n <- 800; T_ <- 3
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(2, Inf), n, replace = TRUE,
              prob = c(0.4, 0.6))
  y <- matrix(rnorm(n * T_), n, T_)

  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                  learners = list(list(what = ols),
                                  list(what = ols)),
                  ensemble_type = "ols",
                  sample_folds = 2, cv_folds = 3,
                  silent = TRUE)

  d <- diagnostics(fit)
  expect_s3_class(d, "ddml_diagnostics")
  # Should have 3 equations (y_X_D0, D_X, D) per cell
  C_ <- nrow(fit$cell_info)
  expect_equal(length(d$tables), 3 * C_)
  # All equation names follow the cell:eq pattern
  expect_true(all(grepl("^ATT\\(", names(d$tables))))
  capture_output({print(d)}, print = FALSE)
})#TEST_THAT

test_that("ddml_attgt fitted pass-through works", {
  set.seed(42)
  n <- 800; T_ <- 3
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(2, Inf), n, replace = TRUE,
              prob = c(0.4, 0.6))
  y <- matrix(rnorm(n * T_), n, T_)

  learners <- list(list(what = ols), list(what = ols))
  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                  learners = learners,
                  ensemble_type = "average",
                  sample_folds = 2,
                  silent = TRUE)

  fit2 <- ddml_attgt(y, X, t = 1:T_, G = G,
                   learners = learners,
                   ensemble_type = "average",
                   sample_folds = 2,
                   silent = TRUE,
                   fitted = fit$fitted,
                   splits = fit$splits)
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-6)

  # Error when fitted provided without splits
  expect_error(
    ddml_attgt(y, X, t = 1:T_, G = G,
             learners = learners,
             ensemble_type = "average",
             sample_folds = 2,
             silent = TRUE,
             fitted = fit$fitted),
    "must be supplied when 'fitted' is supplied"
  )
})#TEST_THAT

test_that("confint.ddml works for ddml_attgt", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, Inf), n, replace = TRUE,
              prob = c(0.4, 0.6))
  y <- matrix(rnorm(n * T_), n, T_)

  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                  learners = list(what = ols),
                  sample_folds = 2, cv_folds = 3,
                  silent = TRUE)

  ci <- confint(fit)
  expect_equal(nrow(ci), length(coef(fit)))
  expect_equal(ncol(ci), 2)
  # Lower < upper
  expect_true(all(ci[, 1] < ci[, 2]))
})#TEST_THAT

test_that("confint uniform produces wider bands for ddml_attgt", {
  set.seed(42)
  n <- 800; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  G <- sample(c(3, 4, Inf), n, replace = TRUE,
              prob = c(0.3, 0.3, 0.4))
  y <- matrix(rnorm(n * T_), n, T_)

  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                  learners = list(what = ols),
                  sample_folds = 2, cv_folds = 3,
                  silent = TRUE)

  ci_pw <- confint(fit)

  set.seed(1)
  ci_uf <- confint(fit, uniform = TRUE,
                    bootstraps = 499)

  # Uniform bands are wider (multiple cells -> p > 1)
  width_pw <- ci_pw[, 2] - ci_pw[, 1]
  width_uf <- ci_uf[, 2] - ci_uf[, 1]
  expect_true(all(width_uf >= width_pw))

  # Critical value attribute
  expect_false(is.null(attr(ci_uf, "crit_val")))
  expect_true(attr(ci_uf, "crit_val") > qnorm(0.975))

  # Pointwise has Gaussian quantile
  expect_equal(attr(ci_pw, "crit_val"), qnorm(0.975), tolerance = 1e-6)
})#TEST_THAT

test_that("ddml_attgt skips cells with empty control pool", {
  set.seed(42)
  n <- 2000; T_ <- 4
  X <- matrix(rnorm(n * 2), n, 2)
  # All units treated at period 2, 3, or 4 — no never-treated.
  # With notyettreated, the last period will have no controls.
  G <- sample(c(2, 3, 4), n, replace = TRUE,
              prob = c(0.3, 0.4, 0.3))
  y <- matrix(rnorm(n * T_), n, T_)

  # Should not error — empty-control cells are skipped
  fit <- ddml_attgt(y, X, t = 1:T_, G = G,
                   learners = list(what = ols),
                   control_group = "notyettreated",
                   sample_folds = 2, cv_folds = 3,
                   silent = TRUE)

  expect_s3_class(fit, "ddml_attgt")
  # All retained cells must have positive control counts
  expect_true(all(fit$cell_info$n_control > 0))
  expect_true(length(coef(fit)) > 0)
})#TEST_THAT
