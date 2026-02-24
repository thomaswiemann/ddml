test_that("cvc_pairwise detects a dominant learner", {
  set.seed(42)
  n <- 500
  # Learner 1 has small residuals, learner 2 has large residuals
  resid1 <- rnorm(n, sd = 0.5)
  resid2 <- rnorm(n, sd = 2.0)
  fid <- rep(seq_len(5), each = n / 5)

  # zeta = resid1^2 - resid2^2 < 0 when learner 1 is better
  # → Tx < 0 → p-value near 1 (cannot reject that 1 is worse)
  pval_12 <- cvc_pairwise(resid1, resid2, fid, bootnum = 500)
  expect_true(pval_12 > 0.9)

  # Reversed: learner 2 worse → p-value near 0
  pval_21 <- cvc_pairwise(resid2, resid1, fid, bootnum = 500)
  expect_true(pval_21 < 0.1)
})

test_that("cvc_pairwise moderate for equal learners", {
  set.seed(123)
  n <- 1000
  resid1 <- rnorm(n)
  resid2 <- rnorm(n)
  fid <- rep(seq_len(5), each = n / 5)

  pval <- cvc_pairwise(resid1, resid2, fid, bootnum = 500)

  expect_true(pval > 0.05 && pval < 0.95)
})

test_that("cvc_pairwise returns NA for identical residuals", {
  n <- 100
  resid <- rnorm(n)
  fid <- rep(seq_len(5), each = n / 5)

  pval <- cvc_pairwise(resid, resid, fid, bootnum = 100)
  expect_true(is.na(pval))
})

test_that("cvc_confidence_set includes the best learner", {
  set.seed(42)
  n <- 500
  resid_mat <- cbind(
    rnorm(n, sd = 0.5),
    rnorm(n, sd = 1.5),
    rnorm(n, sd = 2.0))
  fid <- rep(seq_len(5), each = n / 5)

  cs <- cvc_confidence_set(resid_mat, fid,
                           bootnum = 500, alpha = 0.10)

  expect_length(cs$pvalues, 3)
  expect_length(cs$confidence_set, 3)
  # Best learner (column 1) should be in the confidence set
  expect_true(cs$confidence_set[1])
})

test_that("cvc_test returns correct structure", {
  set.seed(42)
  n <- 200
  resid_mat <- cbind(rnorm(n), rnorm(n), rnorm(n))
  colnames(resid_mat) <- c("ols", "lasso", "ridge")
  subsamples <- list(1:100, 101:200)

  res <- cvc_test(resid_mat, subsamples,
                  bootnum = 100, alpha = 0.05)

  expect_true(is.list(res))
  expect_equal(dim(res$pairwise), c(3, 3))
  expect_true(all(is.na(diag(res$pairwise))))
  expect_length(res$pvalues, 3)
  expect_length(res$confidence_set, 3)
  expect_identical(rownames(res$pairwise),
                   c("ols", "lasso", "ridge"))
})

test_that("cvc_test handles unequal fold sizes", {
  set.seed(42)
  n <- 103
  resid_mat <- cbind(rnorm(n), rnorm(n))
  subsamples <- list(1:51, 52:103)

  res <- cvc_test(resid_mat, subsamples, bootnum = 100)

  expect_equal(dim(res$pairwise), c(2, 2))
  expect_length(res$pvalues, 2)
})
