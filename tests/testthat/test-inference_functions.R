test_that("compute_ddml_variance computes HC1 sandwich variance correctly", {
  set.seed(42)
  nobs <- 100
  p <- 2
  scores <- matrix(rnorm(nobs * p), nobs, p)
  J <- diag(2, p)

  meat <- crossprod(scores) / nobs
  J_inv <- solve(J)
  expected_v <- J_inv %*% meat %*% t(J_inv) * nobs / (nobs - p) / nobs

  V <- compute_ddml_variance(scores, J)

  expect_true(is.matrix(V))
  expect_equal(dim(V), c(p, p))
  expect_equal(V, expected_v)
})

test_that("compute_ddml_variance handles cluster aggregation correctly", {
  set.seed(42)
  nobs <- 100
  p <- 2
  scores <- matrix(rnorm(nobs * p), nobs, p)
  J <- diag(2, p)
  cluster_var <- rep(1:50, each = 2)

  V_clustered <- compute_ddml_variance(
    scores, J, cluster_variable = cluster_var)

  sc_sum <- rowsum(scores, cluster_var)
  n_eff <- 50
  meat <- crossprod(sc_sum) / n_eff
  J_inv <- solve(J)
  expected_v <- J_inv %*% meat %*% t(J_inv) *
    n_eff / (n_eff - p) / n_eff

  expect_equal(V_clustered, expected_v)
})

test_that("compute_ddml_variance computes HC0 correctly", {
  set.seed(42)
  nobs <- 100
  p <- 2
  scores <- matrix(rnorm(nobs * p), nobs, p)
  J <- diag(2, p)

  meat <- crossprod(scores) / nobs
  J_inv <- solve(J)
  expected_v <- J_inv %*% meat %*% t(J_inv) / nobs

  V <- compute_ddml_variance(scores, J, type = "HC0")

  expect_equal(V, expected_v)
})

test_that("HC1 equals HC0 scaled by n/(n-p)", {
  set.seed(42)
  nobs <- 100
  p <- 2
  scores <- matrix(rnorm(nobs * p), nobs, p)
  J <- diag(2, p)

  V_hc0 <- compute_ddml_variance(scores, J, type = "HC0")
  V_hc1 <- compute_ddml_variance(scores, J, type = "HC1")

  expect_equal(V_hc1, V_hc0 * nobs / (nobs - p))
})

test_that("compute_ddml_variance computes HC3 correctly with leverage", {
  set.seed(42)
  nobs <- 100
  p <- 2
  scores <- matrix(rnorm(nobs * p), nobs, p)
  J <- diag(2, p)
  h <- rep(p / nobs, nobs)

  sc_adj <- scores / (1 - h)
  meat <- crossprod(sc_adj) / nobs
  J_inv <- solve(J)
  expected_v <- J_inv %*% meat %*% t(J_inv) / nobs

  V <- compute_ddml_variance(scores, J, type = "HC3",
                             leverage = h)

  expect_true(is.matrix(V))
  expect_equal(dim(V), c(p, p))
  expect_equal(V, expected_v)
})

test_that("HC3 with hat-matrix leverage matches manual computation", {
  set.seed(42)
  nobs <- 200
  p <- 3
  X <- matrix(rnorm(nobs * p), nobs, p)
  h <- rowSums((X %*% solve(crossprod(X))) * X)

  expect_equal(sum(h), p, tolerance = 1e-10)
  expect_true(all(h >= 0))
  expect_true(all(h < 1))
  expect_equal(mean(h), p / nobs, tolerance = 1e-10)
})

test_that("HC3 with cluster aggregation and leverage", {
  set.seed(42)
  nobs <- 100
  p <- 2
  scores <- matrix(rnorm(nobs * p), nobs, p)
  J <- diag(2, p)
  cluster_var <- rep(1:50, each = 2)
  h_obs <- rep(p / nobs, nobs)

  V <- compute_ddml_variance(
    scores, J,
    cluster_variable = cluster_var,
    type = "HC3",
    leverage = h_obs)

  # Manual: aggregate scores and leverage to cluster level
  sc_sum <- rowsum(scores, cluster_var)
  h_clust <- as.vector(tapply(h_obs, cluster_var, sum))
  n_eff <- 50
  sc_adj <- sc_sum / (1 - h_clust)
  meat <- crossprod(sc_adj) / n_eff
  J_inv <- solve(J)
  expected_v <- J_inv %*% meat %*% t(J_inv) / n_eff

  expect_equal(V, expected_v)
})

test_that("type rejects invalid values", {
  scores <- matrix(rnorm(20), 10, 2)
  J <- diag(2)

  expect_error(
    compute_ddml_variance(scores, J, type = "HC2"),
    "arg")
})

test_that("HC3 >= HC1 >= HC0 for typical data", {
  set.seed(42)
  nobs <- 200
  p <- 1
  scores <- matrix(rnorm(nobs * p), nobs, p)
  J <- matrix(-1, 1, 1)
  h <- rep(p / nobs, nobs)

  se_hc0 <- sqrt(diag(compute_ddml_variance(
    scores, J, type = "HC0")))
  se_hc1 <- sqrt(diag(compute_ddml_variance(
    scores, J, type = "HC1")))
  se_hc3 <- sqrt(diag(compute_ddml_variance(
    scores, J, type = "HC3", leverage = h)))

  expect_true(all(se_hc1 >= se_hc0))
  expect_true(all(se_hc3 >= se_hc0))
})
