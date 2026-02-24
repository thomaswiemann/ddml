test_that("compute_ddml_variance computes HC1 sandwich variance correctly", {
  # Mock data
  set.seed(42)
  nobs <- 100
  p <- 2
  scores <- matrix(rnorm(nobs * p), nobs, p)
  J <- diag(2, p) # Simple J matrix
  
  # Manual computation matching HC1 specification
  meat <- crossprod(scores) / nobs
  J_inv <- solve(J)
  expected_v <- J_inv %*% meat %*% t(J_inv) * nobs / (nobs - p) / nobs
  
  # Test function output
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
  cluster_var <- rep(1:50, each = 2) # 50 clusters
  
  # By clustering, the scores should be rowsummed to n_eff = 50
  V_clustered <- compute_ddml_variance(scores, J, cluster_variable = cluster_var)
  
  sc_sum <- rowsum(scores, cluster_var)
  n_eff <- 50
  meat <- crossprod(sc_sum) / n_eff
  J_inv <- solve(J)
  expected_v <- J_inv %*% meat %*% t(J_inv) * n_eff / (n_eff - p) / n_eff
  
  expect_equal(V_clustered, expected_v)
})
