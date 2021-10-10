library(ddml)
context("Testing postLassoTSLS.")

sim_dat <- function(nobs) {
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  Z <- matrix(rnorm(nobs*10), nobs, 10) # overidentified
  UV <- matrix(rnorm(2*nobs), nobs, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2, 2))
  D <-  X %*% runif(40) + Z %*% c(1, runif(9)) + UV[, 1]
  y <- D + X %*% runif(40) + UV[, 2]
  # Organize and return output
  output <- list(y = y, D = D, Z = Z, X = X)
  return(output)
}#SIM_DAT

test_that("postLassoTSLS computes", {
  # Simulate small dataset
  dat <- sim_dat(100)
  # Define arguments
  pLTSLS_fit <- postLassoTSLS(dat$y, dat$X, dat$D, dat$Z,
                              splitSample = FALSE,
                              splitSample.IVChoice = FALSE,
                              penaltyMethod = 'r',
                              K.resid = 1, d = 5,
                              heteroskedasticity = FALSE)
  # Check output with expectations
  expect_equal(length(pLTSLS_fit$coef), 1)
})#TEST_THAT
