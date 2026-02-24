test_that("standard S3 generic methods work correctly", {
  # Generate a minimal fitted object to test the plumbing
  set.seed(42)
  nobs <- 100
  X <- cbind(1, matrix(rnorm(nobs * 4), nobs, 4))
  D_tld <- X %*% runif(5) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(5) + rnorm(nobs)
  
  learners <- list(list(fun = ols), list(fun = ols))
  suppressWarnings({
    fit <- ddml_ate(y, D, X,
                    learners = learners,
                    ensemble_type = c("ols", "nnls"),
                    cv_folds = 2,
                    sample_folds = 2,
                    silent = TRUE)
  })
  
  # coef
  cf <- coef(fit)
  expect_true(is.numeric(cf))
  expect_length(cf, 2)
  
  # vcov
  V <- vcov(fit, ensemble_idx = 2)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(1, 1))
  
  # confint
  ci <- confint(fit, ensemble_idx = 1)
  expect_true(is.matrix(ci))
  expect_equal(dim(ci), c(1, 2))
  expect_identical(colnames(ci), c(" 2.5 %", "97.5 %"))
  
  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.ddml")
  expect_equal(s$nobs, nobs)
  
  # print.summary
  out <- capture_output(print(s))
  expect_true(grepl("DDML estimation", out))
  expect_true(grepl("Average Treatment Effect", out))
})
