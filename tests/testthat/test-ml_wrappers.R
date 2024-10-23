test_that("mdl_glmnet with cv is working", {
  # Simulate a small dataset
  nobs <- 100
  X <- matrix(rnorm(nobs*50), nobs, 50) # Simulate features
  y <- 1 + X %*% (10*runif(50) * (runif(50) < 0.1)) + rnorm(nobs)
  # Estimate the learner
  mdl_fit <- mdl_glmnet(y, X)
  # Check methods predict()
  fitted <- predict(mdl_fit, newdata = X)
  # Check output with expectations
  expect_equal(length(fitted), 100)
})#TEST_THAT

test_that("mdl_glmnet logit predicts probabilities", {
  # Simulate a small dataset
  nobs <- 100
  X <- matrix(rnorm(nobs*50), nobs, 50) # Simulate features
  y <- 1 * (X %*% (10*runif(50) * (runif(50) < 0.1)) + rnorm(nobs) > 0.5)
  # Estimate the learner
  mdl_fit <- mdl_glmnet(y, X, family = binomial)
  # Check methods predict()
  fitted <- predict(mdl_fit, newdata = X)
  # Check output with expectations
  expect_true(max(fitted) <= 1)
  expect_true(min(fitted) >= 0)
})#TEST_THAT

test_that("mdl_glmnet w/o cv is working", {
  # Simulate a small dataset
  nobs <- 100
  X <- matrix(rnorm(nobs*50), nobs, 50) # Simulate features
  y <- 1 + X %*% (10*runif(50) * (runif(50) < 0.1)) + rnorm(nobs)
  # Estimate the learner
  mdl_fit <- mdl_glmnet(y, X, cv = FALSE)
  # Check methods predict()
  fitted <- predict(mdl_fit, newdata = X)
  # Check output with expectations
  expect_equal(length(fitted), 100)
})#TEST_THAT

test_that("mdl_xgboost is working", {
  # Simulate a small dataset
  nobs <- 100
  X <- matrix(rnorm(nobs*50), nobs, 50) # Simulate features
  y <- 1 + X %*% (10*runif(50) * (runif(50) < 0.1))
  y <- 1 * (y - mean(y) >= rnorm(nobs))
  # Estimate the learner
  mdl_fit_reg <- mdl_xgboost(y, X)
  mdl_fit_probability <- mdl_xgboost(y, X, objective = "binary:logistic")
  # Check methods predict()
  fitted_ref <- predict(mdl_fit_reg, newdata = X)
  fitted_probability <- predict(mdl_fit_probability, newdata = X)
  # Check output with expectations
  expect_equal(length(fitted_ref), 100)
  expect_equal(length(fitted_probability), 100)
})#TEST_THAT

test_that("mdl_ranger is working", {
  # Simulate a small dataset
  nobs <- 100
  X <- matrix(rnorm(nobs*50), nobs, 50) # Simulate features
  y <- 1 + X %*% (10*runif(50) * (runif(50) < 0.1))
  y <- 1 * (y - mean(y) >= rnorm(nobs))
  # Estimate learners
  mdl_fit_reg <- mdl_ranger(y, X)
  mdl_fit_probability <- mdl_ranger(y, X, probability = TRUE)
  mdl_fit_classification <- mdl_ranger(y, X, classification = TRUE)
  # Check methods predict()
  fitted_reg <- predict(mdl_fit_reg, newdata = X)
  fitted_probability <- predict(mdl_fit_probability, newdata = X)
  expect_warning({
    fitted_classification <- predict(mdl_fit_classification, newdata = X)
  })
  # Check output with expectations
  expect_equal(length(fitted_reg), 100)
  expect_equal(length(fitted_probability), 100)
})#TEST_THAT
