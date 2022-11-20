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
  y <- 1 + X %*% (10*runif(50) * (runif(50) < 0.1)) + rnorm(nobs)
  # Estimate the learner
  mdl_fit <- mdl_xgboost(y, X)
  # Check methods predict()
  fitted <- predict(mdl_fit, newdata = X)
  # Check output with expectations
  expect_equal(length(fitted), 100)
})#TEST_THAT

test_that("mdl_randomForest is working", {
  # Simulate a small dataset
  nobs <- 100
  X <- matrix(rnorm(nobs*50), nobs, 50) # Simulate features
  y <- 1 + X %*% (10*runif(50) * (runif(50) < 0.1)) + rnorm(nobs)
  # Estimate the learner
  mdl_fit <- mdl_randomForest(as.vector(y), X)
  # Check methods predict()
  fitted <- predict(mdl_fit, newdata = X)
  # Check output with expectations
  expect_equal(length(fitted), 100)
})#TEST_THAT

test_that("mdl_grf is working", {
  # Simulate a small dataset
  nobs <- 100
  X <- matrix(rnorm(nobs*50), nobs, 50) # Simulate features
  y <- 1 + X %*% (10*runif(50) * (runif(50) < 0.1)) + rnorm(nobs)
  # Estimate the learner
  mdl_fit <- mdl_grf(y, X)
    # Check methods predict()
  fitted <- predict(mdl_fit, newdata = X)
  # Check output with expectations
  expect_equal(length(fitted), 100)
})#TEST_THAT

# Skip this for now -- takes to long
if(FALSE) {
test_that("mdl_keras is working", {
  # Simulate a small dataset
  nobs <- 100
  X <- matrix(rnorm(nobs*50), nobs, 50) # Simulate features
  y <- 1 + X %*% (10*runif(50) * (runif(50) < 0.1)) + rnorm(nobs)
  # Build simple neural net
  model <- keras::keras_model_sequential() %>%
    keras::layer_dense(units = 10, activation = "relu",
                       input_shape = dim(X)[[2]]) %>%
    keras::layer_dense(units = 10, activation = "relu") %>%
    keras::layer_dense(units = 1)
  # Estimate model
  mdl_fit <- mdl_keras(y, X, model,
                       epochs = 10)
  # Check methods predict()
  fitted <- predict(mdl_fit, newdata = X)
  # Check output with expectations
  expect_equal(length(fitted), 1000)
})#TEST_THAT
}
