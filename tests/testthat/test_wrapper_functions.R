library(ddml)
context("Testing wrapper functions for imported ML procedures.")

sim_data <- function(N = 100) {
  X <- matrix(rnorm(N*50), N, 50) # Simulate features
  y <- 1 + X %*% (10*runif(50) * (runif(50) < 0.1)) + rnorm(N)
  return(list(y = y, X = X))
}#SIM_DATA

test_that("mdl_glmnet with cv is working", {
  # Simulate small dataset and fit the model
  dat <- sim_data()
  mdl_fit <- mdl_glmnet(dat$y, dat$X)
  # Check methods predict() and any_iv()
  fitted <- predict(mdl_fit, newdata = dat$X)
  iv_selected <- any_iv(mdl_fit, index_iv = c(11:20))
  # Check output with expectations
  expect_is(mdl_fit, "mdl_glmnet")
  expect_equal(length(fitted), 100)
  expect_is(iv_selected, "logical")
})#TEST_THAT

test_that("mdl_glmnet w/o cv is working", {
  # Simulate small dataset and fit the model
  dat <- sim_data()
  mdl_fit <- mdl_glmnet(dat$y, dat$X, cv = FALSE)
  # Check methods predict() and any_iv()
  fitted <- predict(mdl_fit, newdata = dat$X)
  iv_selected <- any_iv(mdl_fit, index_iv = c(11:20))
  # Check output with expectations
  expect_is(mdl_fit, "mdl_glmnet")
  expect_equal(length(fitted), 100)
  expect_is(iv_selected, "logical")
})#TEST_THAT

test_that("mdl_xgboost is working", {
  # Simulate small dataset and fit the model
  dat <- sim_data()
  mdl_fit <- mdl_xgboost(dat$y, dat$X)
  # Check methods predict() and any_iv()
  fitted <- predict(mdl_fit, newdata = dat$X)
  iv_selected <- any_iv(mdl_fit, index_iv = c(11:20),
                        names_iv = colnames(dat$X))
  # Check output with expectations
  expect_is(mdl_fit, "mdl_xgboost")
  expect_equal(length(fitted), 100)
  expect_is(iv_selected, "logical")
})#TEST_THAT

test_that("mdl_randomForest is working", {
  # Simulate small dataset and fit the model
  dat <- sim_data()
  mdl_fit <- mdl_randomForest(dat$y, dat$X)
  # Check methods predict() and any_iv()
  fitted <- predict(mdl_fit, newdata = dat$X)
  iv_selected <- any_iv(mdl_fit, index_iv = c(11:20))
  # Check output with expectations
  expect_is(mdl_fit, "mdl_randomForest")
  expect_equal(length(fitted), 100)
  expect_is(iv_selected, "logical")
})#TEST_THAT

test_that("mdl_grf is working", {
  # Simulate small dataset and fit the model
  dat <- sim_data()
  mdl_fit <- mdl_grf(dat$y, dat$X)
    # Check methods predict() and any_iv()
  fitted <- predict(mdl_fit, newdata = dat$X)
  iv_selected <- any_iv(mdl_fit, index_iv = c(11:20))
  # Check output with expectations
  expect_is(mdl_fit, "mdl_grf")
  expect_equal(length(fitted), 100)
  expect_is(iv_selected, "logical")
})#TEST_THAT

# Skip this for now -- takes to long
if(FALSE) {
test_that("mdl_keras is working", {
  # Simulate small dataset
  dat <- sim_data(N = 1000)
  # Build simple neural net
  model <- keras::keras_model_sequential() %>%
    keras::layer_dense(units = 10, activation = "relu",
                       input_shape = dim(dat$X)[[2]]) %>%
    keras::layer_dense(units = 10, activation = "relu") %>%
    keras::layer_dense(units = 1)
  # Estimate model
  mdl_fit <- mdl_keras(dat$y, dat$X, model,
                       epochs = 10)
  # Check methods predict() and any_iv()
  fitted <- predict(mdl_fit, newdata = dat$X)
  iv_selected <- any_iv(mdl_fit, index_iv = c(11:20))
  # Check output with expectations
  expect_is(mdl_fit, "mdl_keras")
  expect_equal(length(fitted), 1000)
  expect_is(iv_selected, "logical")
})#TEST_THAT
}
