# glmnet =======================================================================
#' Wrapper function for \code{glmnet}.
#'
#' Wrapper function for \code{glmnet}.
#'
#' @export mdl_glmnet
mdl_glmnet <- function(y, X,
                       alpha = 1, lambda = NULL,
                       standardize = TRUE, intercept = TRUE,
                       family = "gaussian",
                       cv = TRUE, nfolds = 10){
  # Either copute glmnet with given lambda or determine lambda with cv.
  if (cv) {
    mdl_fit <- glmnet::cv.glmnet(x = X, y = y,
                                 family = family,
                                 lambda = lambda,
                                 standardize = standardize,
                                 intercept = intercept,
                                 nfolds = nfolds)
  } else {
    mdl_fit <- glmnet::glmnet(x = X, y = y,
                              family = family,
                              lambda = lambda,
                              standardize = standardize,
                              intercept = intercept)
  }#IFELSE

  # Set custom S3 class
  class(mdl_fit) <- c("mdl_glmnet", class(mdl_fit))
  return(mdl_fit)
}#MDL_GLMNET

#' Predict method for mdl_glmnet fits.
#'
#' Predict method for mdl_glmnet fits.
#'
#' @export predict.mdl_glmnet
predict.mdl_glmnet <- function(obj, newdata = NULL){
  # Check whether cv.glmnet was run
  cv <- "cv.glmnet" %in% class(obj)
  # Compute predictions
  if (cv) {
    # Determine mse-minimizing lambda
    which_lambda <- which.min(obj$cvm)
    # Predict using glmnet prediction method
    fitted <- glmnet::predict.glmnet(obj$glmnet.fit, newdata,
                                     obj$lambda[which_lambda])
  } else {
    # Determine least regularizing lambda
    which_lambda <- length(obj$lambda)
    # Predict using glmnet prediction method
    fitted <- glmnet::predict.glmnet(obj, newdata,
                                     obj$lambda[which_lambda])
  }#IFELSE

  return(fitted)
}#PREDICT.MDL_GLMNET

#' Instrument selection for mdl_glmnet fits.
#'
#' Instrument selection for mdl_glmnet fits.
#'
#' @export any_iv.mdl_glmnet
any_iv.mdl_glmnet <- function(obj, index_iv, ...){
  # Check whether cv.glmnet was run
  cv <- "cv.glmnet" %in% class(obj)
  if (cv) {
    # Determine mse-minimizing lambda and get beta matrix
    which_lambda <- which.min(obj$cvm)
    beta <- obj$glmnet.fit$beta
  } else {
    # Determine least regularizing lambda and get beta matrix
    which_lambda <- length(obj$lambda)
    beta <- obj$beta
  }#IFELSE
  # Check whether any instruments are retained
  retained <- c(1, which((beta[, which_lambda] != 0)))
  length(intersect(retained, index_iv)) > 0
}#ANY_IV.MDL_GLMNET

# xgboost ======================================================================
#' Wrapper function for \code{xgboost}.
#'
#' Wrapper function for \code{xgboost}.
#'
#' @export mdl_xgboost
mdl_xgboost <- function(y, X,
                        num_parallel_tree = 1, min_child_weight = 1,
                        colsample_bytree = 0.7, subsample = 0.7,
                        colsample_bylevel = 1, colsample_bynode = 1,
                        max.depth = .Machine$integer.max, eta = 0.3, gamma = 0,
                        nrounds = 3, objective = "reg:squarederror",
                        verbose = 0){
  # Compute xgboost
  mdl_fit <- xgboost::xgboost(data = X, label = y,
                              num_parallel_tree = num_parallel_tree,
                              min_child_weight = min_child_weight,
                              colsample_bytree = colsample_bytree,
                              subsample = subsample,
                              colsample_bylevel = colsample_bylevel,
                              colsample_bynode = colsample_bynode,
                              max.depth = max.depth, eta = eta,
                              gamma = gamma, nrounds = nrounds,
                              objective = objective, verbose = verbose)
  # Set custom S3 class
  class(mdl_fit) <- c("mdl_xgboost", class(mdl_fit))
  return(mdl_fit)
}#MDL_XGBOOST

#' Predict method for mdl_xgboost fits.
#'
#' Predict method for mdl_xgboost fits.
#'
#' @export predict.mdl_xgboost
predict.mdl_xgboost <- function(obj, newdata = NULL){
  # Predict using xgb.Booster prediction method. Note that 'predict.xgb.Booster'
  #     is not an exported object from 'namespace:xgboost', hence the less ideal
  #     fix.
  class(obj) <- class(obj)[2] #
  predict(obj, newdata)
}#PREDICT.MDL_XGBOOST

#' Instrument selection for mdl_xgboost fits.
#'
#' Instrument selection for mdl_xgboost fits.
#'
#' @export any_iv.mdl_xgboost
any_iv.mdl_xgboost <- function(obj, index_iv, names_iv, ...){
  # Check whether names_iv is NULL, in which case a replacement is needed
  if (is.null(names_iv)) names_iv <- index_iv
  # Check whether instruments have non-zero varibale importance
  vselected <- xgboost::xgb.importance(model = obj)[, 1]
  ivselected <- intersect(sapply(vselected, function(x) x), names_iv)
  length(ivselected) > 0
}#ANY_IV.MDL_XGBOOST

# randomForest =================================================================
#' Wrapper function for \code{randomForest}.
#'
#' Wrapper function for \code{randomForest}.
#'
#' @export mdl_randomForest
mdl_randomForest <- function(y, X,
                             ntree = 100, nodesize = 1, maxnodes = NULL,
                             colsample_bytree = 0.6, subsample = 0.7,
                             replace = FALSE){
  # Compute randomForest
  if(!("matrix" %in% class(X))) X <- Matrix::as.matrix(X)
  mdl_fit <- randomForest::randomForest(X, y,
                                        ntree = ntree, nodesize = nodesize,
                                        maxnodes = maxnodes,
                                        mtry = ceiling(colsample_bytree *
                                                         ncol(X)),
                                        sampsize = ceiling(subsample *
                                                             length(y)),
                                        replace = replace)
  # Set custom S3 class
  class(mdl_fit) <- c("mdl_randomForest", class(mdl_fit))
  return(mdl_fit)
}#MDL_RANDOMFOREST

#' Predict method for mdl_randomForest fits.
#'
#' Predict method for mdl_randomForest fits.
#'
#' @export predict.mdl_randomForest
predict.mdl_randomForest <- function(obj, newdata = NULL){
  # Predict using xgb.Booster prediction method. Note that
  #     'predict.randomForest' is not an exported object from
  #     'namespace:randomForest', hence the less ideal fix.
  class(obj) <- class(obj)[2]
  # Predict using randomForest prediction method
  predict(obj, newdata)
}#PREDICT.MDL_RANDOMFOREST

#' Instrument selection for mdl_xgboost fits.
#'
#' Instrument selection for mdl_xgboost fits.
#'
#' @export any_iv.mdl_randomForest
any_iv.mdl_randomForest <- function(obj, index_iv, ...){
  # Check whether instruments have non-zero varibale importance
  vselected <- which(obj$importance > 0)
  length(intersect(vselected, index_iv)) > 0
}#ANY_IV.MDL_RANDOMFOREST

# grf ==========================================================================
#' Wrapper function for \code{grf}'s \code{regression_forest}.
#'
#' Wrapper function for \code{grf}'s \code{regression_forest}.
#'
#' @export mdl_grf
mdl_grf <- function(y, X,
                    num.trees = 100,
                    colsample_bytree = 0.6,
                    sample.fraction = 0.7,
                    min.node.size = 1,
                    honesty = TRUE){
  # Compute regression_forest
  if(!("matrix" %in% class(X))) X <- as.matrix(X)
  mdl_fit <- grf::regression_forest(X, y,
                                    num.trees = num.trees,
                                    mtry = ceiling(colsample_bytree * ncol(X)),
                                    sample.fraction = sample.fraction,
                                    min.node.size = min.node.size,
                                    honesty = honesty,
                                    ci.group.size = 1,
                                    compute.oob.predictions = F)
  # Organize and return output
  class(mdl_fit) <- c("mdl_grf", class(mdl_fit))
  return(mdl_fit)
}#MDL_GRF

#' Predict method for mdl_grf fits.
#'
#' Predict method for mdl_grf fits.
#'
#' @export predict.mdl_grf
predict.mdl_grf <- function(obj, newdata = NULL){
  # Check for new data
  #if(is.null(newdata)) newdata <- obj$X
  class(obj) <- class(obj)[2]
  # Predict data and output as matrix
  as.numeric(predict(obj, newdata)$predictions) # don't return a data.frame
}#PREDICT.MDL_GRF

#' Instrument selection for mdl_grf fits.
#'
#' Instrument selection for mdl_grf fits.
#'
#' @export any_iv.mdl_grf
any_iv.mdl_grf <- function(obj, index_iv, ...){
  # Check whether instruments have non-zero varibale importance
  vselected <- which(grf::variable_importance(obj) > 0)
  ivselected <- intersect(vselected, index_iv)
  length(ivselected) > 0
}#ANY_IV.MDL_GRF

# keras ========================================================================
#' Wrapper function for \code{keras}'s neural net implementation.
#'
#' Wrapper function for \code{keras}'s neural net implementation. Tensorflow is
#'     used as a backend for computation.
#'
#' @export mdl_keras
mdl_keras <- function(y, X,
                      model,
                      optimizer_fun = "rmsprop",
                      loss = "mse",
                      epochs = 10,
                      batch_size = min(1000, length(y)),
                      validation_split = 0,
                      callbacks = NULL,
                      steps_per_epoch = NULL,
                      metrics = c("mae"),
                      verbose = 0) {
  # Compile model
  model %>% keras::compile(optimizer = optimizer_fun,
                           loss = loss,
                           metrics = metrics)

  # Fit neural net
  model %>% keras::fit(X, y,
                       epochs = epochs,
                       batch_size = batch_size,
                       validation_split = validation_split,
                       callbacks = callbacks,
                       steps_per_epoch = steps_per_epoch,
                       verbose = verbose)
  # Return fit
  class(model) <- c("mdl_keras", class(model)) # amend class
  return(model)
}#MDL_KERAS

#' Predict method for mdl_keras fits.
#'
#' Predict method for mdl_keras fits.
#'
#' @export predict.mdl_keras
predict.mdl_keras <- function(obj, newdata = NULL){
  # Check for new data
  #if(is.null(newdata)) newdata <- obj$X
  # Predict data and output as matrix
  class(obj) <- class(obj)[-1] # Not a pretty solution...
  as.numeric(predict(obj, newdata))
}#PREDICT.MDL_KERAS

#' Instrument selection for mdl_keras fits.
#'
#' Instrument selection for mdl_keras fits.
#'
#' @export any_iv.mdl_keras
any_iv.mdl_keras <- function(obj, index_iv, ...){
  TRUE
}#ANY_IV.MDL_KERAS

# kcmeans + keras ==============================================================
#' Wrapper function for \code{keras}'s neural net implementation with
#'     initialization at the KCMeans solution.
#'
#' Wrapper function for \code{keras}'s neural net implementation with
#'     initialization at the KCMeans solution. Tensorflow is used as a backend
#'     for computation.
#'
#' @export mdl_keras_kcmeans
mdl_keras_kcmeans <- function(y, X,
                              model,
                              K,
                              optimizer_fun = "rmsprop",
                              loss = "mse",
                              epochs = 10,
                              batch_size = min(1000, length(y)),
                              validation_split = 0,
                              callbacks = NULL,
                              steps_per_epoch = NULL,
                              metrics = c("mae"),
                              verbose = 0) {
  # Compile model
  model %>% keras::compile(optimizer = optimizer_fun,
                           loss = loss,
                           metrics = metrics)

  # Initialize with the KCMeans solution
  model <- initialize_with_kcmeans(model, y,  X[, -1], X[, 1], K)

  # Fit neural net
  model %>% keras::fit(X[, -1], y,
                       epochs = epochs,
                       batch_size = batch_size,
                       validation_split = validation_split,
                       callbacks = callbacks,
                       steps_per_epoch = steps_per_epoch,
                       verbose = verbose)
  # Return fit
  class(model) <- c("mdl_keras_kcmeans", class(model)) # amend class
  return(model)
}#MDL_KERAS

#' Predict method for mdl_keras_kcmeans fits.
#'
#' Predict method for mdl_keras_kcmeans fits.
#'
#' @export predict.mdl_keras_kcmeans
predict.mdl_keras_kcmeans <- function(obj, newdata = NULL){
  # Check for new data
  #if(is.null(newdata)) newdata <- obj$X
  # Predict data and output as matrix
  class(obj) <- class(obj)[-1] # Not a pretty solution...
  as.numeric(predict(obj, newdata[, -1]))
}#PREDICT.any_iv.mdl_keras_kcmeans

#' Instrument selection for mdl_keras fits.
#'
#' Instrument selection for mdl_keras fits.
#'
#' @export any_iv.mdl_keras_kcmeans
any_iv.mdl_keras_kcmeans <- function(obj, index_iv, ...){
  TRUE
}#ANY_IV.any_iv.mdl_keras_kcmeans

#' Internal function for KCMeans initialization.
#'
initialize_with_kcmeans <- function(model, y, X,
                                    x_kcmeans, K) {
  # Estimate a single epoch to initialize weights
  model %>% keras::fit(X, y,
                       epochs = 1,
                       batch_size = 100,
                       verbose = F)

  # Calculate KCMeans solution
  kcmeans_fit <- kcmeans(y, x_kcmeans, K)
  # Set weights to KCMeans solution
  w_kcmeans <- get_weights(model)
  nhidden <- length(w_kcmeans) - 2
  for (layer in 1:nhidden) {
    if (layer == 1) {
      for (v in 1:K) {
        indx_v <- kcmeans_fit$cluster_map[[v]]
        w_kcmeans[[layer]][, v] <- 0 # runif(length(w_kcmeans[[1]][, v]), -1, 1)
        w_kcmeans[[layer]][indx_v, v] <- 1
      }#FOR
    } else {
      w_kcmeans[[layer]] <- diag(1, dim(w_kcmeans[[1]])[2])
    }#IF
  }#FOR
  w_kcmeans[[nhidden + 1]][, 1] <- kcmeans_fit$alpha - mean(y)
  w_kcmeans[[nhidden + 2]] <- array(mean(y), dim = c(1))
  set_weights(model, w_kcmeans)

  # Return initialized model
  return(model)
}#INITIALIZE_WITH_KCMEANS

