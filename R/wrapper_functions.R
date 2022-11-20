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
                       cv = TRUE, nfolds = 10,
                       thres = 1e-7, ...){
  # Either copute glmnet with given lambda or determine lambda with cv.
  if (cv) {
    mdl_fit <- glmnet::cv.glmnet(x = X, y = y,
                                 family = family,
                                 lambda = lambda,
                                 standardize = standardize,
                                 intercept = intercept,
                                 nfolds = nfolds,
                                 thres = thres, ...)
  } else {
    mdl_fit <- glmnet::glmnet(x = X, y = y,
                              family = family,
                              lambda = lambda,
                              standardize = standardize,
                              intercept = intercept,
                              thres = thres, ...)
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
#' @export
predict.mdl_glmnet <- function(object, newdata = NULL, ...){
  # Check whether cv.glmnet was run
  cv <- "cv.glmnet" %in% class(object)
  # Compute predictions
  if (cv) {
    # Determine mse-minimizing lambda
    which_lambda <- which.min(object$cvm)
    # Predict using glmnet prediction method
    fitted <- glmnet::predict.glmnet(object$glmnet.fit, newdata,
                                     object$lambda[which_lambda], ...)
  } else {
    # Determine least regularizing lambda
    which_lambda <- length(object$lambda)
    # Predict using glmnet prediction method
    fitted <- glmnet::predict.glmnet(object, newdata,
                                     object$lambda[which_lambda], ...)
  }#IFELSE

  return(fitted)
}#PREDICT.MDL_GLMNET

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
                        nrounds = 500, objective = "reg:squarederror",
                        interaction_constraints = list(),
                        verbose = 0, ...){
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
                              objective = objective,
                              interaction_constraints = interaction_constraints,
                              verbose = verbose, ...)
  # Set custom S3 class
  class(mdl_fit) <- c("mdl_xgboost", class(mdl_fit))
  return(mdl_fit)
}#MDL_XGBOOST

#' Predict method for mdl_xgboost fits.
#'
#' Predict method for mdl_xgboost fits.
#'
#' @export predict.mdl_xgboost
#' @export
predict.mdl_xgboost <- function(object, newdata = NULL, ...){
  # Predict using xgb.Booster prediction method. Note that 'predict.xgb.Booster'
  #     is not an exported object from 'namespace:xgboost', hence the less ideal
  #     fix.
  class(object) <- class(object)[2] #
  stats::predict(object, newdata, ...)
}#PREDICT.MDL_XGBOOST

# randomForest =================================================================
#' Wrapper function for \code{randomForest}.
#'
#' Wrapper function for \code{randomForest}.
#'
#' @export mdl_randomForest
mdl_randomForest <- function(y, X,
                             ntree = 100, nodesize = 1, maxnodes = NULL,
                             colsample_bytree = 0.6, subsample = 0.7,
                             replace = FALSE, ...){
  # Compute randomForest
  if(!("matrix" %in% class(X))) X <- Matrix::as.matrix(X)
  mdl_fit <- randomForest::randomForest(X, y,
                                        ntree = ntree, nodesize = nodesize,
                                        maxnodes = maxnodes,
                                        mtry = ceiling(colsample_bytree *
                                                         ncol(X)),
                                        sampsize = ceiling(subsample *
                                                             length(y)),
                                        replace = replace, ...)
  # Set custom S3 class
  class(mdl_fit) <- c("mdl_randomForest", class(mdl_fit))
  return(mdl_fit)
}#MDL_RANDOMFOREST

#' Predict method for mdl_randomForest fits.
#'
#' Predict method for mdl_randomForest fits.
#'
#' @export predict.mdl_randomForest
#' @export
predict.mdl_randomForest <- function(object, newdata = NULL, ...){
  # Predict using xgb.Booster prediction method. Note that
  #     'predict.randomForest' is not an exported object from
  #     'namespace:randomForest', hence the less ideal fix.
  class(object) <- class(object)[2]
  # Predict using randomForest prediction method
  stats::predict(object, newdata, ...)
}#PREDICT.MDL_RANDOMFOREST

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
                    honesty = TRUE, ...){
  # Compute regression_forest
  if(!("matrix" %in% class(X))) X <- as.matrix(X)
  mdl_fit <- grf::regression_forest(X, y,
                                    num.trees = num.trees,
                                    mtry = ceiling(colsample_bytree * ncol(X)),
                                    sample.fraction = sample.fraction,
                                    min.node.size = min.node.size,
                                    honesty = honesty,
                                    ci.group.size = 1,
                                    compute.oob.predictions = F, ...)
  # Organize and return output
  class(mdl_fit) <- c("mdl_grf", class(mdl_fit))
  return(mdl_fit)
}#MDL_GRF

#' Predict method for mdl_grf fits.
#'
#' Predict method for mdl_grf fits.
#'
#' @export predict.mdl_grf
#' @export
predict.mdl_grf <- function(object, newdata = NULL, ...){
  # Check for new data
  #if(is.null(newdata)) newdata <- object$X
  class(object) <- class(object)[2]
  # Predict data and output as matrix
  as.numeric(stats::predict(object, newdata, ...)$predictions)
}#PREDICT.MDL_GRF
