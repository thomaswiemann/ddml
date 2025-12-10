# glmnet =======================================================================

#' Wrapper for [glmnet::glmnet()].
#'
#' @family ml_wrapper
#'
#' @seealso [glmnet::glmnet()],[glmnet::cv.glmnet()]
#'
#' @description Simple wrapper for [glmnet::glmnet()] and [glmnet::cv.glmnet()].
#'
#' @param y The outcome variable.
#' @param X The (sparse) feature matrix.
#' @param cv Boolean to indicate use of lasso with cross-validated penalty.
#' @param ... Additional arguments passed to \code{glmnet}. See
#'     [glmnet::glmnet()] and [glmnet::cv.glmnet()] for a complete list of
#'     arguments.
#'
#' @return \code{mdl_glmnet} returns an object of S3 class \code{mdl_glmnet} as
#'     a simple mask of the return object of [glmnet::glmnet()] or
#'     [glmnet::cv.glmnet()].
#' @export
#'
#' @references
#' Friedman J, Hastie T, Tibshirani R (2010). "Regularization Paths for
#'     Generalized Linear Models via Coordinate Descent." Journal of Statistical
#'     Software, 33(1), 1–22.
#'
#' Simon N, Friedman J, Hastie T, Tibshirani R (2011). "Regularization Paths for
#'     Cox's Proportional Hazards Model via Coordinate Descent." Journal of
#'     Statistical Software, 39(5), 1–13.
#'
#' @examples
#' glmnet_fit <- mdl_glmnet(rnorm(100), matrix(rnorm(1000), 100, 10))
#' class(glmnet_fit)
mdl_glmnet <- function(y, X,
                       cv = TRUE,
                       ...){
  # Either copute glmnet with given lambda or determine lambda with cv.
  if (cv) {
    mdl_fit <- glmnet::cv.glmnet(x = X, y = y, ...)
  } else {
    mdl_fit <- glmnet::glmnet(x = X, y = y, ...)
  }#IFELSE

  # Set custom S3 class
  class(mdl_fit) <- c("mdl_glmnet", class(mdl_fit))
  return(mdl_fit)
}#MDL_GLMNET

#' @exportS3Method
predict.mdl_glmnet <- function(object, newdata = NULL, ...){
  # Check whether cv.glmnet was run
  cv <- "cv.glmnet" %in% class(object)
  class(object) <- class(object)[-1]
  # Compute predictions
  if (cv) {
    # Determine mse-minimizing lambda
    which_lambda <- which.min(object$cvm)
    # Predict using glmnet prediction method
    fitted <- stats::predict(object$glmnet.fit, newx = newdata,
                             s = object$lambda[which_lambda],
                             type = "response", ...)
  } else {
    # Determine least regularizing lambda
    which_lambda <- length(object$lambda)
    # Predict using glmnet prediction method
    fitted <- stats::predict(object, newx = newdata,
                             s = object$lambda[which_lambda],
                             type = "response", ...)
  }#IFELSE

  return(fitted)
}#PREDICT.MDL_GLMNET

# xgboost ======================================================================

#' Wrapper for [xgboost::xgboost()].
#'
#' @family ml_wrapper
#'
#' @seealso [xgboost::xgboost()]
#'
#' @description Simple wrapper for [xgboost::xgboost()] with some changes to the
#'     default arguments.
#'
#' @inheritParams xgboost::xgboost
#' @param y The outcome variable.
#' @param X The (sparse) feature matrix.
#' @param ... Additional arguments passed to \code{xgboost}. See
#'     [xgboost::xgboost()] for a complete list of arguments.
#'
#' @return \code{mdl_xgboost} returns an object of S3 class \code{mdl_xgboost}
#'     as a simple mask to the return object of [xgboost::xgboost()].
#' @export
#'
#' @references
#' Chen T, Guestrin C (2011). "Xgboost: A Scalable Tree Boosting System."
#' Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge
#' Discovery and Data Mining, 785–794.
#'
#' @examples
#' xgboost_fit <- mdl_xgboost(rnorm(50), matrix(rnorm(150), 50, 3),
#'                            nrounds = 1)
#' class(xgboost_fit)
mdl_xgboost <- function(y, X,
                        nrounds = 500, verbosity = 0,
                        ...){
  # Compute xgboost
  mdl_fit <- xgboost::xgboost(x = X, y = y,
                              nrounds = nrounds,
                              verbosity = verbosity, ...)
  # Set custom S3 class
  class(mdl_fit) <- c("mdl_xgboost", class(mdl_fit))
  return(mdl_fit)
}#MDL_XGBOOST

#' @exportS3Method
predict.mdl_xgboost <- function(object, newdata = NULL, ...){
  # Predict using xgb.Booster prediction method.
  class(object) <- class(object)[-1]
  stats::predict(object, newdata, ...)
}#PREDICT.MDL_XGBOOST

# ranger =======================================================================

#' Wrapper for [ranger::ranger()].
#'
#' @family ml_wrapper
#'
#' @seealso [ranger::ranger()]
#'
#' @description Simple wrapper for [ranger::ranger()]. Supports regression
#'     (default) and probability forests (set \code{probability = TRUE}).
#'
#' @param y The outcome variable.
#' @param X The feature matrix.
#' @param ... Additional arguments passed to \code{ranger}. See
#'     [ranger::ranger()] for a complete list of arguments.
#'
#' @return \code{mdl_ranger} returns an object of S3 class \code{ranger} as a
#'     simple mask of the return object of [ranger::ranger()].
#' @export
#'
#' @references
#' Wright M N, Ziegler A (2017). "ranger: A fast implementation of random
#'     forests for high dimensional data in C++ and R." Journal of Statistical
#'     Software 77(1), 1-17.
#'
#' @examples
#' ranger_fit <- mdl_ranger(rnorm(100), matrix(rnorm(1000), 100, 10))
#' class(ranger_fit)
mdl_ranger <- function(y, X, ...){
  # Assign columnames to X if none are given
  if (is.null(colnames(X))) {
    colnames(X) <- seq(dim(X)[2])
  }#IF
  # Compute ranger
  mdl_fit <- ranger::ranger(y = y, x = X, ...)
  # Set custom S3 class
  class(mdl_fit) <- c("mdl_ranger", class(mdl_fit))
  return(mdl_fit)
}#MDL_RANGER

#' @exportS3Method
predict.mdl_ranger <- function(object, newdata = NULL, ...){
  # Assign column names to newdata if none are given
  if (is.null(colnames(newdata))) {
    colnames(newdata) <- seq(dim(newdata)[2])
  }#IF
  class(object) <- class(object)[2]
  # Predict using randomForest prediction method
  if (object$treetype == "Probability estimation") {
    #stats::predict(object, data = newdata, ...)$predictions[, 2]
    stats::predict(object, data = newdata, ...)$predictions[, 2]
  } else if (object$treetype == "Regression") {
    #stats::predict(object, data = newdata, ...)$predictions
    stats::predict(object, data = newdata, ...)$predictions
  } else {
    warning("mdl_ranger is only designed for regression and probability forests")
    stats::predict(object, data = newdata, ...)$predictions
  }#IFELSE
}#PREDICT.MDL_RANGER

# glm ==========================================================================

#' Wrapper for [stats::glm()].
#'
#' @family ml_wrapper
#'
#' @seealso [stats::glm()]
#'
#' @description Simple wrapper for [stats::glm()].
#'
#' @param y The outcome variable.
#' @param X The feature matrix.
#' @param ... Additional arguments passed to \code{glm}. See
#'     [stats::glm()] for a complete list of arguments.
#'
#' @return \code{mdl_glm} returns an object of S3 class \code{mdl_glm} as a
#'     simple mask of the return object of [stats::glm()].
#' @export
#'
#' @examples
#' glm_fit <- mdl_glm(sample(0:1, 100, replace = TRUE),
#'                    matrix(rnorm(1000), 100, 10))
#' class(glm_fit)
mdl_glm <- function(y, X, ...) {
  df <- data.frame(y, X) # transform data from matrices to data.frame
  glm_fit <- stats::glm(y ~ ., data = df, ...) # fit glm
  class(glm_fit) <- c("mdl_glm", class(glm_fit)) # append class
  return(glm_fit) # return fitted glm object
}#MDL_GLM

#' @exportS3Method
predict.mdl_glm <- function(object, newdata, ...) {
  df <- data.frame(newdata) # transform data from matrices to data.frame
  stats::predict.glm(object, df, type = "response", ...)
}#PREDICT.MDL_GLM
