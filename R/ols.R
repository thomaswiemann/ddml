#' Title
#'
#' @param y abc
#' @param X abc
#' @param const abc
#' @param w abc
#'
#' @return output
#' @export
#'
#' @examples
#' ols(rnorm(100), cbind(rnorm(100), rnorm(100)))
ols <- function(y, X,
                const = FALSE,
                w = NULL) {
  # Add constant (optional)
  if (const) X <- cbind(1, X)

  # Data parameters
  calc_wls <- !is.null(w)

  # Compute OLS coefficient
  if (!calc_wls) {
    XX_inv <- csolve(as.matrix(Matrix::crossprod(X)))
    coef <- XX_inv %*% Matrix::crossprod(X, y)
  } else { # Or calculate WLS coefficient whenever weights are specified
    Xw <- X * w # weight rows
    XX_inv <- csolve(as.matrix(Matrix::crossprod(Xw, X)))
    coef <- XX_inv %*% Matrix::crossprod(Xw, y)
  }#IFELSE
  # Return estimate
  coef <- as.matrix(coef)
  try(rownames(coef) <- colnames(X)) # assign coefficient names
  output <- list(coef = coef, y = y, X = X,
                 const = const, w = w)
  class(output) <- "ols" # define S3 class
  return(output)
}#OLS

# Complementary methods ========================================================


#' Title
#'
#' @param object abc
#' @param newdata abc
#' @param ... abc
#'
#' @return output
#' @export predict.ols
#' @export
#'
#' @examples
#' ols_fit <- ols(rnorm(100), cbind(rnorm(100), rnorm(100)))
#' predict(ols_fit, cbind(rnorm(100), rnorm(100)))
predict.ols <- function(object, newdata = NULL, ...){
  # Obtain datamatrix
  if (is.null(newdata)) {
    newdata <- object$X
  } else if (object$const) {
    newdata <- cbind(1, newdata)
  }#IFELSE
  # Calculate and return fitted values with the OLS coefficient
  fitted <- newdata%*%object$coef
  return(fitted)
}#PREDICT.OLS
