#' Compute (weighted) least squares estimator.
#'
#' Compute (weighted) least squares estimator.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param const A boolean indicating inclusion of a constant.
#' @param w An optional weights vector.
#'
#' @return \code{ols} returns an object of S3 class "\code{ols}".
#'
#' The function \code{predict} computes fitted values for a trained model of
#'     this class. The function \code{summary} computes the corresponding
#'     standard errors, t-statistics, and p-values.
#'
#' An object of class "\code{ols}" is a list containig the following components:
#' \describe{
#' \item{\code{coef}}{A vector of least squares coefficients.}
#' \item{\code{y}}{A response vector.}
#' \item{\code{X}}{A feature matrix.}
#' \item{\code{const}}{A boolean indicating inclusion of a constant.}
#' \item{\code{w}}{An optional weights vector.}
#' }
#'
#' @export ols
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
#' Predict method for ols fits.
#'
#' Predict method for ols fits.
#'
#' @export predict.ols
#' @export
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
