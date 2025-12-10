#' Ordinary least squares.
#'
#' @family ml_wrapper
#'
#' @description Simple implementation of ordinary least squares that computes
#'     with sparse feature matrices.
#'
#' @param y The outcome variable.
#' @param X The feature matrix.
#' @param const Boolean equal to \code{TRUE} if a constant should be included.
#' @param w A vector of weights for weighted least squares.
#'
#' @return \code{ols} returns an object of S3 class
#'     \code{ols}. An object of class \code{ols} is a list containing
#'     the following components:
#'     \describe{
#'         \item{\code{coef}}{A vector with the regression coefficents.}
#'         \item{\code{y}, \code{X}, \code{const}, \code{w}}{Pass-through of the
#'             user-provided arguments. See above.}
#'     }
#' @export
#'
#' @examples
#' ols_fit <- ols(rnorm(100), cbind(rnorm(100), rnorm(100)), const = TRUE)
#' ols_fit$coef
ols <- function(y, X,
                const = TRUE,
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
  output
}#OLS

# Complementary methods ========================================================

#' @exportS3Method
predict.ols <- function(object, newdata = NULL, ...){
  # Obtain datamatrix
  if (is.null(newdata)) {
    newdata <- object$X
  } else if (object$const) {
    newdata <- cbind(1, newdata)
  }#IFELSE
  # Calculate and return fitted values with the OLS coefficient
  newdata %*% object$coef
}#PREDICT.OLS
