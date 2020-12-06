#' Compute two stage least squares estimator.
#'
#' Compute two stage least squares estimator.
#'
#' @param y A response vector.
#' @param D An endogeneous variable vector or matrix.
#' @param Z An instrumental variable vector or matrix.
#' @param X A control matrix.
#'
#' @return \code{tsls} returns an object of S3 class "\code{tsls}".
#'
#' The function \code{predict} computes fitted values for a trained model of
#'     this class. The function \code{summary} computes the corresponding
#'     standard errors, t-statistics, and p-values.
#'
#' An object of class "\code{tsls}" is a list containig the following
#'     components:
#' \describe{
#' \item{\code{coef}}{A vector of least squares coefficients.}
#' \item{\code{y}}{A response vector.}
#' \item{\code{X_}}{A matrix combining inputs \code{D} and \code{X}.}
#' \item{\code{Z_}}{A matrix combining inputs \code{Z} and \code{X}.}
#' \item{\code{FS}}{The first stage coefficient matrix.}
#' }
#'
#' @examples
#' # Add example here
#'
#' @export tsls
tsls <- function(y, D, Z, X = matrix(1, nobs)) {
  # Data parameters
  nobs <- length(y)
  # Define sample matrices
  Z_ <- cbind(Z, X)
  X_ <- cbind(D, X)
  # Calculate matrix products
  ZZ <- as.matrix(crossprod(Z_))
  XZ <- crossprod(X_, Z_)
  ZY <- crossprod(Z_, y)
  FS <- tcrossprod(csolve(ZZ), XZ)
  # Calculate TSLS coefficient
  coef <- crossprod(csolve(as.matrix(XZ %*%FS)),
                    crossprod(FS, ZY))
  coef <- as.matrix(coef)
  # Organize and return output
  try(rownames(coef) <- colnames(X_)) # variable names
  output <- list(coef = coef, y = y, X_ = X_, Z_ = Z_, FS = FS)
  class(output) <- "tsls"
  return(output)
}#TSLS

# Complementary methods ========================================================
#' Predict method for tsls fits.
#'
#' Predict method for tsls fits.
#'
#' @export predict.tsls
predict.tsls <- function(obj, newdata = NULL) {
  # Obtain datamatrix
  if(is.null(newdata)) newdata <- obj$X_
  # Calculate and return fitted values with the TSLS coefficient
  fitted <- newdata %*% obj$coef
  return(fitted)
}#PREDICT.TSLS

#' Inference for tsls fits.
#'
#' Inference for tsls fits.
#'
#' @export summary.tsls
summary.tsls <- function(obj, type = "const") {
  # Data parameters
  nobs <- length(obj$y)
  ncol_X <- ncol(obj$X_); ncol_Z <- ncol(obj$Z_)
  # Calculate residuals
  resid <- (obj$y - predict(obj))[, 1]
  # Calculate matrix products
  PZ <- obj$Z_ %*% obj$FS
  PZZP_inv <- csolve(as.matrix(crossprod(obj$X_, obj$Z_) %*% obj$FS))
  if (type == "const") {
    se <- sqrt(diag(sum(resid^2) * PZZP_inv) / (nobs-ncol_Z))
  } else if (type == "HC1") {
    XuuX <- crossprod(PZ * resid^2, PZ) * (nobs/(nobs-ncol_Z))
    S1 <- PZZP_inv %*% XuuX
    se <- sqrt(diag(S1 %*% PZZP_inv)) # in two steps for numerical accuracy
  }#IF
  t_stat <- obj$coef / se
  p_val <- 2 * pnorm(-abs(t_stat))
  # Compile estimate and se
  res <- cbind(obj$coef, se, t_stat, p_val)
  rownames(res) <- rownames(coef)
  colnames(res) <- c("Coef.", "S.E.", "t Stat.", "p-val.")
  # Compute R^2
  R2 <- 1 - var(resid)/var(obj$y)
  # Organize and return output
  output <- list(res = res, nobs = nobs, R2 = R2)
  return(output)
}#SUMMARY.TSLS
