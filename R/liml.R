#' Compute limited information maximum likelihood.
#'
#' Compute limited information maximum likelihood based on Hansen, Hausman, &
#'      Newey (2008).
#'
#' @param y A response vector.
#' @param D An endogeneous variable vector or matrix.
#' @param Z An instrumental variable vector or matrix.
#' @param X A control matrix.
#'
#' @return \code{liml} returns an object of S3 class "\code{liml}".
#'
#' @examples
#' # Add example here
#'
#' @export liml
liml <- function(y, D, Z, X,
                 zero.threshold = 1e-10){
  # Combine data vectors
  Z_ <- cbind(X, Z)
  X_ <- cbind(D, X)
  # Define joint response-control vector
  W <- cbind(y, X_)
  # Calculate matrix products
  ZZ <- as.matrix(Matrix::crossprod(Z_))
  WW <- as.matrix(Matrix::crossprod(W))
  ZZ_inv <- csolve(ZZ)
  WW_inv <-csolve(WW)
  # Calculate remaining LIML matrix
  WZ <- Matrix::crossprod(W, Z_)
  M <- WW_inv %*% WZ %*% tcrossprod(ZZ_inv, WZ)

  # Calculate LIML as normalized eigenvector of the smallest non-zero eigenvalue
  eigen_res <- eigen(M)
  index <- max(which(as.numeric(eigen_res$values) > zero.threshold))
  coef <- as.numeric(eigen_res$vectors[,index])
  coef <- as.matrix((coef / coef[1] * (-1))[-1]) # normalized

  # Organize and return output
  try(rownames(coef) <- colnames(X_)) # variable names
  output <- list(coef = coef, y = y, X_ = X_, Z_ = Z_, ZZ_inv = ZZ_inv)
  class(output) <- "liml"
  return(output)
}#LIML

# Complementary methods ========================================================
#' Predict method for liml fits.
#'
#' Predict method for liml fits.
#'
#' @export predict.liml
#' @export
predict.liml <- function(obj, data = NULL){
  # Obtain datamatrix
  if(is.null(data)) data <- obj$X_
  # Calculate and return fitted values with the TSLS coefficient
  fitted <- data %*% obj$coef
  return(fitted)
}#PREDICT.LIML

#' Inference for liml fits.
#'
#' Inference for liml fits. CSE denote the standard errors in Hansen, Hausman, &
#'      Newey (2008).
#'
#' @export summary.liml
#' @export
summary.liml <- function(obj, type = "CSE") {
  # Data parameters
  nobs <- length(obj$y)
  nX <- dim(obj$X_)[2]
  nZ <- dim(obj$Z_)[2]

  # Calculate matrix products
  u = (obj$y - predict(obj))[, 1]
  uu = sum(u^2)
  ZZZ_inv = obj$Z_ %*% obj$ZZ_inv
  PX = ZZZ_inv %*% crossprod(obj$Z_, obj$X_)
  uPu = (tcrossprod(crossprod(u, ZZZ_inv), obj$Z_) %*% u / uu)[1]

  # Calculate standard errors
  H = crossprod(obj$X_, PX) - uPu  * crossprod(obj$X_)
  H_inv = csolve(H)
  sgm2 <- uu / (nobs - nX)
  if (type == "const") {
    # Standard errors under homoskedastic normal errors
    se <- sqrt(sgm2 * diag(H_inv))
  } else if (type == "CSE") {
    # Calculate additional matrix products
    X_tld <- obj$X_ - u %*% crossprod(u, obj$X_) / uu
    PX_tld <- ZZZ_inv %*% crossprod(obj$Z_, X_tld)
    V <- X_tld - PX_tld
    p_vec <- matrix(t(ZZZ_inv), nZ * nobs) *
      matrix(t(obj$Z_), nZ * nobs) # compute diags w/o redundant crossproducts
    p_vec <- colSums(matrix(p_vec, nZ, nobs))
    kappa <- sum(p_vec^2) / nZ
    tau <- nZ / nobs
    SgmB <- ((1 - uPu)^2 * crossprod(X_tld, PX_tld) +
      uPu^2 * crossprod(X_tld, V)) * sgm2
    A <- tcrossprod(colSums((p_vec - tau) * PX), colSums(u^2 * V)) / nobs
    # Ugly loop to compute B :(
    B <- 0
    for (i in 1:nobs) {
      B <- B + tcrossprod(V[i, ]) * (u[i]^2 - sgm2)
    }#FOR
    B <- B * nZ * (kappa - tau) / (nobs * (1 -  2 * tau + kappa * tau))
    # Combine and calculate corrected standard errors
    Sgm <- SgmB + A + t(A) + B
    Var <- H_inv %*% Sgm %*% H_inv
    se <- sqrt(diag(Var))
  }#IF
  t_stat <- obj$coef / se
  p_val <- 2 * pnorm(-abs(t_stat))
  # Compile estimate and se
  res <- cbind(obj$coef, se, t_stat, p_val)
  rownames(res) <- rownames(coef)
  colnames(res) <- c("Coef.", "S.E.", "t Stat.", "p-val.")
  # Compute R^2
  R2 <- 1 - var(u)/var(obj$y)
  # Organize and return output
  output <- list(res = res, nobs = nobs, R2 = R2)
  return(output)
}#SUMMARY.LIML
