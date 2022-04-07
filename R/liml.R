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
  M <- WW_inv %*% WZ %*% Matrix::tcrossprod(ZZ_inv, WZ)

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
#' To do: give an option for memory-intensive CSE calculation as alternative to
#'     the memory-efficient (but runtime-intensive) current implementation.
#'
#' @export summary.liml
#' @export
summary.liml <- function(obj, type = "CSE") {
  # Data parameters
  nobs <- length(obj$y)
  nX <- dim(obj$X_)[2]
  nZ <- dim(obj$Z_)[2]

  # Calculate matrix products
  u <- (obj$y - predict(obj))[, 1]
  uu <- sum(u^2)
  XX <- Matrix::crossprod(obj$X_)
  XZ <- Matrix::crossprod(obj$X_, obj$Z_)
  uZ <- Matrix::crossprod(u, obj$Z_)
  uPu <- (Matrix::tcrossprod(uZ %*% obj$ZZ_inv, uZ) / uu)[1]
  XPX <- Matrix::tcrossprod(XZ %*% obj$ZZ_inv, XZ)

  # Calculate standard errors
  H <- as.matrix(XPX - uPu * XX)
  H_inv <- csolve(H)
  sgm2 <- uu / (nobs - nX)
  if (type == "const") {
    # Standard errors under homoskedastic normal errors
    se <- sqrt(sgm2 * Matrix::diag(H_inv))
  } else if (type == "CSE") {
    # Calculate additional matrix products
    uX <- Matrix::crossprod(u, obj$X_)
    XuuX <- Matrix::crossprod(uX)
    uPX <- Matrix::tcrossprod(uZ %*% obj$ZZ_inv, XZ)
    XuuPX <- Matrix::crossprod(uX, uPX)
    tau <- nZ / nobs
    XX_tld <- XX - XuuX / uu
    XPX_tld <- XPX + uPu * XuuX / uu - (XuuPX + Matrix::t(XuuPX)) / uu
    # Compute SgmB
    SgmB <- sgm2 *((1-uPu)^2 * XPX_tld + uPu * (XX_tld - XPX_tld))
    # Compute A and B
    p_vec <- matrix(0, nobs, 1)
    A_LHS <- A_RHS <- matrix(0, 1, nX)
    B <- matrix(0, nX, nX)
    ZuuX_ZX <- (Matrix::crossprod(uZ, uX) / uu - Matrix::t(XZ))
    pb_se <- txtProgressBar(min = 0, max = nobs, style = 3) # progress bar
    for (i in 1:nobs) {
      # Calculate observation-specific terms
      ZZZ_inv_i <- obj$Z_[i, ] %*% obj$ZZ_inv
      Upsilon_i <- Matrix::tcrossprod(ZZZ_inv_i, XZ)
      p_vec[i] <- ZZZ_inv_i %*% obj$Z_[i, ]
      V_i <- obj$X_[i, ] - u[1] * uX / uu +
        ZZZ_inv_i %*% ZuuX_ZX
      # Compute A, B
      A_LHS <- A_LHS + (p_vec[i] - tau) * Upsilon_i
      A_RHS <- A_RHS + u[i]^2 * V_i
      B <- B + (u[i]^2 - sgm2) * Matrix::crossprod(V_i)
      # Update progress bar.
      setTxtProgressBar(pb_se, i)
    }#FOR
    close(pb_se)
    A <- Matrix::crossprod(A_LHS, A_RHS) / nobs
    kappa <- sum(p_vec^2) / nobs
    B <- B * nZ * (kappa - tau) / (nobs * (1 -  2 * tau + kappa * tau))
    # Construct Sgm
    Sgm <- SgmB + A + Matrix::t(A) + B
    # Calculate corrected standard errors
    Var <- H_inv %*% Sgm %*% H_inv
    se <- sqrt(Matrix::diag(Var))
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
