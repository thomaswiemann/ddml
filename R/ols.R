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
#' @examples
#' X <- matrix(rnorm(100*3), 100, 3) # Simulate features
#' y <- 1 + X %*% c(-1, 1, 0) + rnorm(100) # Simulate linear model
#' ols(y, X, const = T) # Compute least squares fit
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
predict.ols <- function(obj, newdata = NULL){
  # Obtain datamatrix
  if (is.null(newdata)) {
    newdata <- obj$X
  } else if (obj$const) {
    newdata <- cbind(1, newdata)
  }#IFELSE
  # Calculate and return fitted values with the OLS coefficient
  fitted <- newdata%*%obj$coef
  return(fitted)
}#PREDICT.OLS

#' Inference for ols fits.
#'
#' Inference for ols fits.
#'
#' To do: implement additional HC and HAC types, including clustered se.
#'     Implement se for wls.
#'
#' @export summary.ols
#' @export
summary.ols <- function(obj,
                        type = "const",
                        cluster = NULL) {
  # Data parameters
  nobs <- length(obj$y); ncol_X <- ncol(obj$X)
  calc_wls <- !is.null(obj$w)
  # Calculate standard errors, t-statistic, and p-value
  resid <- as.numeric(obj$y - predict(obj))
  if (!calc_wls) { # OLS
    XX_inv <- csolve(as.matrix(Matrix::crossprod(obj$X)))
    if (type == "const") {
      se <- sqrt(diag(sum(resid^2) * XX_inv) / (nobs - ncol_X))
    } else if (type == "HC1"){
      XuuX <- Matrix::crossprod(obj$X *(resid^2), obj$X)
      S1 <- XX_inv %*% XuuX * (nobs/(nobs - ncol_X))
      se <- sqrt(Matrix::diag(S1 %*% XX_inv))
    } else if (type == "cluster" && !is.null(cluster)) {
      cl <- unique(cluster)
      ncl <- length(cl)
      Xu_m <- t(sapply(cl, function(x, cluster, Xu){
        if(sum(cluster==x)<2){
          Xu[cluster==x,]
        } else {
          colSums(Xu[cluster==x,])
        }#IFELSE
      }, cluster = cluster, Xu = obj$X * resid))
      XuuX <- crossprod(Xu_m) * ((nobs-1)/(nobs-ncol_X)) * (ncl/(ncl-1))
      S1 <- XX_inv %*% (XuuX)
      se <- sqrt(diag(S1%*%XX_inv))
    }#IFELSE
  }#IF
  t_stat <- obj$coef / se
  p_val <- 2 * pnorm(-abs(t_stat))
  # Compile estimate and se
  res <- cbind(obj$coef, se, t_stat, p_val)
  rownames(res) <- rownames(obj$coef)
  colnames(res) <- c("Coef.", "S.E.", "t Stat.", "p-val.")
  # Compute R^2
  R2 <- 1 - var(resid) / var(obj$y)
  # Compile output
  output <- list(res = res, nobs = nobs, R2 = R2)
  return(output)
}#SUMMARY.OLS

#' Instrument selection for ols fits.
#'
#' Instrument selection for ols fits. Always returns \code{TRUE}.
#'
#' @export any_iv.ols
#' @export
any_iv.ols <- function(obj, ...){
  TRUE
}#ANY_IV.OLS
