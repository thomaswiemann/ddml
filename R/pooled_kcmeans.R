#' Compute pooled K conditional means estimator.
#'
#' Compute pooled K conditional means estimator.
#'
#' @param y A response vector.
#' @param X A feature vector.
#' @param K Number of conditional means.
#' @param alpha_0 K dimensional vector of initial  conditional means. When set
#'     to NULL, the KMeans++ initialization procedure is used.
#'
#' @return \code{kcmeans} returns an object of S3 class "\code{kcmeans}".
#'
#' The function \code{predict} computes fitted values for a trained model of
#'     this class.
#'
#' An object of class "\code{kcmeans}" is a list containig the following
#'     components:
#' \describe{
#' \item{\code{coef}}{A vector of least squares coefficients.}
#' \item{\code{y}}{A response vector.}
#' \item{\code{X}}{A feature matrix.}
#' \item{\code{const}}{A boolean indicating inclusion of a constant.}
#' \item{\code{w}}{An optional weights vector.}
#' }
#'
#' @examples
#'
#' @export pooled_kcmeans
pooled_kcmeans <- function(y, X,
                           K, ncluster,
                           subsample = 0.7, replace = F,
                           eps = 0, max_iter = 500,
                           silent = T) {
  # Data parameters
  nobs <- length(y)
  # Compute kcmeans estimator on subsamples
  cluster_maps <- replicate(ncluster, list(0))
  alpha_mat <- matrix(0, ncluster, K)
  for (n in 1:ncluster) {
    if (!silent) print(n)
    # Draw subsample
    sample_n <- sample(c(1:nobs), floor(nobs * subsample), replace = replace)
    # Estimate kcmeans
    fit_kcmeans <- kcmeans(y[sample_n], X[sample_n], K,
                           eps = eps, max_iter = max_iter)
    # Store clustermap and conditional means
    cluster_maps[[n]] <- fit_kcmeans$cluster_map
    alpha_mat[n, ] <- fit_kcmeans$alpha
  }#FOR
  # Organize and return output
  output <-list(cluster_maps = cluster_maps, alpha_mat = alpha_mat,
                y = y, X = X, K = K, ncluster = ncluster)
  class(output) <- "pooled_kcmeans"
  return(output)
}#POOLED_KCMEANS

# Complementary methods ========================================================
#' Predict method for pooled kcmeans.
#'
#' Predict method for pooled kcmeans.
#'
#' @export predict.pooled_kcmeans
#' @export
predict.pooled_kcmeans <- function(obj, newdata = NULL){
  # Obtain datamatrix
  if (is.null(newdata)) {
    newdata <- obj$X
  }#IF
  # Cycle trough predictions via kcmeans prediction method
  fitted <- 0
  for (n in 1:obj$ncluster) {
    fake_obj <- list(alpha = obj$alpha_mat[n,  ],
                     cluster_map = obj$cluster_maps[[n]],
                     K = obj$K, y = obj$y)
    fitted <- fitted + predict.kcmeans(fake_obj, newdata)
  }#FOR
  fitted <- fitted / obj$ncluster
  # Return fitted values
  return(fitted)
}#PREDICT.POOLED_KCMEANS

#' Instrument selection for pooled kcmeans fits.
#'
#' Instrument selection for pooled kcmeans fits. Always returns \code{TRUE}.
#'
#' @export any_iv.pooled_kcmeans
#' @export
any_iv.pooled_kcmeans <- function(obj, ...){
  TRUE
}#ANY_IV.POOLED_KCMEANS
