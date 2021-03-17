#' Compute K conditional means estimator.
#'
#' Compute K conditional means estimator.
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
#' @export kcmeans
kcmeans <- function(y, X,
                    K, alpha_0 = init_kcmeans(y, X, K),
                    eps = 0, max_iter = 500) {
  # Data parameters
  nobs <- length(y)
  unq_x <- sort(unique(X))
  J <- length(unq_x)
  # Get subsamples
  indx_j <- lapply(unq_x, function(j) which(j == X))
  # Run K conditional means algorithm
  alpha <- alpha_0
  for (n in 1:max_iter) {
    # Obtain new clustering
    cluster_map <- rep(list(NULL), K)
    for (j in 1:J) {
      # Calculate distances
      dist_k <- sapply(c(1:K), function(k) mean((y[indx_j[[j]]] - alpha[k])^2))
      # assign to new cluster
      cluster_map[[which.min(dist_k)]] <- c(cluster_map[[which.min(dist_k)]], j)
    }#FOR
    # Obtain new cluster means
    for (k in 1:K) {
      indx_jk <- unlist(indx_j[cluster_map[[k]]])
      alpha[k] <- mean(y[indx_jk])
    }#FOR
    # Check convergence
    if (all(abs(alpha - alpha_0) <= eps)) break
    alpha_0 <- alpha
  }#FOR
  # Organize and return output
  output <- list(alpha = alpha,
                 cluster_map = cluster_map,
                 y = y, X = X, K = K)
  class(output) <- "kcmeans"
  return(output)
}#KCMEANS

# Complementary methods ========================================================
#' Predict method for kcmeans.
#'
#' Predict method for kcmeans.
#'
#' @export predict.kcmeans
#' @export
predict.kcmeans <- function(obj, newdata = NULL){
  # Obtain datamatrix
  if (is.null(newdata)) {
    newdata <- obj$X
  }#IF
  nobs <- length(newdata)
  # Calculate and return fitted values
  unq_x <- sort(unique(newdata))
  J <- length(unq_x)
  # Get subsamples
  indx_j <- lapply(unq_x, function(j) which(j == newdata))
  # Get cluster assignment
  fitted <- matrix(0, nobs, 1)
  for (k in 1:obj$K) {
    indx_jk <- unlist(indx_j[obj$cluster_map[[k]]])
    fitted[indx_jk] <- obj$alpha[k]
  }#FOR
  # Return fitted values
  return(fitted)
}#PREDICT.kcmeans

#' Instrument selection for kcmeans fits.
#'
#' Instrument selection for kcmeans fits. Always returns \code{TRUE}.
#'
#' @export any_iv.kcmeans
#' @export
any_iv.kcmeans <- function(obj, ...){
  TRUE
}#ANY_IV.KCMEANS

# Complementary functions ======================================================
#' Initial cluster means for kcmeans.
#'
#' Initial cluster means for kcmeans.
#'
#' @examples
#' # Add example here.
#'
#' @export init_kcmeans
init_kcmeans <- function(y, X, K) {
  # Data parameters
  nobs <- length(y)
  unq_x <- sort(unique(X))
  J <- length(unq_x)
  # Get subsamples
  indx_j <- lapply(unq_x, function(j) which(j == X))
  # Initialize centers, means, and distances
  centers <- alpha_0 <- NULL
  dist_c <- matrix(1, J, 1)
  # Initialize clusters
  free_j <- c(1:J)
  for (k in 1:K) {
    # Sample free center \propto dist_k^2
    j <- sample(free_j, 1, prob = dist_c^2)
    # Assign new center and corresponding mean
    centers <- c(centers, j)
    free_j <- setdiff(free_j, j)
    alpha_0 <- c(alpha_0, mean(y[indx_j[[j]]]))
    # Calculate new minimum distances
    J_c <- length(free_j)
    dist_c <- matrix(0, J_c, 1)
    for (j in 1:J_c) {
      dist_xk <- sapply(c(1:length(centers)),
                       function(k) {
                         mean((y[indx_j[[free_j[j]]]] - alpha_0[k])^2)
                       })#SAPPLY
      dist_c[j] <- min(dist_xk)
    }#FOR
  }#FOR
  # Return inital conditional means
  return(alpha_0)
}#INIT_KCMEANS
