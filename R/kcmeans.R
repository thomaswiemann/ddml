#' Compute K conditional means estimator.
#'
#' Compute K conditional means estimator.
#'
#' @param y A response vector.
#' @param X A feature vector.
#' @param W A matrix of regressors.
#' @param K Number of conditional means.
#' @param alpha_0 K dimensional vector of initial  conditional means. When set
#'     to NULL, the KMeans++ initialization procedure is used.
#' @param eps Convergence tolerance.
#' @param max_iter Maximum number of iterations until algorithm is terminated.
#'
#' @return \code{kcmeans} returns an object of S3 class "\code{kcmeans}".
#'
#' The function \code{predict} computes fitted values for a trained model of
#'     this class.
#'
#' An object of class "\code{kcmeans}" is a list containig the following
#'     components:
#' \describe{
#' \item{\code{alpha}}{A vector of conditional means.}
#' \item{\code{cluster_map}}{A list of sets of indices, denoting which values
#'     of \code{X} correspond to which conditional means.}
#' \item{\code{y}}{The outcome vector.}
#' \item{\code{X}}{The feature vector.}
#' \item{\code{K}}{The number of conditional means.}
#' }
#'
#' @examples
#'
#' @export kcmeans
kcmeans <- function(y, X, K,
                    alpha_0 = NULL, beta_0 = NULL,
                    eps = 0, max_iter = 500) {
  # Check inputs
  nX <- ncol(X)
  if (!is.null(nX) && nX == 1) {
    x <- X[, 1] # convert to numeric
    W <- NULL
  } else if (!is.null(nX) && nX > 1) {
    x <- X[, 1]
    W <- X[, -1]
  } else if (is.null(nX)) {
    x <- X
    W <- NULL
  }#IFELSE

  # Data parameters
  nobs <- length(y)
  unq_x <- sort(unique(x))
  J <- length(unq_x)

  # Initialize parameters
  if (is.null(beta_0)) {
    if (is.null(W)) W <- matrix(0, nobs, 1)
    fit_ols <- ols(y, W)
    beta_0 <- fit_ols$coef
  }#IF
  if (is.null(alpha_0)) {
    resid <- y - W %*% beta_0
    alpha_0 <- init_kcmeans(resid, x, K)
  }#IF

  # Get subsamples
  indx_j <- lapply(unq_x, function(j) which(j == X))
  names(indx_j) <- unq_x

  # Run K conditional means algorithm
  alpha <- alpha_0
  beta <- beta_0
  for (n in 1:max_iter) {
    # Obtain new clustering
    cluster_map <- rep(list(NULL), K)
    for (j in 1:J) {
      # Calculate distances
      dist_k <- sapply(c(1:K), function(k) {
        resid <- y[indx_j[[j]]] - alpha[k] - W[indx_j[[j]], , drop = F] %*% beta
        mean(resid^2)
        })#SAPPLY
      # assign to new cluster
      min_k <- which.min(dist_k)[1]
      cluster_map[[min_k]] <- c(cluster_map[[min_k]], unq_x[j])
    }#FOR

    # Obtain new cluster means or intercepts
    x2 <- factor(rep(1, nobs), levels = c(1:K))
    for (k in 1:K) {
      indx_k <- unlist(indx_j[cluster_map[[k]]])
      x2[indx_k] <- k
    }#FOR
    x2_mat <- Matrix::sparse.model.matrix(~ 0 + x2)
    fit_ols <- ols(y, cbind(x2_mat, W))
    alpha <- fit_ols$coef[1:K]
    beta <- fit_ols$coef[-c(1:K)]

    # Assign random center to empty clusters
    is_empty <- alpha == 0
    alpha[is_empty] <- sample(alpha[!is_empty], sum(is_empty)) +
      (1 - 2*runif(sum(is_empty))) # add noise

    # Check convergence
    if (all(abs(c(alpha - alpha_0, beta - beta_0)) <= eps)) break
    alpha_0 <- alpha; beta_0 <- beta
  }#FOR
  # Organize and return output
  output <- list(alpha = alpha, beta = beta,
                 cluster_map = cluster_map,
                 y = y, X = cbind(x, W), K = K)
  class(output) <- "kcmeans"
  return(output)
}#KCMEANS

# Complementary methods ========================================================
#' Predict method for objects of type \code{kcmeans}.
#'
#' Predict method for objects of type \code{kcmeans}.
#'
#' @export predict.kcmeans
#' @export
predict.kcmeans <- function(obj, newdata = NULL){

  # Obtain datamatrix
  if (is.null(newdata)) {
    newdata <- obj$X
  }#IF

  # Check inputs
  nX <- ncol(newdata)
  if (!is.null(nX) && nX == 1) {
    x <- newdata[, 1] # convert to numeric
    W <- matrix(0, length(x), 1)
  } else if (!is.null(nX) && nX > 1) {
    x <- newdata[, 1]
    W <- newdata[, -1]
  } else if (is.null(nX)) {
    x <- newdata
    W <- matrix(0, length(x), 1)
  }#IFELSE
  nobs = length(x)

  # Calculate and return fitted values
  unq_x <- sort(unique(x))
  J <- length(unq_x)

  # Get subsamples
  indx_j <- lapply(unq_x, function(j) which(j == x))
  names(indx_j) <- unq_x

  # Calculate mean intercept
  mean_intercept <- mean(obj$y - obj$X[, -1, drop = F] %*% obj$beta)

  # Get cluster assignment
  fitted <- matrix(mean_intercept, nobs, 1)
  for (k in 1:obj$K) {
    # Get unq_x associated with cluster k
    is_in_k <- sapply(unq_x, function(v) any(v == obj$cluster_map[[k]]))
    # Associate observations with cluster k
    indx_jk <- unlist(indx_j[is_in_k])
    fitted[indx_jk] <- obj$alpha[k]
  }#FOR
  fitted <- fitted + W %*% obj$beta

  # Return fitted values
  return(fitted)
}#PREDICT.KCMEANS

#' Instrument selection for objects of type \code{kcmeans}.
#'
#' Instrument selection for objects of type \code{kcmeans}. Always returns
#'     \code{TRUE}.
#'
#' @export any_iv.kcmeans
#' @export
any_iv.kcmeans <- function(obj, ...){
  TRUE
}#ANY_IV.KCMEANS

# Complementary functions ======================================================
#' Initial cluster means for \code{kcmeans}.
#'
#' Initial cluster means for \code{kcmeans}.
#'
#' @examples
#' # Add example here.
#'
#' @export init_kcmeans
init_kcmeans <- function(y, X, K) {
  # Check inputs
  if (!is.null(ncol(X))) X <- X[, 1] # convert to numeric
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
