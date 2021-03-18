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
    alpha <- get_cmeans(y, indx_j, K, cluster_map)
    # Assign random center to empty clusters
    is_NaN <- is.nan(alpha)
    alpha[is_NaN] <- sample(alpha[!is_NaN], 1) + (1 - 2*runif(1)) # add noise
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

#' Compute K conditional means estimator via Variable Neighborhood Search.
#'
#' Compute K conditional means estimator via Variable Neighborhood Search..
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
#' @export kcmeans_vns
kcmeans_vns <- function(y, X,
                        K, alpha_0 = init_kcmeans(y, X, K), N = 2,
                        eps = 0, max_iter = 10, max_neighborhood = 10,
                        max_iter_kcmeans = 100) {
  # Data parameters
  nobs <- length(y)
  unq_x <- sort(unique(X))
  J <- length(unq_x)
  # Get subsamples
  indx_j <- lapply(unq_x, function(j) which(j == X))
  # Get initial cluster map and corresponding loss
  cluster_map <- kcmeans(y, X, K, alpha_0 = alpha_0, max_iter = 1)$cluster_map
  loss <- get_loss(y, indx_j, K, alpha_0, cluster_map)

  # Run variable neighborhood search algorithm
  alpha <- alpha_0
  for (iter in 1:max_iter) {
    n = 1
    while(n <= max_neighborhood) {
      # Relocate n randomly selected values
      rnd_j <- sample(c(1:J), n)
      cl_j <- matrix(sapply(cluster_map, function (c) !(rnd_j %in% c)), n, K)
      cluster_map_n <- cluster_map
      for (j in 1:n) {
        new_cl <- sample(c(1:K)[cl_j[j, ]], 1)
        old_cl <- which(!cl_j[j, ])
        cluster_map_n <- reassign_j(rnd_j[j], old_cl, new_cl, cluster_map)
      }#FOR

      # Obtain new cluster means
      alpha <- get_cmeans(y, indx_j, K, cluster_map_n)

      # Run kcmeans until convergence
      fit_kcmeans <- kcmeans(y, X, K, alpha_0 = alpha,
                             max_iter = max_iter_kcmeans)
      alpha <- fit_kcmeans$alpha
      cluster_map_n <- fit_kcmeans$cluster_map

      # Local search
      loss_n <- get_loss(y, indx_j, K, alpha, cluster_map_n)
      set_j <- sample(c(1:J))

      for (i in 1:J) {
        j_i <- set_j[i]
        # Find current cluster
        cl_j <- which(sapply(cluster_map, function (c) (j_i %in% c)))
        not_cl_j <- setdiff(c(1:K), cl_j)
        # Check loss associated with each reassignment
        for (k in 1:(K-1)) {
          # Reassign from cl_j to not_cl_j[k]
          cluster_map_k <- reassign_j(j_i, cl_j, not_cl_j[k], cluster_map_n)
          # Compute new cluster means for the two clusters
          alpha_k <- alpha
          alpha_k[cl_j] <- mean(y[unlist(indx_j[cluster_map_k[[cl_j]]])])
          alpha_k[not_cl_j[k]] <-
            mean(y[unlist(indx_j[cluster_map_k[[not_cl_j[k]]]])])
          # Get new loss
          loss_k <- get_loss(y, indx_j, K, alpha_k, cluster_map_k)
          # Update if improvement
          if (loss_k < loss_n) {
            loss_n <- loss_k
            alpha <- alpha_k
            cluster_map_n <- cluster_map_k
          }#IF
        }#FOR
      }#FOR

      # Check for improvement
      if (loss_n < loss) {
        cluster_map <- cluster_map_n
        loss <- loss_n
        n = 1
      } else {
        n <- n + 1
      }#IFELSE
    }#WHILE
  }#FOR

  # Organize and return output
  output <- list(alpha = alpha,
                 cluster_map = cluster_map,
                 y = y, X = X, K = K)
  class(output) <- "kcmeans"
  return(output)
}#KCMEANS_VNS

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
}#PREDICT.KCMEANS

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

#' Computation of conditional means given cluster map
get_cmeans <- function(y, indx_j, K,  cluster_map){
  alpha <- c(1:K)
  for (k in 1:K) {
    indx_jk <- unlist(indx_j[cluster_map[[k]]])
    alpha[k] <- mean(y[indx_jk])
  }#FOR
  return(alpha)
}#GET_CMEANS

#' Reassign values to different cluster
reassign_j <- function(j, old_cl, new_cl, cluster_map) {
  cluster_map[[new_cl]] <- c(cluster_map[[new_cl]], j)
  cluster_map[[old_cl]] <- setdiff(cluster_map[[old_cl]], j)
  return(cluster_map)
}#REASSIGN_J

#' Calculate current loss
get_loss <- function(y, indx_j, K, alpha, cluster_map){
  # Get fitted values
  fitted <- matrix(0, length(y), 1)
  for (k in 1:K) {
    indx_jk <- unlist(indx_j[cluster_map[[k]]])
    fitted[indx_jk] <- alpha[k]
  }#FOR
  # Get and return loss
  loss <- mean((y-fitted)^2)
  return(loss)
}#GET_LOSS
