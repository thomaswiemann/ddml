#' Compute (generalized) inverse of a square matrix.
#'
#' @param X A square matrix.
#'
csolve <- function(X) {
  # Attmept inversion
  X_inv <- tryCatch(solve(X), error = function(e) NA)
  # If inversion failed, calculate generalized inverse
  if (any(is.na(X_inv))) {
    X_inv <- MASS::ginv(X)
    #warning("Inversion failed. The generalized inverse was calulated instead.")
  }#IF
  #Return (generalized) inverse
  return(X_inv)
}#CSOLVE

#' Compute covariance matrix by column.
#'
#' @param X A matrix.
#' @param dof_adjust A number for degree of freedom correction.
#'
ccov <- function(X, dof_adjust = 1) {
  # Data parameters
  nrow_X <- nrow(X)
  # Calculate covariance
  colmeans_X <- Matrix::colMeans(X)
  cov_XX <- ( Matrix::crossprod(X) - nrow_X *  Matrix::tcrossprod(colmeans_X)) /
    (nrow_X - dof_adjust)
  # Return covariance
  return(cov_XX)
}#CCOV

#' Compute correlation matrix by column.
#'
#' Compute correlation matrix by column. If y is specified, the correlation
#' between X and y is calculated
#'
#' To do: If y is specified, calculate covariance only between X and y, not
#'     among Xs.
#'
#' @param X A matrix.
#' @param y An optional vector.
#'
ccor <- function(X, y = NULL) {
  # When y is specifed
  if (!is.null(y)) {
    ncol_y <- ncol(y)
    if (is.null(ncol_y)) ncol_y <- 1
    X <- cbind(y, X)
  }#IF
  # Calculate correlation matrix
  cov_XX <- ccov(X)
  sd_XX <-  Matrix::tcrossprod(sqrt(Matrix::diag(cov_XX)))
  cor_XX <- cov_XX / sd_XX
  # Select only partial correlation vector if y is not NULL
  if (!is.null(y)) {
    cor_XX <- as.matrix(cor_XX[-c(1:ncol_y), c(1:ncol_y)])
  }#IF
  # Return correlations
  return(cor_XX)
}#CCOR

#' Generate subsamples by cell.
#'
#' Generate subsamples by cell.
#'
#' @export
subsample_by_cell <- function(X, folds = 2, cell_threshold = 2) {
  # Data parameters
  nX <- ncol(X)

  # Subsample by cell
  subsamples <- rep(list(integer(0)), folds)
  for (j in 1:nX) {
    # Define cell
    cell_x <- which(X[, j] == 1)
    n_cell_x <- length(cell_x)

    # Check whether cell has sufficiently many observations
    if (n_cell_x < cell_threshold) next

    # Sample by cell
    cell_split <- split(c(1:n_cell_x),
                        sample(rep(c(1:folds),
                                   ceiling(n_cell_x / folds)))[1:n_cell_x])

    for (k in 1:folds) {
      subsamples[[k]] <- c(subsamples[[k]],
                           cell_x[cell_split[[toString(k)]]])
    }#FOR
  }#FOR

  # Return subsample list
  return(subsamples)
}#SUBSAMPLE_BY_CELL
