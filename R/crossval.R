#' Title
#'
#' @param y abc
#' @param X abc
#' @param Z abc
#' @param learners abc
#' @param cv_folds abc
#' @param cv_subsamples abc
#' @param silent abc
#' @param progress abc
#'
#' @return object
#' @export
#'
#' @examples
#' 1 + 1
crossval <- function(y, X, Z = NULL,
                     learners,
                     cv_folds = 5,
                     cv_subsamples = NULL,
                     silent = F, progress = NULL) {
  # Data parameters
  nobs <- length(y)
  nlearners <- length(learners)

  # Create cv sample fold tuple
  if (is.null(cv_subsamples)) {
    cv_subsamples <- generate_subsamples(nobs, cv_folds)
  }#IF
  cv_folds <- length(cv_subsamples)
  nobs <- length(unlist(cv_subsamples)) # In case subsamples are user-provided

  # Compute out-of-sample errors
  cv_res <- sapply(c(1:(cv_folds * nlearners)), function(x) {
    # Select model and cv-fold for this job
    j <- ceiling(x / cv_folds) # jth model
    i <- x - cv_folds * (ceiling(x / cv_folds) - 1) # ith CV fold
    fold_x <- cv_subsamples[[i]]
    # Print progress
    if (!silent) {
      cat(paste0("\r", progress,
                 "learner ", j, "/", nlearners,
                 ", cv fold ", i, "/", cv_folds))
    }#IF
    # Compute model for this fold
    crossval_compute(test_sample = fold_x,
                     learner = learners[[j]],
                     y, X, Z)
  })#SAPPLY

  # Compile residual matrix
  oos_resid <- unlist(cv_res)
  oos_resid <- matrix(oos_resid, nobs, nlearners)

  # Organize and return output
  output <- list(oos_resid = oos_resid)
  return(output)
}#CROSSVAL

# Complementary functions ======================================================
crossval_compute <- function(test_sample, learner,
                             y, X, Z = NULL) {
  # Check whether X, Z assignment has been specified. If not, include all.
  if (is.null(learner$assign_X)) learner$assign_X <- c(1:ncol(X))
  if (is.null(learner$assign_Z) & !is.null(Z)) learner$assign_Z <- c(1:ncol(Z))

  # Extract model arguments
  mdl_fun <- list(what = learner$fun, args = learner$args)
  assign_X <- learner$assign_X
  assign_Z <- learner$assign_Z

  # Compute model for this fold
  #     Note: this is effectively copying the data -- improvement needed.
  mdl_fun$args$y <- y[-test_sample];
  mdl_fun$args$X <- cbind(X[-test_sample, assign_X, drop = F],
                          Z[-test_sample, assign_Z, drop = F])
  mdl_fit <- do.call(do.call, mdl_fun)

  # Compute out of sample residuals
  oos_fitted <- stats::predict(mdl_fit,
                               cbind(X[test_sample, assign_X, drop = F],
                                     Z[test_sample, assign_Z, drop = F]))
  oos_resid <- y[test_sample] - methods::as(oos_fitted, "matrix")

  # Return residuals and cv_Z
  return(oos_resid)
}#CROSSVAL_COMPUTE
