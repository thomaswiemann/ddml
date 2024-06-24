#' Estimator of the Mean Squared Prediction Error using Cross-Validation.
#'
#' @family utilities
#'
#' @description Estimator of the mean squared prediction error of
#'     different learners using cross-validation.
#'
#' @inheritParams ddml_plm
#' @param y The outcome variable.
#' @param X A (sparse) matrix of predictive variables.
#' @param Z Optional additional (sparse) matrix of predictive variables.
#' @param learners \code{learners} is a list of lists, each containing four
#'     named elements:
#'     \itemize{
#'         \item{\code{fun} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}.}
#'         \item{\code{assign_X} An optional vector of column indices
#'             corresponding to variables in \code{X} that are passed to
#'             the base learner.}
#'         \item{\code{assign_Z} An optional vector of column indices
#'             corresponding to variables in \code{Z} that are passed to the
#'             base learner.}
#'     }
#'     Omission of the \code{args} element results in default arguments being
#'     used in \code{fun}. Omission of \code{assign_X} (and/or \code{assign_Z})
#'     results in inclusion of all predictive variables in \code{X} (and/or
#'     \code{Z}).
#' @param cv_folds Number of folds used for cross-validation.
#' @param cv_subsamples List of vectors with sample indices for
#'     cross-validation.
#' @param progress String to print before learner and cv fold progress.
#'
#' @return \code{crossval} returns a list containing the following components:
#'     \describe{
#'         \item{\code{mspe}}{A vector of MSPE estimates,
#'             each corresponding to a base learners (in chronological order).}
#'         \item{\code{oos_resid}}{A matrix of out-of-sample prediction errors,
#'             each column corresponding to a base learners (in chronological
#'             order).}
#'         \item{\code{cv_subsamples}}{Pass-through of \code{cv_subsamples}.
#'             See above.}
#'     }
#' @export
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' X = AE98[, c("morekids", "age","agefst","black","hisp","othrace","educ")]
#'
#' # Compare ols, lasso, and ridge using 4-fold cross-validation
#' cv_res <- crossval(y, X,
#'                    learners = list(list(fun = ols),
#'                                    list(fun = mdl_glmnet),
#'                                    list(fun = mdl_glmnet,
#'                                         args = list(alpha = 0))),
#'                    cv_folds = 4,
#'                    silent = TRUE)
#' cv_res$mspe
crossval <- function(y, X, Z = NULL,
                     learners,
                     cv_folds = 5,
                     cv_subsamples = NULL,
                     silent = FALSE,
                     progress = NULL) {
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
  cv_res <- sapply(1:(cv_folds * nlearners), function(x) {
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
  oos_resid <- oos_resid[order(unlist(cv_subsamples)), , drop = FALSE]

  # Compute MSPE by learner
  mspe <- colMeans(oos_resid^2)

  # Organize and return output
  output <- list(mspe = mspe,
                 oos_resid = oos_resid,
                 cv_subsamples = cv_subsamples)
  return(output)
}#CROSSVAL

# Complementary functions ======================================================
crossval_compute <- function(test_sample, learner,
                             y, X, Z = NULL) {
  # Check whether X, Z assignment has been specified. If not, include all.
  if (is.null(learner$assign_X)) learner$assign_X <- 1:ncol(X)
  if (is.null(learner$assign_Z) & !is.null(Z)) learner$assign_Z <- 1:ncol(Z)

  # Extract model arguments
  mdl_fun <- list(what = learner$fun, args = learner$args)
  assign_X <- learner$assign_X
  assign_Z <- learner$assign_Z

  # Compute model for this fold
  #     Note: this is effectively copying the data -- improvement needed.
  mdl_fun$args$y <- y[-test_sample]
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
