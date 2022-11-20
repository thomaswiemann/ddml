#' Compute out-of-sample residuals via cross-validation.
#'
#' Compute out-of-sample residuals via cross-validation.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param Z An optional instrument matrix.
#' @param learners A list of lists, each containing three named elements:
#'     \itemize{
#'         \item{\code{fun} The function used to trained the model. The
#'             function must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{assign_X} A vector of indices corresponding to features
#'             in \code{X} that should be used for training.}
#'         \item{\code{assign_Z} An optional vector of indices corresponding to
#'             instruments in \code{Z} that should be used for training.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}}
#'     }
#' @param cv_folds The number for cross-validation folds.
#' @param cv_subsamples An optional list of vectors, each containing indices of
#'     a test-sample. If not used-provided, the split sample folds are randomly
#'     drawn.
#' @param setup_parallel An list containing two named elements:
#'     \itemize{
#'         \item{\code{type} A string of value \code{"static"} or
#'             \code{"dynamic"}, indicating whether job scheduling should be
#'             static (via \code{parSapply}) or dynamic (via \code{foreach}).}
#'         \item{\code{cores} The number of processor units utilized. Note that
#'             if \code{cores == 1}, \code{crossval} will not be computed in
#'             parallel.
#'         }
#'     }
#' @param silent A boolean indicating whether current learners and folds should be
#'     printed to the console.
#'
#' @return \code{crossval} a list containig the following components:
#' \describe{
#' \item{\code{oos_resid}}{A matrix of out-of-sample residuals, each column
#'     corresponding to a model in \code{learners}.}
#' \item{\code{cv_Z}}{A vector of model indices corresponding to learners that
#'     selected at least one instrument in each cross-validation fold.}
#' }
#'
#' @examples
#'
#' # Simple example w/o instruments ============================================
#' # Simulate small dataset
#' X <- cbind(1, matrix(rnorm(100*99), 100, 99)) # Simulate features
#' y <- X %*% (10*runif(100) * (runif(100) < 0.05)) + rnorm(100)
#' # Define arguments
#' learners <- list(list(fun = rlasso,
#'                     args = list(include = NULL,
#'                                 iter_resid = 1, d = 5),
#'                     assign_X = c(1:ncol(X))), # rlasso
#'                list(fun = ols,
#'                     args = list(),
#'                     assign_X = c(1:ncol(X))), # ols w/ all features
#'                list(fun = ols,
#'                     args = list(),
#'                     assign_X = c(1:50)) # ols w/ some arbitrary features
#' # Compute cross-validation residuals
#' cv_res <- crossval(y, X,
#'                    learners = learners,
#'                    cv_folds = 3,
#'                    silent = T)
#' # Compute MSPE by model
#' colMeans(cv_res$oos_resid)
#'
#' @export crossval
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
#' Compute results for a single cross-validation instance.
#'
#' Compute results for a single cross-validation instance.
#'
#' @export crossval_compute
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
  oos_fitted <- predict(mdl_fit, cbind(X[test_sample, assign_X, drop = F],
                                       Z[test_sample, assign_Z, drop = F]))
  oos_resid <- y[test_sample] - as(oos_fitted, "matrix")

  # Return residuals and cv_Z
  return(oos_resid)
}#CROSSVAL_COMPUTE
