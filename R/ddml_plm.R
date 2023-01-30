#' Estimator for the Partially Linear Model.
#'
#' @family ddml
#'
#' @description descripton goes here.
#'
#' @details details go here.
#'
#' @param y The outcome variable.
#' @param D The endogeneous variable.
#' @param X A (sparse) matrix of control variables.
#' @param learners May take one of two forms, depending on whether a single
#'     learner or stacking with multiple learners is used for estimation of the
#'     conditional expectation function(s) \eqn{E[Y|X]} (and \eqn{E[D|X]} if
#'     \code{learners_DX} is not specified).
#'     If a single learner is used, \code{learners} is a list with two named
#'     elements:
#'     \itemize{
#'         \item{\code{what} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}.}
#'     }
#'     If stacking with multiple learners is used, \code{learners} is a list of
#'     lists, each containing four named elements:
#'     \itemize{
#'         \item{\code{fun} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}.}
#'         \item{\code{assign_X} An optional vector of indices corresponding to
#'             features in \code{X} that are passed to the corresponding base
#'             learner.}
#'     }
#' @param learners_DX Optional argument to allow for different estimators of
#'     \eqn{E[D|X]}. Setup is identical to \code{learners}.
#' @param sample_folds Number of cross-fitting folds.
#' @param ensemble_type Ensemble method to combine base learners into final
#'     estimate of the conditional expectation functions. Possible values are:
#'     \itemize{
#'         \item{\code{"nnls"} Non-negative least squares.}
#'         \item{\code{"nnls1"} Non-negative least squares with the constraint
#'             that all weights sum to one.}
#'         \item{\code{"singlebest"} Select base learner with minimum MSPE.}
#'         \item{\code{"ols"} Ordinary least squares.}
#'         \item{\code{"average"} Simple average over base learners.}
#'     }
#'     Multiple ensemble types may be passed as a vector of strings.
#' @param shortstack Boolean to use short-stacking.
#' @param cv_folds Number of folds used for cross-validation in ensemble
#'     construction.
#' @param subsamples List of vectors with sample indices for cross-fitting.
#' @param cv_subsamples_list List of lists, each corresponding to a subsample
#'     containing vectors with vectors subsample indices for cross-validation.
#' @param silent Boolean to silence estimation updates.
#'
#' @return \code{ddml_plm} returns an object of S3 class
#'     \code{ddml_plm}. An object of class \code{ddml_plm} is a list containing
#'     the following components:
#'     \describe{
#'         \item{\code{coef}}{A vector with the PLM coefficents.}
#'         \item{\code{weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{ols_fit}}{Object of class \code{lm} from the second
#'             stage regression of \eqn{Y - E[Y|X]} on \eqn{D - E[D|X]}.}
#'         \item{\code{learners},\code{learners_DX},\code{subsamples},
#'             \code{cv_subsamples_list},\code{ensemble_type}}{Pass-through of
#'             selected user-provided arguments. See above.}
#'     }
#' @export
#'
#' @examples
#' # Construct data from the included BLP_1995 data
#' y <- log(BLP_1995$share) - log(BLP_1995$outshr)
#' D <- BLP_1995$price
#' X <- as.matrix(subset(BLP_1995, select = c(air, hpwt, mpd, mpg, space)))
#' # Estimate the partially linear model using a single base learners: lasso.
#' plm_fit <- ddml_plm(y, D, X,
#'                     learners = list(what = mdl_glmnet),
#'                     sample_folds = 5,
#'                     silent = TRUE)
#' plm_fit$coef
ddml_plm <- function(y, D, X,
                    learners,
                    learners_DX = learners,
                    sample_folds = 2,
                    ensemble_type = "nnls",
                    shortstack = FALSE,
                    cv_folds = 5,
                    subsamples = NULL,
                    cv_subsamples_list = NULL,
                    silent = FALSE) {
  # Data parameters
  nobs <- length(y)

  # Create sample fold tuple
  if (is.null(subsamples)) {
    subsamples <- generate_subsamples(nobs, sample_folds)
  }#IF
  sample_folds <- length(subsamples)

  # Create cv-subsamples tuple
  if (is.null(cv_subsamples_list) & !shortstack) {
    cv_subsamples_list <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      nobs_k <- nobs - length(subsamples[[k]])
      cv_subsamples_list[[k]] <- generate_subsamples(nobs_k, cv_folds)
    }# FOR
  }#IF

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|X]
  y_X_res <- get_CEF(y, X,
                     learners = learners,
                     ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     subsamples = subsamples,
                     cv_subsamples_list = cv_subsamples_list,
                     silent = silent, progress = "E[Y|X]: ")

  # Compute estimates of E[D|X].
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX,
                     ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     subsamples = subsamples,
                     cv_subsamples_list = cv_subsamples_list,
                     silent = silent, progress = "E[D|X]: ")

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- length(ensemble_type) > 1

  # If a single ensemble is calculated, no loops are required.
  if (!multiple_ensembles) {

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - D_X_res$oos_fitted

    # Compute OLS estimate with constructed variables
    ols_fit <- stats::lm(y_r ~ D_r)

    # Organize complementary ensemble output
    coef <- stats::coef(ols_fit)[2]
    weights <- list(y_X = y_X_res$weights,
                    D_X = D_X_res$weights)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (multiple_ensembles) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    nensb <- length(ensemble_type)
    coef <- matrix(0, 1, nensb)
    mspe <- ols_fit <- rep(list(1), nensb)
    nlearners <- length(learners); nlearners_DX <- length(learners_DX)
    # Ensemble weights
    weights <- list()
    weights[[1]] <- y_X_res$weights; weights[[2]] <- D_X_res$weights
    # Assign names for more legible output
    colnames(coef) <- names(mspe) <- names(ols_fit) <- ensemble_type
    names(weights) <- c("y_X", "D_X")
    # Compute coefficients for each ensemble
    for (j in 1:nensb) {
      # Residualize
      D_r <- D - D_X_res$oos_fitted[, j]

      # Residualize y
      y_r <- y - y_X_res$oos_fitted[, j]

      # Compute OLS estimate with constructed variables
      ols_fit_j <- stats::lm(y_r ~ D_r)

      # Organize complementary ensemble output
      coef[j] <- stats::coef(ols_fit_j)[2]
      ols_fit[[j]] <- ols_fit_j
    }#FOR
    # Name output appropriately by ensemble type
    names(ols_fit) <- ensemble_type
  }#IF

  # Store complementary ensemble output
  mspe <- list(y_X = y_X_res$mspe,
               D_X = D_X_res$mspe)

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights, mspe = mspe,
                   learners = learners,
                   learners_DX = learners_DX,
                   ols_fit = ols_fit,
                   subsamples = subsamples,
                   cv_subsamples_list = cv_subsamples_list,
                   ensemble_type = ensemble_type)

  # Print estimation progress
  if (!silent) cat("DDML estimation completed. \n")

  # Amend class and return
  class(ddml_fit) <- c("ddml_plm")
  return(ddml_fit)
}#DDML_PLM
