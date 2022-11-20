#' Compute DDML PLM estimators.
#'
#' Compute DDML PLM estimators.
#'
#' @param y A response vector.
#' @param D An endogeneous variable vector.
#' @param X A control matrix.
#' @param learners May take one of two forms, depending on whether a single model
#'     should be used for residualization, or whether an ensemble procedure
#'     should be employed.
#'     If a single model should be used, \code{learners} is a list with two named
#'     elements:
#'     \itemize{
#'         \item{\code{what} The function used to trained the model. The
#'             function must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}.}
#'     }
#'     If an ensemble should be used, \code{learners} is a list of lists, each
#'     containing four named elements:
#'     \itemize{
#'         \item{\code{fun} The function used to trained the model. The
#'             function must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{assign_X} A vector of indices corresponding to features
#'             in \code{X} that should be used for training.}
#'         \item{\code{assign_Z} An optional vector of indices corresponding to
#'             instruments in \code{Z} that should be used for training.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}}
#'     }
#' @param ensemble_type A string indicating the type of ensemble. Multiple types may
#'     be passed in form of a vector of strings.
#' @param cv_folds The number for cross-validation folds.
#' @param sample_folds The number of split sample folds used for calculation of
#'     the out-of-sample predictions.
#' @param subsamples An optional list of vectors, each containing indices of
#'     a test-sample. If not used-provided, the split sample folds are randomly
#'     drawn.
#' @param setup_parallel An list containing details on the parallelization of
#'    \code{crossval}.
#' @param silent A boolean indicating whether current learners and folds should be
#'     printed to the console.
#'
#' @return \code{ddml_plm} returns an object of S3 class
#'     \code{ddml_plm}.
#'
#' An object of class \code{ddml_plm}is a list containig the
#'     following components:
#' \describe{
#' \item{\code{coef}}{A vector with the DDML PLM coefficent in the
#'     first entry.}
#'     \item{\code{...}}{...}
#' }
#'
#' @export ddml_plm
ddml_plm <- function(y, D, X,
                    learners,
                    learners_DX = learners,
                    sample_folds = 2,
                    ensemble_type = c("average"),
                    cv_folds = 5,
                    subsamples = NULL,
                    cv_subsamples_list = NULL,
                    silent = F) {
  # Data parameters
  nobs <- length(y)

  # Create sample fold tuple
  if (is.null(subsamples)) {
    subsamples <- generate_subsamples(nobs, sample_folds)
  }#IF
  sample_folds <- length(subsamples)

  # Create cv-subsamples tuple
  if (is.null(cv_subsamples_list)) {
    cv_subsamples_list <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      nobs_k <- nobs - length(subsamples[[k]])
      cv_subsamples_list[[k]] <- generate_subsamples(nobs_k, cv_folds)
    }# FOR
  }#IF

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|X]
  y_X_res <- crosspred(y, X,
                       learners = learners, ensemble_type = ensemble_type,
                       cv_subsamples_list = cv_subsamples_list,
                       subsamples = subsamples,
                       silent = silent, progress = "E[Y|X]: ")
  update_progress(silent)

  # Compute estimates of E[D|X].
  D_X_res <- crosspred(D, X, Z = NULL,
                       learners = learners_DX, ensemble_type = ensemble_type,
                       cv_subsamples_list = cv_subsamples_list,
                       subsamples = subsamples,
                       silent = silent, progress = "E[D|X]: ")
  update_progress(silent)

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
    weights <- list(array(0, dim = c(nlearners, nensb, sample_folds)),
                    array(0, dim = c(nlearners_DX, nensb, sample_folds)))
    weights[[1]] <- y_X_res$weights; weights[[2]] <- D_X_res$weights
    # Assign names for more legible output
    colnames(coef) <- names(mspe) <- names(ols_fit) <- ensemble_type
    names(weights) <- c("y_X", "D_X")
    for (j in 1:2) {
      dimnames(weights[[j]]) <- list(NULL, ensemble_type, NULL)
    }#FOR
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
      weights[[2]][, j, ] <- D_X_res$weights[, j, ]
    }#FOR
    # Name output appropriately by ensemble type
    names(ols_fit) <- ensemble_type
  }#IF

  # Store complementary ensemble output
  mspe <- list(D_X = D_X_res$mspe,
               y_X = y_X_res$mspe)

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
