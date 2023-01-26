#' Title
#'
#' @param y abc
#' @param D abc
#' @param X abc
#' @param learners abc
#' @param learners_DX abc
#' @param sample_folds abc
#' @param ensemble_type abc
#' @param shortstack abc
#' @param cv_folds abc
#' @param subsamples abc
#' @param cv_subsamples_list abc
#' @param silent abc
#'
#' @return abc
#' @export
#'
#' @examples
#' 1 + 1
ddml_plm <- function(y, D, X,
                    learners,
                    learners_DX = learners,
                    sample_folds = 2,
                    ensemble_type = "average",
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
  if (shortstack) {
    y_X_res <- shortstacking(y, X,
                             learners = learners,
                             ensemble_type = ensemble_type,
                             subsamples = subsamples,
                             silent = silent, progress = "E[Y|X]: ")
  } else {
    y_X_res <- crosspred(y, X,
                         learners = learners,
                         ensemble_type = ensemble_type,
                         cv_subsamples_list = cv_subsamples_list,
                         subsamples = subsamples,
                         silent = silent, progress = "E[Y|X]: ")
  }#IFELSE
  update_progress(silent)

  # Compute estimates of E[D|X].
  if (shortstack) {
    D_X_res <- shortstacking(D, X, Z = NULL,
                             learners = learners_DX,
                             ensemble_type = ensemble_type,
                             subsamples = subsamples,
                             silent = silent, progress = "E[D|X]: ")
  } else {
    D_X_res <- crosspred(D, X, Z = NULL,
                         learners = learners_DX,
                         ensemble_type = ensemble_type,
                         cv_subsamples_list = cv_subsamples_list,
                         subsamples = subsamples,
                         silent = silent, progress = "E[D|X]: ")
  }#IFELSE
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
    # Initialize weights
    if (shortstack) {
      weights <- list(array(0, dim = c(nlearners, nensb)),
                      array(0, dim = c(nlearners_DX, nensb)))
    } else {
      weights <- list(array(0, dim = c(nlearners, nensb, sample_folds)),
                      array(0, dim = c(nlearners_DX, nensb, sample_folds)))
    }#IFELSE
    weights[[1]] <- y_X_res$weights; weights[[2]] <- D_X_res$weights
    # Assign names for more legible output
    colnames(coef) <- names(mspe) <- names(ols_fit) <- ensemble_type
    names(weights) <- c("y_X", "D_X")
    for (j in 1:2) {
      if (shortstack) {
        dimnames(weights[[j]]) <- list(NULL, ensemble_type)
      } else {
        dimnames(weights[[j]]) <- list(NULL, ensemble_type, NULL)
      }#IFELSE
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
