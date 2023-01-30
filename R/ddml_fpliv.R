#' Title
#'
#' @family ddml
#'
#' @description abc
#'
#' @details abc
#'
#' @param y abc
#' @param D abc
#' @param Z abc
#' @param X abc
#' @param learners abc
#' @param learners_DXZ abc
#' @param learners_DX abc
#' @param sample_folds abc
#' @param ensemble_type abc
#' @param shortstack abc
#' @param cv_folds abc
#' @param enforce_LIE abc
#' @param subsamples abc
#' @param cv_subsamples_list abc
#' @param silent abc
#'
#' @return object
#' @export
#'
#' @examples
#' 1 + 1
ddml_fpliv <- function(y, D, Z, X,
                       learners,
                       learners_DXZ = learners,
                       learners_DX = learners,
                       sample_folds = 2,
                       ensemble_type = "average",
                       shortstack = FALSE,
                       cv_folds = 5,
                       enforce_LIE = TRUE,
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
  y_X_res <- get_CEF(y, X, Z = NULL,
                     learners = learners, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     subsamples = subsamples,
                     cv_subsamples_list = cv_subsamples_list,
                     compute_insample_predictions = F,
                     silent = silent, progress = "E[Y|X]: ")

  # Compute estimates of E[D|X,Z]. Also calculate in-sample predictions when
  #     the LIE is enforced.
  D_XZ_res <- get_CEF(D, X, Z,
                      learners = learners_DXZ, ensemble_type = ensemble_type,
                      shortstack = shortstack,
                      subsamples = subsamples,
                      cv_subsamples_list = cv_subsamples_list,
                      compute_insample_predictions = enforce_LIE,
                      silent = silent, progress = "E[D|Z,X]: ")

  # When the LIE is not enforced, estimating E[D|X] is straightforward.
  if (!enforce_LIE) {
    D_X_res <- get_CEF(D, X, Z = NULL,
                       learners = learners_DX, ensemble_type = ensemble_type,
                       shortstack = shortstack,
                       subsamples = subsamples,
                       cv_subsamples_list = cv_subsamples_list,
                       compute_insample_predictions = F,
                       silent = silent, progress = "E[D|X]: ")
  }#IF

  # Check whether multiple ensembles are computed simultaneously.
  multiple_ensembles <- length(ensemble_type) > 1

  # If a single ensemble is calculated, no loops are required.
  if (!multiple_ensembles) {
    # Check whether the law of iterated expectations (LIE) should be enforced.
    #     When the LIE is enforced (recommended), the estimates of E[D|X,Z] are
    #     used for the calculation of the estimates of E[D|X].
    if (enforce_LIE) {
      # Compute LIE-conform estimates of E[D|X]
      D_X_res <- get_CEF(D_XZ_res$is_fitted, X, Z = NULL,
                         learners = learners_DX,
                         ensemble_type = ensemble_type,
                         shortstack = shortstack,
                         cv_subsamples_list = cv_subsamples_list,
                         subsamples = subsamples,
                         compute_insample_predictions = F,
                         silent = silent, progress = "E[D|X]: ",
                         shortstack_y = D_XZ_res$oos_fitted)
    }#IFELSE

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - D_X_res$oos_fitted
    V_r <- D_XZ_res$oos_fitted - D_X_res$oos_fitted

    # Compute IV estimate with constructed variables
    iv_fit <- AER::ivreg(y_r ~ D_r | V_r)

    # Organize complementary ensemble output
    coef <- stats::coef(iv_fit)[2]
    weights <- list(y_X = y_X_res$weights,
                    D_X = D_X_res$weights,
                    D_XZ = D_XZ_res$weights)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (multiple_ensembles) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    nensb <- length(ensemble_type)
    coef <- matrix(0, 1, nensb)
    iv_fit <- rep(list(1), nensb)
    nlearners <- length(learners)
    nlearners_DX <- length(learners_DX); nlearners_DXZ <- length(learners_DXZ)
    weights <- list(array(0, dim = c(nlearners, nensb, sample_folds)),
                    array(0, dim = c(nlearners_DX, nensb, sample_folds)),
                    array(0, dim = c(nlearners_DXZ, nensb, sample_folds)))
    weights[[1]] <- y_X_res$weights; weights[[3]] <- D_XZ_res$weights
    # Assign names for more legible output
    colnames(coef) <- names(iv_fit) <- ensemble_type
    names(weights) <- c("y_X", "D_X", "D_XZ")
    # Compute coefficients for each ensemble
    for (j in 1:nensb) {
      # When the LIE is enforced, compute LIE-conform estimates of E[D|X].
      #     Otherwise use the previously calculated estimates of E[D|X].
      if (enforce_LIE) {
        progress_j <- paste0("E[D|X] (", ensemble_type[j], "): ")
        D_X_res <- get_CEF(D_XZ_res$is_fitted[[j]], X, Z = NULL,
                           learners = learners_DX,
                           ensemble_type = ensemble_type[j],
                           shortstack = shortstack,
                           cv_subsamples_list = cv_subsamples_list,
                           subsamples = subsamples,
                           compute_insample_predictions = F,
                           silent = silent,
                           progress = progress_j,
                           shortstack_y = D_XZ_res$oos_fitted[, j])
      }#IF

      # Residualize
      if (enforce_LIE) {
        D_r <- D - D_X_res$oos_fitted
        V_r <- D_XZ_res$oos_fitted[, j] - D_X_res$oos_fitted
      } else {
        D_r <- D - D_X_res$oos_fitted[, j]
        V_r <- D_XZ_res$oos_fitted[, j] - D_X_res$oos_fitted[, j]
      }#IFELSE

      # Residualize y
      y_r <- y - y_X_res$oos_fitted[, j]

      # Compute IV estimate with constructed variables
      iv_fit_j <- AER::ivreg(y_r ~ D_r | V_r)

      # Organize complementary ensemble output
      coef[j] <- stats::coef(iv_fit_j)[2]
      iv_fit[[j]] <- iv_fit_j
      if (enforce_LIE) {
        weights[[2]][, j, ] <- D_X_res$weights
      } else {
        weights[[2]][, j, ] <- D_X_res$weights[, j, ]
      }#IFELSE

    }#FOR
    # Name output appropriately by ensemble type
    names(iv_fit) <- ensemble_type
  }#IF

  # Store complementary ensemble output
  mspe <- list(y_X = y_X_res$mspe,
               D_X = D_X_res$mspe,
               D_XZ = D_XZ_res$mspe)

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights, mspe = mspe,
                   learners = learners,
                   learners_DXZ = learners_DXZ,
                   learners_DX = learners_DX,
                   iv_fit = iv_fit,
                   subsamples = subsamples,
                   cv_subsamples_list = cv_subsamples_list,
                   ensemble_type = ensemble_type,
                   enforce_LIE = enforce_LIE)

  # Print estimation progress
  if (!silent) cat("DDML estimation completed. \n")

  # Amend class and return
  class(ddml_fit) <- c("ddml_epliv")
  return(ddml_fit)
}#DDML_FPLIV
