#' Title
#'
#' @param y abc
#' @param D abc
#' @param Z abc
#' @param X abc
#' @param learners abc
#' @param learners_ZX abc
#' @param learners_DX abc
#' @param sample_folds abc
#' @param ensemble_type abc
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
ddml_pliv <- function(y, D, Z, X,
                      learners,
                      learners_ZX = learners,
                      learners_DX = learners,
                      sample_folds = 2,
                      ensemble_type = c("average"),
                      cv_folds = 5,
                      subsamples = NULL,
                      cv_subsamples_list = NULL,
                      silent = F) {
  # Data parameters
  nobs <- length(y)
  nlearners <- length(learners)

  # Draw samples if not user-supplied
  if (is.null(subsamples)) {
    subsamples <- split(c(1:nobs),
                        sample(rep(c(1:sample_folds),
                                   ceiling(nobs / sample_folds)))[1:nobs])
  }#IF
  sample_folds <- length(subsamples)

  # Compute estimates of E[y|X]
  y_X_res <- crosspred(y, X,
                       learners = learners, ensemble_type = ensemble_type,
                       cv_subsamples_list = cv_subsamples_list,
                       subsamples = subsamples,
                       compute_insample_predictions = F,
                       silent = silent, progress = "E[Y|X]: ")
  update_progress(silent)

  # Compute estimates of E[Z|X].
  Z_X_res <- crosspred(Z, X,
                       learners = learners_ZX, ensemble_type = ensemble_type,
                       cv_subsamples_list = cv_subsamples_list,
                       subsamples = subsamples,
                       compute_insample_predictions = F,
                       silent = silent, progress = "E[Z|X]: ")
  update_progress(silent)

  # Compute estimates of E[D|X].
  D_X_res <- crosspred(D, X,
                       learners = learners_DX, ensemble_type = ensemble_type,
                       cv_subsamples_list = cv_subsamples_list,
                       subsamples = subsamples,
                       compute_insample_predictions = F,
                       silent = silent, progress = "E[D|X]: ")
  update_progress(silent)

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- length(ensemble_type) > 1

  # If a single ensemble is calculated, no loops are required.
  if (!multiple_ensembles) {

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - D_X_res$oos_fitted
    V_r <- Z - Z_X_res$oos_fitted

    # Compute IV estimate with constructed variables
    iv_fit <- AER::ivreg(y_r ~ D_r | V_r)

    # Organize complementary ensemble output
    coef <- stats::coef(iv_fit)[2]
    weights <- list(y_X = y_X_res$weights,
                    D_X = D_X_res$weights,
                    Z_X = Z_X_res$weights)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (multiple_ensembles) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    nensb <- length(ensemble_type)
    coef <- matrix(0, 1, nensb)
    mspe <- iv_fit <- rep(list(1), nensb)
    nlearners <- length(learners)
    nlearners_DX <- length(learners_DX); nlearners_ZX <- length(learners_ZX)
    weights <- list(array(0, dim = c(nlearners, nensb, sample_folds)),
                    array(0, dim = c(nlearners_DX, nensb, sample_folds)),
                    array(0, dim = c(nlearners_ZX, nensb, sample_folds)))
    weights[[1]] <- y_X_res$weights; weights[[3]] <- Z_X_res$weights
    # Assign names for more legible output
    colnames(coef) <- names(mspe) <- names(iv_fit) <- ensemble_type
    names(weights) <- c("y_X", "D_X", "Z_X")
    for (j in 1:3) {
      dimnames(weights[[j]]) <- list(NULL, ensemble_type, NULL)
    }#FOR
    # Compute coefficients for each ensemble
    for (j in 1:nensb) {
      # Residualize
      y_r <- y - y_X_res$oos_fitted[, j]
      D_r <- D - D_X_res$oos_fitted[, j]
      V_r <- Z - Z_X_res$oos_fitted[, j]

      # Compute IV estimate with constructed variables
      iv_fit_j <- AER::ivreg(y_r ~ D_r | V_r)

      # Organize complementary ensemble output
      coef[j] <- stats::coef(iv_fit_j)[2]
      iv_fit[[j]] <- iv_fit_j
    }#FOR
    # Name output appropriately by ensemble type
    names(iv_fit) <- ensemble_type
  }#IF

  # Store complementary ensemble output
  mspe <- list(y_X = y_X_res$mspe,
               D_X = D_X_res$mspe,
               Z_X = Z_X_res$mspe)

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights, mspe = mspe,
                   learners = learners,
                   learners_ZX = learners_ZX,
                   learners_DX = learners_DX,
                   iv_fit = iv_fit,
                   subsamples = subsamples,
                   cv_subsamples_list = cv_subsamples_list,
                   ensemble_type = ensemble_type)

  # Amend class and return
  class(ddml_fit) <- c("ddml_pliv")
  return(ddml_fit)
}#DDML_PLIV
