#' Title
#'
#' @param y abc
#' @param X abc
#' @param Z abc
#' @param learners abc
#' @param sample_folds abc
#' @param ensemble_type abc
#' @param compute_insample_predictions abc
#' @param subsamples abc
#' @param silent abc
#' @param progress abc
#' @param auxilliary_X abc
#' @param shortstack_y abc
#'
#' @return object
#' @export
#'
#' @examples
#' 1 + 1
shortstacking <- function (y, X, Z = NULL,
                           learners,
                           sample_folds = 2,
                           ensemble_type,
                           compute_insample_predictions = FALSE,
                           subsamples = NULL,
                           silent = FALSE,
                           progress = NULL,
                           auxilliary_X = NULL,
                           shortstack_y = y) {

  # Data parameters
  nobs <- nrow(X)
  nlearners <- length(learners)
  calc_ensemble <- !("what" %in% names(learners))
  # Create sample fold tuple
  if (is.null(subsamples)) {
    subsamples <- generate_subsamples(nobs, sample_folds)
  }#IF
  sample_folds <- length(subsamples)
  nensb <- length(ensemble_type)

  # Compute out-of-sample predictions for each learner
  res <- crosspred(y, X, Z,
                   learners = learners,
                   ensemble_type = "average",
                   compute_insample_predictions = compute_insample_predictions,
                   compute_predictions_bylearner = TRUE,
                   subsamples = subsamples,
                   silent = silent, progress = progress,
                   auxilliary_X = auxilliary_X)

  # Compute ensemble weights via subsample cross-fitted residual
  fakecv <- list()
  fakecv$oos_resid <- kronecker(shortstack_y, t(rep(1, nlearners))) -
    res$oos_fitted_bylearner
  weights <- ensemble_weights(shortstack_y, X, learners = learners,
                              type = ensemble_type,
                              cv_results = fakecv)$weights

  # Compute predictions
  oos_fitted <- res$oos_fitted_bylearner %*% weights

  # Compute auxilliary predictions (optional)
  auxilliary_fitted <- rep(list(NULL), sample_folds)
  if (!is.null(auxilliary_X)) {
    for (k in 1:sample_folds) {
      auxilliary_fitted[[k]] <- res$auxilliary_fitted_bylearner[[k]] %*% weights
    }#FOR
  }#if

  # Compute in-sample predictions (optional)
  is_fitted <- rep(list(NULL), sample_folds)
  fakecv_k <- list()
  if (compute_insample_predictions) {
    for (k in 1:sample_folds) {
      # Compute shortstacking weights in-sample
      fakecv_k$oos_resid <- kronecker(y[-subsamples[[k]]], t(rep(1, nlearners))) -
        res$is_fitted_bylearner[[k]]
      weights_k <- ensemble_weights(y[-subsamples[[k]]], X[-subsamples[[k]], ],
                                    learners = learners,
                                    type = ensemble_type,
                                    cv_results = fakecv_k)$weights
      # Combine base learners
      is_fitted[[k]] <- res$is_fitted_bylearner[[k]] %*% weights_k
    }#FOR

    # When multiple ensembles are computed, need to reorganize is_fitted
    if (nensb > 1) {
      # Loop over each ensemble type to creat list of is_fitted's
      new_is_fitted <- rep(list(rep(list(1), sample_folds)), nensb)
      for (i in 1:nensb) {
        for (k in 1:sample_folds) {
          new_is_fitted[[i]][[k]] <- is_fitted[[k]][, i, drop = F]
        }#FOR
      }#FOR
      is_fitted <- new_is_fitted
    }#IF
  }#IF

  # Compute mspe
  mspe <- colMeans((kronecker(shortstack_y, t(rep(1, nensb))) - oos_fitted)^2)

  # Assign names
  colnames(weights) <- names(mspe) <- ensemble_type

  # return shortstacking output
  output <- list(oos_fitted = oos_fitted, weights = weights, mspe = mspe,
                 is_fitted = is_fitted, auxilliary_fitted = auxilliary_fitted)
  return(output)
}#SHORTSTACKING
