#' Title
#'
#' @param y abc
#' @param X abc
#' @param Z abc
#' @param learners abc
#' @param sample_folds abc
#' @param ensemble_type abc
#' @param cv_folds abc
#' @param compute_insample_predictions abc
#' @param subsamples abc
#' @param cv_subsamples_list abc
#' @param silent abc
#' @param progress abc
#' @param auxilliary_X abc
#'
#' @return object
#' @export
#'
#' @examples
#' 1 + 1
shortstacking <- function (y, X, Z = NULL,
                           learners = learners,
                           ensemble_type,
                           subsamples,
                           silent, progress) {

  # Data parameters
  nlearners <- length(learners)
  nensb <- length(ensemble_type)

  # Compute out-of-sample predictions for each learner
  res <- crosspred(y, X, Z,
                   learners = learners,
                   ensemble_type = "average",
                   compute_predictions_bylearner = TRUE,
                   subsamples = subsamples,
                   silent = silent, progress = progress)

  # Compute ensemble_weights
  fakecv <- list()
  fakecv$oos_resid <- kronecker(y, t(rep(1, nlearners))) -
    res$oos_fitted_bylearner
  weights <- ensemble_weights(y, X, learners = learners,
                              type = ensemble_type,
                              cv_results = fakecv)$weights

  # Compute predictions
  oos_fitted <- res$oos_fitted_bylearner %*% weights

  # Compute mspe
  mspe <- colMeans((kronecker(y, t(rep(1, nensb))) - oos_fitted)^2)
  names(mspe) <- ensemble_type

  # return shortstacking output
  return(list(oos_fitted = oos_fitted, weights = weights, mspe = mspe))
}#SHORTSTACKING
