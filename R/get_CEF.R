#' Wrapper for [ddml::crosspred()] and [ddml::shortstacking()].
#'
#' @family utilities
#'
#' @description Wrapper for [ddml::crosspred()] and [ddml::shortstacking()].
#'
#' @inheritParams shortstacking
#' @param shortstack Boolean to use short-stacking.
#'
#' @return \code{get_CEF} returns a list containing the following components:
#'     \describe{
#'         \item{\code{oos_fitted}}{A matrix of out-of-sample predictions,
#'             each column corresponding to an ensemble type (in chronological
#'             order).}
#'         \item{\code{weights}}{An array, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedures.}
#'         \item{\code{is_fitted}}{When \code{compute_insample_predictions = T}.
#'             a list of matrices with in-sample predictions by sample fold.}
#'         \item{\code{auxilliary_fitted}}{When \code{auxilliary_X} is not
#'             \code{NULL}, a list of matrices with additional predictions.}
#'         \item{\code{oos_fitted_bylearner}}{When
#'             \code{compute_predictions_bylearner = T}, a matrix of
#'             out-of-sample predictions, each column corresponding to a base
#'             learner (in chronological order).}
#'         \item{\code{is_fitted_bylearner}}{When
#'             \code{compute_insample_predictions = T} and
#'             \code{compute_predictions_bylearner = T}, a list of matrices with
#'             in-sample predictions by sample fold.}
#'         \item{\code{auxilliary_fitted_bylearner}}{When \code{auxilliary_X} is
#'             not \code{NULL} and \code{compute_predictions_bylearner = T}, a
#'             list of matrices with additional predictions for each learner.}
#'     }
get_CEF <- function(y, X, Z = NULL,
                    learners,
                    ensemble_type,
                    shortstack,
                    compute_insample_predictions = FALSE,
                    compute_predictions_bylearner = FALSE,
                    subsamples,
                    cv_subsamples_list,
                    silent = FALSE,
                    progress = NULL,
                    auxilliary_X = NULL,
                    shortstack_y = y) {
  # Compute CEF
  if (shortstack) {
    res <- shortstacking(y, X, Z,
                         learners = learners,
                         ensemble_type = ensemble_type,
                         compute_insample_predictions =
                           compute_insample_predictions,
                         subsamples = subsamples,
                         silent = silent, progress = progress,
                         auxilliary_X = auxilliary_X,
                         shortstack_y = shortstack_y)
  } else {
    res <- crosspred(y, X, Z,
                     learners = learners,
                     ensemble_type = ensemble_type,
                     compute_insample_predictions =
                       compute_insample_predictions,
                     compute_predictions_bylearner =
                       compute_predictions_bylearner,
                     subsamples = subsamples,
                     cv_subsamples_list = cv_subsamples_list,
                     silent = silent, progress = progress,
                     auxilliary_X = auxilliary_X)
  }#IFELSE
  update_progress(silent)

  # Return estimates
  return(res)
}#GET_CEF

# Utility to print progress to console
update_progress <- function(silent) {
  if (!silent) cat(" -- Done! \n")
}#UPDATE_PROGRESS
