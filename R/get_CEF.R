#' Title
#'
#' @param y abc
#' @param X abc
#' @param Z abc
#' @param learners abc
#' @param ensemble_type abc
#' @param shortstack abc
#' @param compute_insample_predictions abc
#' @param subsamples abc
#' @param cv_subsamples_list abc
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
get_CEF <- function(y, X, Z = NULL,
                    learners,
                    ensemble_type,
                    shortstack,
                    compute_insample_predictions = FALSE,
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
                         compute_insample_predictions = compute_insample_predictions,
                         subsamples = subsamples,
                         silent = silent, progress = progress,
                         auxilliary_X = auxilliary_X,
                         shortstack_y = shortstack_y)
  } else {
    res <- crosspred(y, X, Z,
                     learners = learners,
                     ensemble_type = ensemble_type,
                     compute_insample_predictions = compute_insample_predictions,
                     subsamples = subsamples,
                     cv_subsamples_list = cv_subsamples_list,
                     silent = silent, progress = progress,
                     auxilliary_X = auxilliary_X)
  }#IFELSE
  update_progress(silent)

  # Return estimates
  return(res)
}#GET_CEF

# Utility to print progress ot console
update_progress <- function(silent) {
  if (!silent) cat(" -- Done! \n")
}#UPDATE_PROGRESS
