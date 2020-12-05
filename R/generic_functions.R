#' Generic function to check for instrument selection of a trained model.
#'
#' Generic function to check for instrument selection of a trained model.
#'
#' @param obj A trained model.
#' @param ... Additional arguments for the identification of the instruments.
#'     In most cases, the user should pass a vector of indices corresponding to
#'     instruments among the features used by \code{obj} for prediction. When
#'     \code{class(obj) == "boosting"}, the user should also pass a vector of
#'     instrument names.
#'
#' @return Returns \code{TRUE} if at least one instrument is used for
#'     prediction, and \code{FALSE} otherwise.
#'
#' @export any_iv
any_iv <- function(obj, ...) {
  UseMethod("any_iv")
}#ANY_IV
