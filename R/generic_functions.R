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

#' Generic function to check first stage fit.
#'
#' Generic function to check first stage fit.
#'
#' @export fstage
fstage <- function(obj, ...) {
  UseMethod("fstage")
}#FSTAGE

#' Generic function for exporting fitted ddml model.
#'
#' Generic function for exporting fitted ddml model.
#'
#' @param obj A ddml fit object.
#' @param filename Location and name of the savefiles.
#'
#' @export export_ddml
export_ddml <- function(obj, filename) {
  UseMethod("export_ddml")
}#EXPORT_DDML

