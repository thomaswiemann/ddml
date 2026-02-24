#' Extract Model Coefficients
#'
#' @param object An object of class \code{ddml}.
#' @param ... Currently unused.
#'
#' @return Named vector (single ensemble) or matrix
#'     (multiple ensembles).
#'
#' @examples
#' \donttest{
#' # Fit a PLM and extract coefficients
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = list(what = ols),
#'                     sample_folds = 2, silent = TRUE)
#' coef(plm_fit)
#' }
#'
#' @export
coef.ddml <- function(object, ...) {
  cf <- object$coefficients
  if (is.matrix(cf) && ncol(cf) == 1) {
    cf <- drop(cf)
  }#IF
  cf
}#COEF.DDML

#' Variance-Covariance Matrix for DDML Estimators
#'
#' @param object An object of class \code{ddml}.
#' @param ensemble_idx Integer index of the ensemble type to
#'     use. Defaults to 1 (first ensemble type).
#' @param type Character string specifying the
#'     variance-covariance estimator. One of \code{"HC1"}
#'     (default, degrees-of-freedom corrected), \code{"HC0"}
#'     (uncorrected), or \code{"HC3"} (score-based leverage
#'     adjustment for finite-sample robustness).
#' @param ... Currently unused.
#'
#' @return A p x p variance-covariance matrix.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = list(what = ols),
#'                     sample_folds = 2, silent = TRUE)
#' vcov(plm_fit)
#' vcov(plm_fit, type = "HC3")
#' }
#'
#' @export
vcov.ddml <- function(object, ensemble_idx = 1,
                      type = "HC1", ...) {
  V <- compute_ddml_variance(
    object$scores[[ensemble_idx]],
    object$J[[ensemble_idx]],
    object$cluster_variable,
    type = type)
  rownames(V) <- colnames(V) <- object$coef_names
  V
}#VCOV.DDML

#' Confidence Intervals for DDML Estimators
#'
#' @param object An object of class \code{ddml}.
#' @param parm Not used (included for generic compatibility).
#' @param level Confidence level. Default 0.95.
#' @inheritParams vcov.ddml
#' @param ... Currently unused.
#'
#' @return A matrix with columns for lower and upper bounds.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = list(what = ols),
#'                     sample_folds = 2, silent = TRUE)
#' confint(plm_fit)
#' confint(plm_fit, level = 0.90)
#' }
#'
#' @importFrom stats vcov
#' @export
confint.ddml <- function(object, parm, level = 0.95,
                         ensemble_idx = 1,
                         type = "HC1", ...) {
  cf <- object$coefficients
  if (is.matrix(cf)) {
    cf <- cf[, ensemble_idx]
  } else {
    cf <- cf[ensemble_idx]
  }#IFELSE

  V <- vcov(object, ensemble_idx = ensemble_idx,
            type = type)
  se <- sqrt(diag(V))
  z <- stats::qnorm((1 + level) / 2)
  ci <- cbind(cf - z * se, cf + z * se)
  pct <- c((1 - level) / 2, (1 + level) / 2) * 100
  colnames(ci) <- paste0(format(pct, digits = 3), " %")
  rownames(ci) <- object$coef_names
  ci
}#CONFINT.DDML

#' Summary for DDML Estimators
#'
#' @param object An object of class \code{ddml}.
#' @param type Character string specifying the
#'     variance-covariance estimator. One of \code{"HC1"}
#'     (default, degrees-of-freedom corrected), \code{"HC0"}
#'     (uncorrected), or \code{"HC3"} (score-based leverage
#'     adjustment for finite-sample robustness).
#' @param ... Currently unused.
#'
#' @return An object of class \code{summary.ddml}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = list(what = ols),
#'                     sample_folds = 2, silent = TRUE)
#' summary(plm_fit)
#' summary(plm_fit, type = "HC3")
#' }
#'
#' @export
summary.ddml <- function(object, type = "HC1", ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  single_learner <- ("what" %in% names(object$learners))
  ens_type <- if (single_learner) {
    "single base learner"
  } else {
    object$ensemble_type
  }#IFELSE

  inf <- compute_ddml_inference(
    coefficients = object$coefficients,
    scores = object$scores,
    J = object$J,
    coef_names = object$coef_names,
    ensemble_type = ens_type,
    cluster_variable = object$cluster_variable,
    type = type)

  result <- list(
    inf_results = inf,
    type = type,
    model_type = class(object)[1],
    nobs = object$nobs,
    sample_folds = object$sample_folds,
    shortstack = object$shortstack,
    ensemble_type = ens_type)
  class(result) <- c(
    paste0("summary.", class(object)[1]),
    "summary.ddml")
  result
}#SUMMARY.DDML

#' Print Summary for DDML Estimators
#'
#' @param x An object of class \code{summary.ddml}.
#' @param digits Number of significant digits. Default 3.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.summary.ddml <- function(x, digits = 3, ...) {
  type_labels <- c(
    ddml_plm = "Partially Linear Model",
    ddml_pliv = "Partially Linear IV Model",
    ddml_fpliv =
      "Flexible Partially Linear IV Model",
    ddml_ate = "Average Treatment Effect",
    ddml_att =
      "Average Treatment Effect on the Treated",
    ddml_late = "Local Average Treatment Effect")
  model_name <- type_labels[x$model_type]
  if (is.na(model_name)) model_name <- x$model_type

  cat("DDML estimation:", model_name, "\n")
  cat("Obs:", x$nobs,
      "  Folds:", x$sample_folds)
  if (!is.null(x$shortstack) && x$shortstack) {
    cat("  Stacking: short-stack")
  }#IF
  if (!is.null(x$type) && x$type != "HC1") {
    cat("  SE:", x$type)
  }#IF
  cat("\n\n")

  nensb <- dim(x$inf_results)[3]
  for (j in seq_len(nensb)) {
    if (nensb > 1) {
      cat("Ensemble type:",
          dimnames(x$inf_results)[[3]][j], "\n")
    }#IF
    tbl <- x$inf_results[, , j]
    if (!is.matrix(tbl)) {
      tbl <- matrix(tbl, nrow = 1,
                    dimnames = list(
                      dimnames(x$inf_results)[[1]],
                      dimnames(x$inf_results)[[2]]))
    }#IF
    stats::printCoefmat(tbl, digits = digits,
                        has.Pvalue = TRUE,
                        signif.stars = TRUE)
    if (j < nensb) cat("\n")
  }#FOR

  invisible(x)
}#PRINT.SUMMARY.DDML
