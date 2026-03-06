#' Access Elements of a DDML Object
#'
#' Provides backward compatibility for renamed elements.
#' Accessing \code{$weights} is deprecated; use
#' \code{$ensemble_weights} instead.
#'
#' @param x An object of class \code{ddml}.
#' @param name Element name.
#'
#' @return The requested element.
#'
#' @export
`$.ddml` <- function(x, name) {
  if (name == "weights") {
    message("Note: '$weights' is deprecated for ddml ",
            "objects. Use '$ensemble_weights' instead.")
    return(.subset2(x, "ensemble_weights"))
  }#IF
  .subset2(x, name)
}#`$.DDML`

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
    nm <- rownames(cf)
    cf <- as.vector(cf)
    names(cf) <- nm
  }#IF
  cf
}#COEF.DDML

#' Extract Number of Observations
#'
#' @param object An object of class \code{ddml}.
#' @param ... Currently unused.
#'
#' @return An integer specifying the number of observations.
#'
#' @importFrom stats nobs
#' @export
nobs.ddml <- function(object, ...) {
  object$nobs
}#NOBS.DDML

#' Variance-Covariance Matrix for DDML Estimators
#'
#' @description Computes a heteroskedasticity-robust
#'     variance-covariance matrix for the DDML estimator
#'     \eqn{\hat\theta}.
#'
#' @details All implemented DDML estimators solve a moment
#' condition of the form
#' \eqn{E[m(W; \theta_0, \eta_0)] = 0} where the score
#' decomposes as
#'
#' \eqn{m(W_i; \theta, \eta) = \psi_{b}(W_i; \eta) + \psi_{a}(W_i; \eta)\,\theta.}
#'
#' Write \eqn{m_i = m(W_i; \hat\theta, \hat\eta)} for the
#' evaluated score and
#' \eqn{\hat{J} = n^{-1}\sum_i \psi_a(W_i; \hat\eta)} for
#' the sample Jacobian. Three sandwich estimators are
#' available:
#'
#' \strong{HC0}:
#' \deqn{V_{\textrm{HC0}} = \hat{J}^{-1}
#'   \left(\frac{1}{n}\sum_i m_i m_i'\right)
#'   \hat{J}^{-\top} / n}
#'
#' \strong{HC1} (default):
#' \deqn{V_{\textrm{HC1}} = V_{\textrm{HC0}}
#'   \times \frac{n}{n - p}}
#'
#' where \eqn{p} is the dimension of \eqn{\theta}.
#'
#' \strong{HC3}:
#' \deqn{V_{\textrm{HC3}} = \hat{J}^{-1}
#'   \left(\frac{1}{n}\sum_i
#'   \frac{m_i m_i'}{(1 - h_{ii})^2}\right)
#'   \hat{J}^{-\top} / n}
#'
#' where \eqn{h_{ii}} is the generalized leverage. In the
#' general case,
#' \deqn{h_{ii} = \mathrm{tr}\!\left(
#'   \nabla_\theta m(W_i; \hat\theta, \hat\eta)
#'   \left[\sum_{j=1}^n
#'   \nabla_\theta m(W_j; \hat\theta, \hat\eta)
#'   \right]^{-1}\right)}
#'
#' where \eqn{\nabla_\theta m(W_i; \hat\theta, \hat\eta)}
#' is the derivative of the \eqn{i}-th score with respect
#' to \eqn{\theta} (the nuisance parameters \eqn{\hat\eta}
#' are treated as fixed, having been estimated
#' out-of-sample). For the linear scores implemented in
#' \code{ddml}, the derivative is
#' \eqn{\nabla_\theta m(W_i; \theta, \eta) = \psi_a(W_i; \eta)},
#' so the leverage simplifies to
#'
#' \eqn{h_{ii} = \mathrm{tr}(\psi_a(W_i; \hat\eta)\,(n\hat{J})^{-1}).}
#'
#' @param object An object of class \code{ddml}.
#' @param ensemble_idx Integer index of the ensemble type to
#'     use. Defaults to 1 (first ensemble type).
#' @param type Character string specifying the
#'     variance-covariance estimator. One of \code{"HC1"}
#'     (default), \code{"HC0"}, or \code{"HC3"}.
#' @param ... Currently unused.
#'
#' @return A \eqn{p \times p}{p x p} variance-covariance matrix.
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
  h <- if (type == "HC3") {
    compute_leverage(object, ensemble_idx)
  }#IF
  V <- compute_ddml_variance(
    object$scores[[ensemble_idx]],
    object$J[[ensemble_idx]],
    object$cluster_variable,
    type = type,
    leverage = h)
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
  cf <- object$coefficients[, ensemble_idx]

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

#' Subscript a DDML Summary Object (Deprecated)
#'
#' Delegates \code{[} to the underlying inference array.
#' This method exists for backward compatibility only;
#' use \code{x$coefficients[...]} instead.
#'
#' @param x An object of class \code{summary.ddml}.
#' @param ... Indices passed to \code{[} on the inference
#'     array.
#'
#' @return The subset of the inference array.
#'
#' @export
`[.summary.ddml` <- function(x, ...) {
  message("Note: subscripting a summary.ddml object with ",
          "'[' is deprecated. Use x$coefficients[...] instead.")
  x$coefficients[...]
}#`[.SUMMARY.DDML`

#' Summary for DDML Estimators
#'
#' @description Computes a coefficient table with estimates,
#'     standard errors, t-values, and p-values for all
#'     ensemble types. Standard errors are based on a
#'     heteroskedasticity-robust sandwich variance; see
#'     \code{\link{vcov.ddml}} for the HC0/HC1/HC3 formulas.
#'
#' @param object An object of class \code{ddml}.
#' @param type Character string specifying the
#'     variance-covariance estimator. One of \code{"HC1"}
#'     (default), \code{"HC0"}, or \code{"HC3"}. See
#'     \code{\link{vcov.ddml}} for details.
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
  single_learner <- is_single_learner(object$learners)
  ens_type <- if (single_learner) {
    "single base learner"
  } else {
    object$ensemble_type
  }#IFELSE

  h_list <- if (type == "HC3") {
    lapply(seq_along(object$scores), function(j) {
      compute_leverage(object, j)
    })
  }#IF
  inf <- compute_ddml_inference(
    coefficients = object$coefficients,
    scores = object$scores,
    J = object$J,
    coef_names = object$coef_names,
    ensemble_type = ens_type,
    cluster_variable = object$cluster_variable,
    type = type,
    leverage_list = h_list)

  result <- list(
    coefficients = inf,
    type = type,
    model_type = class(object)[1],
    estimator_name = object$estimator_name,
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
  model_name <- x$estimator_name
  if (is.null(model_name)) model_name <- x$model_type

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

  nensb <- dim(x$coefficients)[3]
  for (j in seq_len(nensb)) {
    if (nensb > 1) {
      cat("Ensemble type:",
          dimnames(x$coefficients)[[3]][j], "\n")
    }#IF
    tbl <- x$coefficients[, , j]
    if (!is.matrix(tbl)) {
      tbl <- matrix(tbl, nrow = 1,
                    dimnames = list(
                      dimnames(x$coefficients)[[1]],
                      dimnames(x$coefficients)[[2]]))
    }#IF
    stats::printCoefmat(tbl, digits = digits,
                        has.Pvalue = TRUE,
                        signif.stars = TRUE)
    if (j < nensb) cat("\n")
  }#FOR

  invisible(x)
}#PRINT.SUMMARY.DDML
