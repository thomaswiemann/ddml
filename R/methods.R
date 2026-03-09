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
#' @description Extracts the estimated coefficients
#'     from a DDML model for the specified or default
#'     ensemble type.
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
#' @family ddml inference
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
#' @description Returns the number of observations
#'     used to fit the DDML model.
#'
#' @param object An object of class \code{ddml}.
#' @param ... Currently unused.
#'
#' @return An integer specifying the number of observations.
#'
#' @family ddml inference
#' @importFrom stats nobs
#' @export
nobs.ddml <- function(object, ...) {
  object$nobs
}#NOBS.DDML

#' Extract Generalized Leverage (Hat Values)
#'
#' @description Computes the generalized leverage (hat values) for a DDML 
#'     estimator. These values are used internally to compute
#'     heteroskedasticity-robust HC3 standard errors.
#'
#' @details See \code{\link{ddml-class}} for the DML framework and
#'     the definition of \eqn{\psi_a} and \eqn{\hat{J}}.
#'     For the linear scores used in \code{ddml},
#'     the generalized leverage simplifies to
#'
#'     \eqn{h_{ii} = \mathrm{tr}(\psi_a(W_i;
#'     \hat\eta)\,(n\hat{J})^{-1}).}
#'
#' @param model An object of class \code{ddml}.
#' @param ensemble_idx Integer index of the ensemble type to extract leverage
#'     values for. Defaults to 1.
#' @param ... Currently unused.
#'
#' @return A numeric vector of generalized leverage values.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = list(what = ols),
#'                     sample_folds = 2, silent = TRUE)
#' h = hatvalues(plm_fit)
#' head(h)
#' }
#' 
#' @seealso \code{\link{vcov.ddml}} for the use of leverage in HC3 standard errors.
#'
#' @family ddml inference
#' @importFrom stats hatvalues
#' @export
hatvalues.ddml <- function(model, ensemble_idx = 1, ...) {
  validate_method_args(model, ensemble_idx = ensemble_idx)
  psi_a_j <- model$psi_a[[ensemble_idx]]
  J_j <- model$J[, , ensemble_idx, drop = FALSE]
  dim(J_j) <- dim(J_j)[1:2]
  J_inv <- csolve(J_j)
  p <- ncol(J_j)
  n <- model$nobs
  dim(psi_a_j) <- c(n, p * p)
  as.vector((psi_a_j %*% as.vector(t(J_inv))) / n)
}#HATVALUES.DDML





#' Variance-Covariance Matrix for DDML Estimators
#'
#' @description Computes a heteroskedasticity-robust
#'     variance-covariance matrix for the DDML estimator
#'     \eqn{\hat\theta}.
#'
#' @details See \code{\link{ddml-class}} for the DML framework,
#'     including the definitions of the score \eqn{m_i},
#'     the Jacobian \eqn{\hat{J}}, and the base sandwich
#'     estimator \eqn{\hat\Sigma}. This function provides
#'     three variants:
#'
#' \strong{HC0}:
#' \deqn{V_{\mathrm{HC0}} = \hat{J}^{-1}
#'   \left(\frac{1}{n}\sum_i m_i m_i'\right)
#'   \hat{J}^{-\top} / n}
#'
#' \strong{HC1} (default):
#' \deqn{V_{\mathrm{HC1}} = V_{\mathrm{HC0}}
#'   \times \frac{n}{n - p}}
#'
#'     where \eqn{p} is the dimension of \eqn{\theta}.
#'
#' \strong{HC3}:
#' \deqn{V_{\mathrm{HC3}} = \hat{J}^{-1}
#'   \left(\frac{1}{n}\sum_i
#'   \frac{m_i m_i'}{(1 - h_{ii})^2}\right)
#'   \hat{J}^{-\top} / n}
#'
#'     where \eqn{h_{ii}} is the generalized leverage;
#'     see \code{\link{hatvalues.ddml}}.
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
#' @seealso \code{\link{hatvalues.ddml}},
#'     \code{\link{confint.ddml}},
#'     \code{\link{summary.ddml}}
#'
#' @family ddml inference
#' @export
vcov.ddml <- function(object, ensemble_idx = 1,
                      type = "HC1", ...) {
  type <- validate_method_args(object,
    ensemble_idx = ensemble_idx, type = type)

  h <- if (type == "HC3") {
    stats::hatvalues(object, ensemble_idx = ensemble_idx)
  }#IF

  sc <- object$scores[, , ensemble_idx, drop = FALSE]
  dim(sc) <- dim(sc)[1:2]
  J_j <- object$J[, , ensemble_idx, drop = FALSE]
  dim(J_j) <- dim(J_j)[1:2]
  p <- ncol(sc)

  # Cluster aggregation: rowsum scores to cluster level
  clustered <- !is.null(object$cluster_variable) &&
    length(unique(object$cluster_variable)) < nrow(sc)
  if (clustered) {
    sc <- rowsum(sc, object$cluster_variable)
  }#IF

  n_eff <- nrow(sc)

  if (type == "HC3") {
    if (clustered) {
      h <- as.vector(tapply(h, object$cluster_variable, sum))
    }#IF
    sc <- sc / (1 - h)
    meat <- crossprod(sc) / n_eff
    J_inv <- csolve(J_j)
    V <- J_inv %*% meat %*% t(J_inv) / n_eff
  } else if (type == "HC1") {
    meat <- crossprod(sc) / n_eff
    J_inv <- csolve(J_j)
    V <- J_inv %*% meat %*% t(J_inv) *
      n_eff / (n_eff - p) / n_eff
  } else {
    # HC0
    meat <- crossprod(sc) / n_eff
    J_inv <- csolve(J_j)
    V <- J_inv %*% meat %*% t(J_inv) / n_eff
  }#IFELSE
  
  rownames(V) <- colnames(V) <- object$coef_names
  V
}#VCOV.DDML

#' Confidence Intervals for DDML Estimators
#'
#' @description Computes confidence intervals for one or more 
#'     parameters in a fitted DDML model.
#'
#' @param object An object of class \code{ddml}.
#' @param parm A specification of which parameters are to be
#'     given confidence intervals, either a vector of numbers
#'     or a vector of names. If missing, all parameters are
#'     considered.
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
#' confint(plm_fit, parm = "D1")
#' confint(plm_fit, level = 0.90)
#' }
#' 
#' @seealso \code{\link{vcov.ddml}}
#'
#' @family ddml inference
#' @importFrom stats vcov
#' @export
confint.ddml <- function(object, parm, level = 0.95,
                         ensemble_idx = 1,
                         type = "HC1", ...) {
  cf <- object$coefficients[, ensemble_idx]
  cf_names <- object$coef_names

  if (missing(parm)) {
    parm <- cf_names
  } else if (is.numeric(parm)) {
    parm <- cf_names[parm]
  } else {
    parm <- intersect(parm, cf_names)
    if (length(parm) == 0) {
      stop("None of the specified 'parm' were found in ",
           "the model coefficients.", call. = FALSE)
    }#IF
  }#IFELSE
  
  cf <- cf[parm]

  V <- vcov(object, ensemble_idx = ensemble_idx,
            type = type)
  se <- sqrt(diag(V))[parm]
  
  z <- stats::qnorm((1 + level) / 2)
  ci <- cbind(cf - z * se, cf + z * se)
  pct <- c((1 - level) / 2, (1 + level) / 2) * 100
  colnames(ci) <- paste0(format(pct, digits = 3), " %")
  rownames(ci) <- parm
  ci
}#CONFINT.DDML

#' @rdname summary.ddml
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
#' @inheritParams vcov.ddml
#' @param ... Currently unused.
#'
#' @return An object of class \code{summary.ddml} with:
#'     \describe{
#'         \item{\code{coefficients}}{A 3-dimensional array
#'             (\eqn{p \times 4 \times}{p x 4 x} nensb) of
#'             estimates, standard errors, z-values, and
#'             p-values.}
#'         \item{\code{type}}{The HC type used.}
#'         \item{\code{nobs}}{Number of observations.}
#'         \item{\code{sample_folds}}{Number of cross-fitting
#'             folds.}
#'         \item{\code{ensemble_type}}{Ensemble type labels.}
#'     }
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
#' @seealso \code{\link{vcov.ddml}}
#'
#' @family ddml inference
#' @export
summary.ddml <- function(object, type = "HC1", ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  
  single_learner <- is_single_learner(object$learners)
  ens_type <- if (single_learner) {
    "single base learner"
  } else {
    object$ensemble_type
  }#IFELSE
  
  nensb <- length(ens_type)
  p <- nrow(object$coefficients)

  inf <- array(0, dim = c(p, 4, nensb))

  for (j in seq_len(nensb)) {
    theta_j <- object$coefficients[, j]

    V <- stats::vcov(object, ensemble_idx = j, type = type)
    se <- sqrt(diag(V))
    z_val <- theta_j / se
    p_val <- 2 * stats::pnorm(abs(z_val),
                               lower.tail = FALSE)

    inf[, 1, j] <- theta_j
    inf[, 2, j] <- se
    inf[, 3, j] <- z_val
    inf[, 4, j] <- p_val
  }#FOR

  dimnames(inf) <- list(
    object$coef_names,
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
    ens_type
  )

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

#' @rdname summary.ddml
#'
#' @param x An object of class \code{summary.ddml}.
#' @param digits Number of significant digits. Default 3.
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
