# RAL: Regular Asymptotically Linear ===========================================
#
# Base class for influence-function-based inference. Estimator-agnostic:
# consumes pre-computed influence functions, does not compute them.
#
# Class hierarchy:
#   ral           — single-fit inference base
#   ddml > ral    — adds DML-specific fields (scores, J, ensemble weights)
#   lincom > ral  — linear combination of parameters
#
# See ral_rep.R for the replicated-fit counterpart.

# Constructor ==================================================================

#' Construct a RAL Inference Object
#'
#' @description Creates a regular asymptotically linear (RAL)
#'     inference object from pre-computed influence functions.
#'     This is the base class for all influence-function-based
#'     inference in \pkg{ddml}.
#'
#' @details A regular asymptotically linear (RAL) estimator
#'     \eqn{\hat\theta} satisfies
#'
#' \deqn{\hat\theta - \theta_0 = \frac{1}{n} \sum_{i=1}^{n}
#'   \phi(W_i; \theta_0) + o_p(n^{-1/2}),}
#'
#' where \eqn{\phi(W_i; \theta_0)} is the \emph{influence
#' function}. This package stores the estimated influence
#' function \eqn{\hat\phi_i \equiv \phi(W_i; \hat\theta)}
#' in the \code{inf_func} slot.
#'
#' When an observation-level derivative
#' \eqn{-n^{-1}\,\partial \hat\phi_i / \partial \theta}
#' is available (stored in \code{dinf_dtheta}), the estimator
#' supports HC3 inference via leverage; see
#' \code{\link{hatvalues.ral}}.
#'
#' The RAL framework is estimator-agnostic: it consumes
#' pre-computed influence functions and does not prescribe how
#' they are obtained. For the specific construction under
#' cross-fitting and Neyman-orthogonal scores, see
#' \code{\link{ddml-intro}}. For linear combinations of
#' \code{ddml} estimators, see \code{\link{lincom}}.
#'
#' @param coefficients A \eqn{p \times}{p x} \code{nfit} matrix
#'     of estimated coefficients. Rows are parameters, columns
#'     are fits (e.g., ensemble types).
#' @param inf_func A 3D array of dimension
#'     \eqn{n \times p \times}{n x p x} \code{nfit}. The
#'     influence function evaluated at each observation.
#' @param dinf_dtheta Optional 4D array of dimension
#'     \eqn{n \times p \times p \times}{n x p x p x}
#'     \code{nfit}. The derivative of the influence function
#'     with respect to \eqn{\theta}, used for HC3 leverage.
#'     If \code{NULL}, HC3 is unavailable.
#' @param nobs Integer number of observations.
#' @param coef_names Character vector of parameter names (length \eqn{p}).
#' @param cluster_variable Optional vector of cluster identifiers (length \eqn{n}).
#'     If non-\code{NULL}, cluster-robust inference is used.
#' @param estimator_name Character string for display.
#' @param subclass Optional character string prepended to the class vector.
#' @param ... Additional named elements stored in the object.
#'
#' @return An object of class \code{ral} (or \code{c(subclass, "ral")}).
#'
#' @export
ral <- function(coefficients,
                inf_func,
                dinf_dtheta = NULL,
                nobs,
                coef_names,
                cluster_variable = NULL,
                estimator_name = "RAL estimator",
                subclass = NULL, ...) {
  # Validate inputs ------------------------------------------------------------
  coefficients <- as.matrix(coefficients)
  p <- nrow(coefficients)
  nfit <- ncol(coefficients)

  if (!is.numeric(inf_func) || length(dim(inf_func)) != 3) {
    stop("'inf_func' must be a 3D numeric array.", call. = FALSE)
  }#IF
  if (dim(inf_func)[1] != nobs || dim(inf_func)[2] != p ||
      dim(inf_func)[3] != nfit) {
    stop("'inf_func' dimensions must be (nobs x p x nfit).", call. = FALSE)
  }#IF
  if (!is.null(dinf_dtheta)) {
    if (!is.numeric(dinf_dtheta) || length(dim(dinf_dtheta)) != 4) {
      stop("'dinf_dtheta' must be a 4D array or NULL.", call. = FALSE)
    }#IF
    if (dim(dinf_dtheta)[1] != nobs ||
        dim(dinf_dtheta)[2] != p ||
        dim(dinf_dtheta)[3] != p ||
        dim(dinf_dtheta)[4] != nfit) {
      stop("'dinf_dtheta' dimensions must be (nobs x p x p x nfit).", 
           call. = FALSE)
    }#IF
  }#IF

  # Assemble the object --------------------------------------------------------
  obj <- c(list(
    coefficients = coefficients,
    inf_func = inf_func,
    dinf_dtheta = dinf_dtheta,
    nobs = nobs,
    nfit = nfit,
    coef_names = coef_names,
    fit_labels = colnames(coefficients),
    cluster_variable = cluster_variable,
    estimator_name = estimator_name), list(...))

  cls <- if (!is.null(subclass)) c(subclass, "ral") else "ral"
  class(obj) <- cls
  obj
}#RAL

# Internal helpers =============================================================

# Validate fit_idx argument for S3 methods.
# Returns validated type invisibly.
validate_fit_idx <- function(object, fit_idx = NULL,
                             type = NULL) {
  if (!is.null(type)) {
    type <- match.arg(type, c("HC0", "HC1", "HC3"))
  }#IF
  if (!is.null(fit_idx)) {
    nfit <- dim(object$inf_func)[3]
    if (is.null(nfit)) nfit <- 1L
    if (fit_idx < 1 || fit_idx > nfit) {
      stop("fit_idx must be between 1 and ", nfit, ".", call. = FALSE)
    }#IF
  }#IF
  invisible(type)
}#VALIDATE_FIT_IDX

# S3 methods ===================================================================

#' Extract Coefficients from a RAL Object
#'
#' @param object An object inheriting from class \code{ral}.
#' @param ... Currently unused.
#'
#' @return Named vector (single fit) or matrix (multiple fits).
#'
#' @method coef ral
#' @export
coef.ral <- function(object, ...) {
  cf <- object$coefficients
  if (is.matrix(cf) && ncol(cf) == 1) {
    nm <- rownames(cf)
    cf <- as.vector(cf)
    names(cf) <- nm
  }#IF
  cf
}#COEF.RAL

#' Number of Observations in a RAL Object
#'
#' @param object An object inheriting from class \code{ral}.
#' @param ... Currently unused.
#'
#' @return Integer.
#'
#' @method nobs ral
#' @importFrom stats nobs
#' @export
nobs.ral <- function(object, ...) {
  object$nobs
}#NOBS.RAL

#' Extract leverage (Hat Values)
#'
#' @description Computes the leverage (hat values)
#'     for a RAL estimator. Used internally for HC3 standard
#'     errors.
#'
#' @details The leverage for observation \eqn{i}
#'     is
#'
#' \deqn{h_i(\theta)
#'   = \mathrm{tr}\!\left(
#'   \frac{1}{n}
#'   \frac{\partial \phi_i(\theta)}
#'   {\partial \theta}
#' \right).}
#'
#' The sample analog replaces \eqn{\phi_i(\theta)} with its
#' estimate \eqn{\hat\phi_i}: \eqn{\hat{h}_i = h_i(\hat\theta).}
#'
#' The derivative \eqn{\partial \hat\phi_i / \partial \theta}
#' is stored in the \code{dinf_dtheta} slot. For the specific
#' form of this derivative in the DML context, see
#' \code{\link{ddml-intro}}. For the leverage of linear
#' combinations, see \code{\link{lincom}}.
#'
#' @param model An object inheriting from class \code{ral}.
#' @param fit_idx Integer index of the fit to extract leverage
#'     values for. Defaults to 1.
#' @param ... Currently unused.
#'
#' @return A numeric vector of leverage values.
#'
#' @importFrom stats hatvalues
#' @method hatvalues ral
#' @export
hatvalues.ral <- function(model, fit_idx = 1, ...) {
  validate_fit_idx(model, fit_idx = fit_idx)
  if (is.null(model$dinf_dtheta)) {
    warning("hatvalues: dinf_dtheta not available; returning NA", 
            call. = FALSE)
    return(rep(NA_real_, nobs(model)))
  }#IF

  n <- model$nobs
  p <- nrow(model$coefficients)
  dinf_j <- model$dinf_dtheta[, , , fit_idx, drop = FALSE]

  h <- rep(0, n)
  for (k in seq_len(p)) h <- h + dinf_j[, k, k, 1]
  h <- h / n

  as.vector(h)
}#HATVALUES.RAL

#' Variance-Covariance Matrix for RAL Estimators
#'
#' @description Computes a heteroskedasticity-robust
#'     variance-covariance matrix.
#'
#' @details Let \eqn{\hat\phi_i} denote the estimated
#'     influence function at observation \eqn{i}. Three
#'     variance estimators are available:
#'
#' \strong{HC0}:
#' \deqn{V_{\mathrm{HC0}} = \frac{1}{n^2}\sum_i
#'   \hat\phi_i\,\hat\phi_i'}
#'
#' \strong{HC1} (default):
#' \deqn{V_{\mathrm{HC1}} = V_{\mathrm{HC0}}
#'   \times \frac{n}{n - p}}
#'
#' \strong{HC3}:
#' \deqn{V_{\mathrm{HC3}} = \frac{1}{n^2}\sum_i
#'   \frac{\hat\phi_i\,\hat\phi_i'}
#'   {(1 - \hat{h}_{\theta,i})^2}}
#'
#' where \eqn{\hat{h}_{\theta,i}} is the leverage;
#' see \code{\link{hatvalues.ral}}.
#'
#' \strong{Cluster-robust inference.} When
#' \code{cluster_variable} is non-\code{NULL} and identifies
#' fewer groups than observations, the observation-level
#' influence functions are aggregated to cluster-level
#' influence functions
#'
#' \deqn{\hat\Phi_g = \frac{G}{n} \sum_{i \in C_g} \hat\phi_i} 
#'
#' and the variance is computed as
#'
#' \deqn{V_{\mathrm{HC0}} = \frac{1}{G^2} \sum_{g=1}^{G}
#'   \hat\Phi_g\,\hat\Phi_g'.}
#'
#' @param object An object inheriting from class \code{ral}.
#' @param fit_idx Integer index of the fit. Defaults to 1.
#' @param type Character. One of \code{"HC1"} (default),
#'     \code{"HC0"}, or \code{"HC3"}.
#' @param ... Currently unused.
#'
#' @return A \eqn{p \times p}{p x p} variance-covariance
#'     matrix.
#'
#' @seealso \code{\link{hatvalues.ral}},
#'     \code{\link{confint.ral}}
#'
#' @method vcov ral
#' @importFrom stats vcov
#' @export
vcov.ral <- function(object, fit_idx = 1,
                     type = "HC1", ...) {
  type <- validate_fit_idx(object, fit_idx = fit_idx,
                           type = type)

  if_j <- object$inf_func[, , fit_idx, drop = FALSE]
  dim(if_j) <- dim(if_j)[1:2]
  p <- ncol(if_j)
  n <- nrow(if_j)

  # Cluster aggregation: rescale to cluster-level influence functions
  clustered <- !is.null(object$cluster_variable) &&
    length(unique(object$cluster_variable)) < n
  if (clustered) {
    if_j <- rowsum(if_j, object$cluster_variable)
    if_j <- if_j * (nrow(if_j) / n)
  }
  n_eff <- nrow(if_j)
  
  if (type == "HC3") {
    h <- stats::hatvalues(object, fit_idx = fit_idx)
    if (clustered) h <- as.vector(tapply(h, object$cluster_variable, sum))
    if_j <- if_j / (1 - h)
  }#IF

  V <- crossprod(if_j) / n_eff^2

  # HC1 degrees-of-freedom correction
  if (type == "HC1") V <- V * n_eff / (n_eff - p)

  rownames(V) <- colnames(V) <- object$coef_names
  V
}#VCOV.RAL

#' Confidence Intervals for RAL Estimators
#'
#' @description Computes confidence intervals for one or more
#'     parameters.
#'
#' @param object An object inheriting from class \code{ral}.
#' @param parm A specification of which parameters are to be
#'     given confidence intervals, either a vector of numbers
#'     or a vector of names. If missing, all parameters are
#'     considered.
#' @param level Confidence level. Default 0.95.
#' @param fit_idx Integer index of the fit. Defaults to 1.
#' @param type Character. HC type. Default \code{"HC1"}.
#' @param uniform Logical. If \code{TRUE}, computes uniform
#'     confidence bands using the multiplier bootstrap.
#'     Default \code{FALSE}.
#' @param bootstraps Integer number of bootstrap draws.
#'     Only used when \code{uniform = TRUE}. Default 999.
#' @param ... Currently unused.
#'
#' @return A matrix with columns for lower and upper bounds.
#'     When \code{uniform = TRUE}, the attribute
#'     \code{"crit_val"} contains the critical value.
#'
#' @references
#' Chernozhukov V, Chetverikov D, Kato K (2013). "Gaussian
#' approximations and multiplier bootstrap for maxima of sums
#' of high-dimensional random vectors." Annals of Statistics,
#' 41(6), 2786-2819.
#'
#' @seealso \code{\link{vcov.ral}}
#'
#' @method confint ral
#' @export
confint.ral <- function(object, parm = NULL, level = 0.95,
                        fit_idx = 1,
                        type = "HC1",
                        uniform = FALSE,
                        bootstraps = 999L, ...) {
  validate_fit_idx(object, fit_idx = fit_idx)
  cf <- object$coefficients[, fit_idx]
  cf_names <- object$coef_names
  names(cf) <- cf_names

  if (is.null(parm)) {
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

  V <- vcov(object, fit_idx = fit_idx, type = type)
  se_all <- sqrt(diag(V))
  se <- se_all[parm]

  if (uniform) {
    # Multiplier bootstrap
    inf_func <- object$inf_func[, , fit_idx, drop = FALSE]
    dim(inf_func) <- dim(inf_func)[1:2]
    cl <- object$cluster_variable
    n <- nrow(inf_func)
    if (!is.null(cl) && length(unique(cl)) < n) {
      inf_func <- rowsum(inf_func, cl)
      inf_func <- inf_func * (nrow(inf_func) / n)
    }#IF
    n_eff <- nrow(inf_func)
    parm_idx <- match(parm, cf_names)
    xi <- matrix(stats::rnorm(bootstraps * n_eff), bootstraps, n_eff)
    bres <- xi %*% inf_func[, parm_idx, drop = FALSE] / sqrt(n_eff)
    sigma <- se_all[parm_idx] * sqrt(n_eff)
    # Exclude degenerate components (sigma=0, e.g., reference period)
    active <- which(sigma > 0)
    if (length(active) == 0) {
      bT <- rep(0, bootstraps)
    } else {
      bT <- apply(bres[, active, drop = FALSE], 1,
                  function(b) max(abs(b / sigma[active])))
    }#IFELSE
    z <- as.numeric(stats::quantile(bT, level, type = 1, names = FALSE))
  } else {
    z <- stats::qnorm((1 + level) / 2)
  }#IFELSE
  ci <- cbind(cf - z * se, cf + z * se)
  pct <- c((1 - level) / 2, (1 + level) / 2) * 100
  colnames(ci) <- paste0(format(pct, digits = 3), " %")
  rownames(ci) <- parm
  attr(ci, "crit_val") <- z
  ci
}#CONFINT.RAL

#' Summary for RAL Estimators
#'
#' @description Computes a coefficient table with estimates,
#'     standard errors, z-values, and p-values.
#'
#' @param object An object inheriting from class \code{ral}.
#' @param type Character. HC type. Default \code{"HC1"}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{summary.ral} with:
#' \describe{
#'     \item{\code{coefficients}}{A 3-dimensional array
#'         (\eqn{p \times 4 \times}{p x 4 x} nfit).}
#'     \item{\code{type}}{The HC type used.}
#'     \item{\code{nobs}}{Number of observations.}
#' }
#'
#' @method summary ral
#' @export
summary.ral <- function(object, type = "HC1", ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))

  fit_labels <- object$fit_labels
  if (is.null(fit_labels)) {
    fit_labels <- paste0("fit", seq_len(object$nfit))
  }#IF

  nfit <- length(fit_labels)
  p <- nrow(object$coefficients)
  inf <- array(0, dim = c(p, 4, nfit))
  for (j in seq_len(nfit)) {
    theta_j <- object$coefficients[, j]

    V <- stats::vcov(object, fit_idx = j, type = type)
    se <- sqrt(diag(V))
    z_val <- theta_j / se
    p_val <- 2 * stats::pnorm(abs(z_val), lower.tail = FALSE)

    inf[, 1, j] <- theta_j
    inf[, 2, j] <- se
    inf[, 3, j] <- z_val
    inf[, 4, j] <- p_val
  }#FOR

  dimnames(inf) <- list(
    object$coef_names,
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
    fit_labels
  )

  result <- list(
    coefficients = inf,
    type = type,
    estimator_name = object$estimator_name,
    nobs = object$nobs,
    fit_labels = fit_labels)
  class(result) <- "summary.ral"
  result
}#SUMMARY.RAL

#' @rdname summary.ral
#'
#' @param x An object of class \code{summary.ral}.
#' @param digits Number of significant digits. Default 3.
#'
#' @method print summary.ral
#' @export
print.summary.ral <- function(x, digits = 3, ...) {
  name <- x$estimator_name
  if (is.null(name)) name <- "RAL estimator"

  cat("RAL estimation:", name, "\n")
  cat("Obs:", x$nobs)
  if (!is.null(x$type) && x$type != "HC1") cat("  SE:", x$type)
  cat("\n\n")

  print_coef_tables(x$coefficients, fit_label = "Fit", digits = digits)

  invisible(x)
}#PRINT.SUMMARY.RAL

#' Tidy a RAL Object
#'
#' Extracts coefficient estimates, standard errors, test
#' statistics, and p-values in a tidy data frame.
#'
#' @param x An object inheriting from class \code{ral}.
#' @param fit_idx Integer index of the fit to report.
#'     Defaults to 1. Set to \code{NULL} for all fits.
#' @param conf.int Logical. Include confidence intervals?
#'     Default \code{FALSE}.
#' @param conf.level Confidence level. Default 0.95.
#' @param type Character. HC type. Default \code{"HC1"}.
#' @param uniform Logical. Uniform confidence bands?
#'     Default \code{FALSE}.
#' @param bootstraps Integer. Bootstrap draws. Default 999.
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns \code{term},
#'     \code{estimate}, \code{std.error}, \code{statistic},
#'     \code{p.value}, and \code{fit_label}.
#'
#' @references
#' Chernozhukov V, Chetverikov D, Kato K (2013). "Gaussian
#' approximations and multiplier bootstrap for maxima of sums
#' of high-dimensional random vectors." Annals of Statistics,
#' 41(6), 2786-2819.
#'
#' @export
#' @method tidy ral
tidy.ral <- function(x, fit_idx = 1, conf.int = FALSE,
                     conf.level = 0.95,
                     type = "HC1",
                     uniform = FALSE,
                     bootstraps = 999L, ...) {
  type <- match.arg(type, c("HC1", "HC0", "HC3"))

  s <- summary(x, type = type)
  inf <- s$coefficients
  nfit <- dim(inf)[3]
  p <- dim(inf)[1]

  if (is.null(fit_idx)) {
    j_seq <- seq_len(nfit)
  } else {
    if (any(fit_idx < 1) || any(fit_idx > nfit)) {
      stop(sprintf("fit_idx must be between 1 and %d", nfit), call. = FALSE)
    }#IF
    j_seq <- fit_idx
  }#IFELSE

  # Build tidy output
  n_rows <- length(j_seq) * p
  term <- rep(dimnames(inf)[[1]], length(j_seq))
  fit_label <- rep(dimnames(inf)[[3]][j_seq], each = p)
  estimate <- std.error <- statistic <- p.value <- numeric(n_rows)
  idx <- 1
  for (j in j_seq) {
    for (k in seq_len(p)) {
      estimate[idx] <- inf[k, 1, j]
      std.error[idx] <- inf[k, 2, j]
      statistic[idx] <- inf[k, 3, j]
      p.value[idx] <- inf[k, 4, j]
      idx <- idx + 1
    }#FOR
  }#FOR

  res <- data.frame(
    term = term,
    estimate = estimate,
    std.error = std.error,
    statistic = statistic,
    p.value = p.value,
    fit_label = fit_label,
    stringsAsFactors = FALSE
  )

  if (conf.int) {
    ci_list <- lapply(j_seq, function(j) {
      stats::confint(x, fit_idx = j, level = conf.level,
                     type = type, uniform = uniform,
                     bootstraps = bootstraps)
    })
    ci_mat <- do.call(rbind, ci_list)
    res$conf.low <- as.numeric(ci_mat[, 1])
    res$conf.high <- as.numeric(ci_mat[, 2])
  }#IF

  res
}#TIDY.RAL

#' Glance at a RAL Object
#'
#' Returns a one-row summary of model-level statistics.
#'
#' @param x An object inheriting from class \code{ral}.
#' @param ... Currently unused.
#'
#' @return A one-row \code{data.frame} with columns
#'     \code{nobs} and \code{estimator_name}.
#'
#' @export
#' @method glance ral
glance.ral <- function(x, ...) {
  data.frame(
    nobs = x$nobs,
    estimator_name = if (is.null(x$estimator_name)) {
      class(x)[1]
    } else {
      x$estimator_name
    },
    stringsAsFactors = FALSE
  )
}#GLANCE.RAL

# Plot methods ================================================================

#' Plot Coefficients from a RAL Estimator
#'
#' @description Plots point estimates with confidence intervals
#'     from an object inheriting from class \code{ral}. 
#'
#' @param x An object inheriting from class \code{ral}.
#' @param parm A specification of which parameters to plot.
#'     Either a vector of names or indices. Default: all.
#' @param level Numeric. Confidence level. Default \code{0.95}.
#' @param uniform Logical. If \code{TRUE}, uses uniform
#'     confidence bands via the multiplier bootstrap. Default
#'     \code{TRUE}.
#' @param fit_idx Integer. Which fit to plot (column index of
#'     \code{coefficients}). Default \code{1}.
#' @param type Character. HC type for standard errors.
#'     Default \code{"HC1"}.
#' @param xlab Character. Label for the x-axis.
#' @param ylab Character. Label for the y-axis.
#' @param main Character. Title for the plot.
#' @param col Color for points and segments.
#'     Default \code{"black"}.
#' @param pch Point character. Default \code{19} (solid dot).
#' @param lwd Line width for confidence interval segments.
#'     Default \code{1.5}.
#' @param ... Additional arguments passed to
#'     \code{\link[graphics]{plot.default}}.
#'
#' @return Invisibly returns a list with components
#'     \code{coefficients}, \code{ci}, and \code{labels}.
#'
#' @examples
#' # Simulate a simple example
#' n <- 200
#' X <- cbind(1, stats::rnorm(n))
#' theta <- c(0.5, -0.3)
#' inf <- matrix(stats::rnorm(n * 2), n, 2)
#' obj <- ral(matrix(theta, 2, 1),
#'            array(inf, c(n, 2, 1)),
#'            nobs = n,
#'            coef_names = c("b1", "b2"))
#' plot(obj)
#'
#' @seealso \code{\link{confint.ral}}, \code{\link{summary.ral}}
#'
#' @importFrom graphics plot
#' @method plot ral
#' @export
plot.ral <- function(x, parm = NULL, level = 0.95,
                     uniform = TRUE,
                     fit_idx = 1,
                     type = "HC1",
                     xlab = NULL, ylab = NULL,
                     main = NULL,
                     col = "black", pch = 19, lwd = 1.5,
                     ...) {
  # Get coefficients and confidence intervals
  ci <- confint(x, parm = parm, level = level, fit_idx = fit_idx,
                type = type, uniform = uniform)
  cf <- x$coefficients[, fit_idx]
  labels <- rownames(ci)
  cf <- cf[labels]

  # Set up plot coordinates
  p <- length(cf)
  idx <- seq_len(p)

  # Default labels
  if (is.null(xlab)) xlab <- ""
  if (is.null(ylab)) ylab <- "Estimate"
  if (is.null(main)) {
    ci_type <- if (uniform) "uniform" else "pointwise"
    main <- paste0(format(level * 100, digits = 3), "% ", ci_type, " CI")
  }#IF

  # Plot
  ylim <- range(ci)
  graphics::plot.default(
    idx, cf, type = "n",
    xlim = c(0.5, p + 0.5), ylim = ylim,
    xaxt = "n", xlab = xlab, ylab = ylab,
    main = main, ...)
  graphics::abline(h = 0, lty = 2, col = "grey50")
  graphics::segments(idx, ci[, 1], idx, ci[, 2], col = col, lwd = lwd)
  graphics::points(idx, cf, pch = pch, col = col)
  graphics::axis(1, at = idx, labels = labels, las = 2)

  invisible(list(coefficients = cf, ci = ci, labels = labels))
}#PLOT.RAL

# List conversion =============================================================

#' Split a RAL Object by Fit
#'
#' Returns a named list of single-fit \code{ral} objects,
#'     one per column of \code{coefficients}. This is the
#'     primary mechanism for passing multi-ensemble results
#'     to \pkg{modelsummary}.
#'
#' @param x An object inheriting from class \code{ral}.
#' @param ... Currently unused.
#'
#' @return A named list of \code{ral} objects, each with
#'     \code{nfit = 1}.
#'
#' @seealso \code{\link{ral}}, \code{\link{lincom}}
#'
#' @method as.list ral
#' @export
as.list.ral <- function(x, ...) {
  nfit <- ncol(x$coefficients)
  labels <- x$fit_labels
  if (is.null(labels)) labels <- paste0("fit", seq_len(nfit))

  # Known ral fields to slice or skip
  slice_fields <- c("coefficients", "inf_func", "dinf_dtheta")
  skip_fields <- c("nfit", "fit_labels")

  out <- vector("list", nfit)
  names(out) <- labels
  for (j in seq_len(nfit)) {
    dinf_j <- if (!is.null(x$dinf_dtheta)) {
      x$dinf_dtheta[, , , j, drop = FALSE]
    }#IF
    obj <- ral(
      coefficients = x$coefficients[, j, drop = FALSE],
      inf_func = x$inf_func[, , j, drop = FALSE],
      dinf_dtheta = dinf_j,
      nobs = x$nobs,
      coef_names = x$coef_names,
      cluster_variable = x$cluster_variable,
      estimator_name = x$estimator_name,
      subclass = setdiff(class(x), "ral")[1])
    # Carry through extra fields
    extras <- setdiff(names(x),
                      c(slice_fields, skip_fields,
                        "nobs", "coef_names",
                        "cluster_variable",
                        "estimator_name"))
    for (nm in extras) {
      if (is.null(obj[[nm]])) obj[[nm]] <- x[[nm]]
    }#FOR
    out[[j]] <- obj
  }#FOR
  out
}#AS.LIST.RAL
