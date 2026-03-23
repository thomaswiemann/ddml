# RAL_REP: Replicated RAL Inference ============================================
#
# Base class for multi-resample inference. Wraps a list of ral objects
# and aggregates coefficients/covariance via inflate-then-median (or
# mean, spectral). Estimator-agnostic.
#
# Class hierarchy:
#   ral_rep               — replicated inference base
#   ddml_rep > ral_rep    — adds DML-specific display
#   lincom_rep > ral_rep  — adds lincom-specific display

# Internal helpers =============================================================

# Spectral-norm median of PSD matrices via SDP (CVXR).
#
# Finds V* = argmin_{V >= 0} sum_r ||V - V_r||_op where
# ||.||_op is the spectral norm (largest singular value).
# For p = 1 this reduces to the standard scalar median.
spectral_median_psd <- function(matrices) {
  p <- nrow(matrices[[1]])
  R <- length(matrices)

  if (p == 1) {
    vals <- vapply(matrices, function(m) m[1, 1], numeric(1))
    return(matrix(stats::median(vals), 1, 1))
  }#IF

  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' is required for spectral ",
         "aggregation. Install it with:\n",
         "  install.packages('CVXR')", call. = FALSE)
  }#IF

  V <- CVXR::Variable(c(p, p), PSD = TRUE)
  obj <- 0
  for (r in seq_len(R)) obj <- obj + CVXR::norm(V - matrices[[r]], "2")

  prob <- CVXR::Problem(CVXR::Minimize(obj))
  result <- CVXR::psolve(prob)

  # CVXR 1.8.x (S7): psolve may return a plain numeric (optimal value)
  # CVXR < 1.8 (S4): psolve returns an object with $status and $getValue
  if (is.atomic(result)) {
    # S7: variable value populated after solve
    V_sol <- tryCatch(CVXR::value(V), error = function(e) V@value)
  } else {
    status <- if (is.list(result)) result$status else result@status
    if (!identical(status, "optimal")) {
      stop("SDP solver returned status '", status, "'.", call. = FALSE)
    }#IF
    V_sol <- if (is.list(result)) result$getValue(V) else result@getValue(V)
  }#IFELSE
  (V_sol + t(V_sol)) / 2
}#SPECTRAL_MEDIAN_PSD

# Core aggregation workhorse.
#
# For each replication r the inflated covariance is
#   V_r = Sigma_r + (theta_r - theta_tilde)(theta_r - theta_tilde)'
# where theta_tilde is the aggregated coefficient vector.
# The three aggregation rules differ only in how they
# summarise {V_1, ..., V_R} into a single matrix.
aggregate_reps <- function(object, aggregation = "median", type = "HC1") {
  aggregation <- match.arg(aggregation, c("median", "mean", "spectral"))
  R <- object$nresamples
  nfit <- object$nfit
  p <- length(object$coef_names)

  # Collect per-replication estimates
  coef_array <- array(0, dim = c(p, nfit, R))
  vcov_array <- array(0, dim = c(p, p, nfit, R))
  for (r in seq_len(R)) {
    fit <- object$fits[[r]]
    coef_array[, , r] <- fit$coefficients
    for (j in seq_len(nfit)) {
      vcov_array[, , j, r] <- stats::vcov(fit, fit_idx = j, type = type)
    }#FOR
  }#FOR

  # Aggregate coefficients
  if (aggregation == "mean") {
    agg_coef <- apply(coef_array, c(1, 2), mean)
  } else {
    agg_coef <- apply(coef_array, c(1, 2), stats::median)
  }#IFELSE

  # Inflate & aggregate covariance
  agg_vcov <- array(0, dim = c(p, p, nfit))
  for (j in seq_len(nfit)) {
    V_list <- vector("list", R)
    for (r in seq_len(R)) {
      bdiff <- coef_array[, j, r] - agg_coef[, j]
      Sigma_r <- matrix(vcov_array[, , j, r], p, p)
      V_list[[r]] <- Sigma_r + tcrossprod(bdiff)
    }#FOR

    if (aggregation == "mean") {
      V_sum <- Reduce(`+`, V_list)
      agg_vcov[, , j] <- V_sum / R
    } else if (aggregation == "spectral") {
      agg_vcov[, , j] <- spectral_median_psd(V_list)
    } else {
      V_arr <- array(unlist(lapply(V_list, as.vector)), dim = c(p, p, R))
      agg_vcov[, , j] <- apply(V_arr, c(1, 2), stats::median)
    }#IFELSE
  }#FOR

  # Standard errors
  agg_se <- matrix(0, nrow = p, ncol = nfit)
  for (j in seq_len(nfit)) {
    V_j <- matrix(agg_vcov[, , j, drop = FALSE], nrow = p, ncol = p)
    agg_se[, j] <- sqrt(diag(V_j))
  }#FOR

  list(coefficients = agg_coef, se = agg_se,
       vcov = agg_vcov, coef_array = coef_array,
       vcov_array = vcov_array)
}#AGGREGATE_REPS

# Constructor ==================================================================

#' Construct a Replicated RAL Inference Object
#'
#' @description Creates a replicated RAL inference object
#'     from a list of \code{ral} objects. Provides
#'     cross-resample aggregation for coefficients and
#'     covariance matrices.
#'
#' @param fits A list of at least 2 objects inheriting from
#'     class \code{"ral"}. All fits must share the same
#'     coefficient names and number of observations.
#' @param subclass Optional character string prepended to
#'     the class vector.
#' @param ... Additional named elements stored in the object.
#'
#' @return An object of class \code{"ral_rep"} (or
#'     \code{c(subclass, "ral_rep")}).
#'
#' @export
ral_rep <- function(fits, subclass = NULL, ...) {
  # Input validation
  if (!is.list(fits) || length(fits) < 2) {
    stop("'fits' must be a list of at least 2 ral objects.", call. = FALSE)
  }#IF
  for (i in seq_along(fits)) {
    if (!inherits(fits[[i]], "ral")) {
      stop("Element ", i, " does not inherit from class 'ral'.", call. = FALSE)
    }#IF
  }#FOR

  ref <- fits[[1]]
  for (i in seq_along(fits)[-1]) {
    if (!identical(fits[[i]]$coef_names, ref$coef_names)) {
      stop("Fit ", i, " has different 'coef_names' than fit 1.", call. = FALSE)
    }#IF
    if (!identical(fits[[i]]$nobs, ref$nobs)) {
      stop("Fit ", i, " has different 'nobs' than fit 1.", call. = FALSE)
    }#IF
  }#FOR

  obj <- c(list(
    fits = fits,
    nresamples = length(fits),
    nfit = ref$nfit,
    coef_names = ref$coef_names,
    fit_labels = ref$fit_labels,
    estimator_name = ref$estimator_name,
    nobs = ref$nobs), list(...))

  cls <- if (!is.null(subclass)) {
    c(subclass, "ral_rep")
  } else {
    "ral_rep"
  }#IFELSE
  class(obj) <- cls
  obj
}#RAL_REP

# S3 methods ===================================================================

#' @method [[ ral_rep
#' @export
`[[.ral_rep` <- function(x, i) x$fits[[i]]

#' @method length ral_rep
#' @export
length.ral_rep <- function(x) x$nresamples

#' @method nobs ral_rep
#' @importFrom stats nobs
#' @export
nobs.ral_rep <- function(object, ...) object$nobs

#' @method print ral_rep
#' @export
print.ral_rep <- function(x, ...) {
  cat("RAL replicated fits:", x$estimator_name, "\n")
  cat("  Resamples:", x$nresamples, "  Obs:", x$nobs, "\n\n")
  cat("Use summary() for aggregated inference.\n")
  cat("Use x[[i]] to access individual fits.\n")
  invisible(x)
}#PRINT.RAL_REP

#' Extract Aggregated Coefficients
#'
#' @param object An object inheriting from class
#'     \code{ral_rep}.
#' @param aggregation Character string: \code{"median"}
#'     (default), \code{"mean"}, or \code{"spectral"}.
#' @param ... Currently unused.
#'
#' @return Named vector (single fit) or matrix (multiple).
#'
#' @method coef ral_rep
#' @export
coef.ral_rep <- function(object,
                         aggregation = c("median", "mean",
                                         "spectral"),
                         ...) {
  aggregation <- match.arg(aggregation)
  agg <- aggregate_reps(object, aggregation = aggregation)
  cf <- agg$coefficients
  rownames(cf) <- object$coef_names
  colnames(cf) <- object$fit_labels
  if (ncol(cf) == 1) cf <- drop(cf)
  cf
}#COEF.RAL_REP

#' Variance-Covariance Matrix for RAL Rep Objects
#'
#' @param object An object inheriting from class
#'     \code{ral_rep}.
#' @param fit_idx Integer index of the fit. Defaults to 1.
#' @param aggregation Character string. Aggregation rule.
#' @param type Character. HC type. Default \code{"HC1"}.
#' @param ... Currently unused.
#'
#' @return A \eqn{p \times p}{p x p} variance-covariance
#'     matrix.
#'
#' @method vcov ral_rep
#' @importFrom stats vcov
#' @export
vcov.ral_rep <- function(object, fit_idx = 1,
                         aggregation = c("median", "mean", "spectral"),
                         type = "HC1", ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  aggregation <- match.arg(aggregation)
  agg <- aggregate_reps(object, aggregation = aggregation, type = type)
  V <- agg$vcov[, , fit_idx]
  p <- length(object$coef_names)
  if (!is.matrix(V)) V <- matrix(V, nrow = p, ncol = p)
  rownames(V) <- colnames(V) <- object$coef_names
  V
}#VCOV.RAL_REP

#' Confidence Intervals for RAL Rep Objects
#'
#' @param object An object inheriting from class
#'     \code{ral_rep}.
#' @param parm Parameter specification (names or indices).
#' @param level Confidence level. Default 0.95.
#' @param fit_idx Integer index of the fit. Defaults to 1.
#' @param aggregation Character string. Aggregation rule.
#' @param type Character. HC type. Default \code{"HC1"}.
#' @param uniform Logical. Uniform bands via multiplier
#'     bootstrap? Default \code{FALSE}.
#' @param bootstraps Integer. Bootstrap draws. Default 999.
#' @param ... Currently unused.
#'
#' @return A matrix with columns for lower and upper bounds.
#'     When \code{uniform = TRUE}, the attribute
#'     \code{"crit_val"} contains the aggregated critical
#'     value.
#'
#' @references
#' Chernozhukov V, Chetverikov D, Kato K (2013). "Gaussian
#' approximations and multiplier bootstrap for maxima of sums
#' of high-dimensional random vectors." Annals of Statistics,
#' 41(6), 2786-2819.
#'
#' @importFrom stats confint coef
#' @method confint ral_rep
#' @export
confint.ral_rep <- function(object, parm = NULL,
                            level = 0.95,
                            fit_idx = 1,
                            aggregation = c("median", "mean",
                                            "spectral"),
                            type = "HC1",
                            uniform = FALSE,
                            bootstraps = 999L, ...) {
  aggregation <- match.arg(aggregation)
  agg <- aggregate_reps(object, aggregation = aggregation, type = type)
  cf <- agg$coefficients[, fit_idx]
  se <- agg$se[, fit_idx]
  cf_names <- object$coef_names
  names(cf) <- cf_names
  names(se) <- cf_names

  # Select parameters
  if (is.null(parm)) {
    parm <- cf_names
  } else if (is.numeric(parm)) {
    parm <- cf_names[parm]
  } else {
    parm <- intersect(parm, cf_names)
    if (length(parm) == 0) {
      stop("None of the specified 'parm' were found ",
           "in the model coefficients.", call. = FALSE)
    }#IF
  }#IFELSE
  cf <- cf[parm]
  se <- se[parm]

  # Construct confidence intervals
  if (uniform) {
    R <- object$nresamples
    crit_vals <- vapply(seq_len(R), function(r) {
      ci_r <- confint(object$fits[[r]],
                       level = level,
                       fit_idx = fit_idx,
                       type = type,
                       uniform = TRUE,
                       bootstraps = bootstraps)
      attr(ci_r, "crit_val")
    }, numeric(1))
    agg_fn <- if (aggregation == "mean") mean else stats::median
    z <- agg_fn(crit_vals)
  } else {
    z <- stats::qnorm((1 + level) / 2)
  }#IFELSE
  ci <- cbind(cf - z * se, cf + z * se)
  pct <- c((1 - level) / 2, (1 + level) / 2) * 100
  colnames(ci) <- paste0(format(pct, digits = 3), " %")
  rownames(ci) <- parm
  attr(ci, "crit_val") <- z
  ci
}#CONFINT.RAL_REP

#' Summary for RAL Rep Objects
#'
#' Aggregates coefficient estimates and covariance matrices
#' across independent replications.
#'
#' @details
#' Let \eqn{\hat\theta_s} and \eqn{\hat\Sigma_s} denote
#' the coefficient vector and sandwich covariance matrix
#' from replication \eqn{s}.
#'
#' \strong{Coefficient aggregation.}
#' For \code{"mean"}:
#' \eqn{\tilde\theta = S^{-1} \sum_{s=1}^{S} \hat\theta_s}.
#' For \code{"median"} and \code{"spectral"}:
#' \eqn{\tilde\theta_j = \mathrm{median}_{s}(\hat\theta_{s,j})}.
#'
#' \strong{Covariance aggregation.}
#' Define the inflated per-replication covariance as
#' \deqn{V_s = \hat\Sigma_s + (\hat\theta_s - \tilde\theta)
#'     (\hat\theta_s - \tilde\theta)^\top.}
#' For \code{"mean"}:
#' \eqn{\tilde\Sigma = S^{-1} \sum_{s=1}^{S} V_s}.
#' For \code{"median"}:
#' \eqn{\tilde\Sigma_{s,ij} = \mathrm{median}_{s}(V_{s,ij})}.
#' For \code{"spectral"}:
#' solved via \pkg{CVXR}, guaranteeing PSD.
#'
#' @references
#' Chernozhukov V, Chetverikov D, Demirer M, Duflo E,
#'     Hansen C B, Newey W, Robins J (2018).
#'     "Double/debiased machine learning for treatment
#'     and structural parameters." The Econometrics
#'     Journal, 21(1), C1-C68.
#'
#' @param object An object inheriting from class
#'     \code{ral_rep}.
#' @param aggregation Character string. Aggregation rule.
#' @param type Character. HC type. Default \code{"HC1"}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"summary.ral_rep"}.
#'
#' @method summary ral_rep
#' @export
summary.ral_rep <- function(object,
                            aggregation = c("median", "mean", "spectral"),
                            type = "HC1", ...) {
  aggregation <- match.arg(aggregation)
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  agg <- aggregate_reps(object, aggregation = aggregation, type = type)

  p <- length(object$coef_names)
  nfit <- object$nfit
  fit_labels <- object$fit_labels
  if (is.null(fit_labels)) fit_labels <- paste0("fit", seq_len(nfit))

  inf_results <- array(0, dim = c(p, 4, nfit))
  for (j in seq_len(nfit)) {
    theta_j <- agg$coefficients[, j]
    se_j <- agg$se[, j]
    t_val <- theta_j / se_j
    p_val <- 2 * stats::pnorm(abs(t_val), lower.tail = FALSE)
    inf_results[, 1, j] <- theta_j
    inf_results[, 2, j] <- se_j
    inf_results[, 3, j] <- t_val
    inf_results[, 4, j] <- p_val
  }#FOR
  dimnames(inf_results) <- list(
    object$coef_names,
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
    fit_labels)

  result <- list(
    coefficients = inf_results,
    type = type,
    estimator_name = object$estimator_name,
    nobs = object$nobs,
    fit_labels = fit_labels,
    nresamples = object$nresamples,
    aggregation = aggregation)
  class(result) <- "summary.ral_rep"
  result
}#SUMMARY.RAL_REP

#' @rdname summary.ral_rep
#'
#' @param x An object of class \code{summary.ral_rep}.
#' @param digits Number of significant digits. Default 3.
#'
#' @method print summary.ral_rep
#' @export
print.summary.ral_rep <- function(x, digits = 3, ...) {
  cat("RAL estimation:", x$estimator_name, "\n")
  cat("Obs:", x$nobs,
      "  Resamples:", x$nresamples,
      "  Aggregation:", x$aggregation)
  if (!is.null(x$type) && x$type != "HC1") cat("  SE:", x$type)
  cat("\n\n")

  print_coef_tables(x$coefficients, fit_label = "Fit", digits = digits)

  invisible(x)
}#PRINT.SUMMARY.RAL_REP

#' Tidy a RAL Rep Object
#'
#' @param x An object inheriting from class \code{ral_rep}.
#' @param fit_idx Integer index of the fit. Defaults to 1.
#'     Set to \code{NULL} for all fits.
#' @param aggregation Character string. Aggregation rule.
#' @param type Character. HC type. Default \code{"HC1"}.
#' @param conf.int Logical. Include CIs? Default
#'     \code{FALSE}.
#' @param conf.level Confidence level. Default 0.95.
#' @param uniform Logical. Uniform CIs? Default
#'     \code{FALSE}.
#' @param bootstraps Integer. Bootstrap draws. Default 999.
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns \code{term},
#'     \code{estimate}, \code{std.error}, \code{statistic},
#'     \code{p.value}, \code{fit_label}, and
#'     \code{aggregation}.
#'
#' @method tidy ral_rep
#' @export
tidy.ral_rep <- function(x, fit_idx = 1,
                         aggregation = c("median", "mean", "spectral"),
                         type = "HC1",
                         conf.int = FALSE,
                         conf.level = 0.95,
                         uniform = FALSE,
                         bootstraps = 999L, ...) {
  aggregation <- match.arg(aggregation)
  s <- summary(x, aggregation = aggregation, type = type)
  inf <- s$coefficients
  nfit <- dim(inf)[3]
  p <- dim(inf)[1]

  if (is.null(fit_idx)) {
    j_seq <- seq_len(nfit)
  } else {
    j_seq <- fit_idx
  }#IFELSE

  rows <- list()
  for (j in j_seq) {
    for (k in seq_len(p)) {
      row <- data.frame(
        term = dimnames(inf)[[1]][k],
        estimate = inf[k, 1, j],
        std.error = inf[k, 2, j],
        statistic = inf[k, 3, j],
        p.value = inf[k, 4, j],
        fit_label = dimnames(inf)[[3]][j],
        aggregation = aggregation,
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1]] <- row
    }#FOR
  }#FOR

  res <- do.call(rbind, rows)
  if (conf.int) {
    ci_list <- lapply(j_seq, function(j) {
      confint(x, level = conf.level, fit_idx = j,
              aggregation = aggregation, type = type,
              uniform = uniform, bootstraps = bootstraps)
    })
    ci_mat <- do.call(rbind, ci_list)
    res$conf.low <- as.numeric(ci_mat[, 1])
    res$conf.high <- as.numeric(ci_mat[, 2])
  }#IF
  res
}#TIDY.RAL_REP

#' Glance at a RAL Rep Object
#'
#' @param x An object inheriting from class \code{ral_rep}.
#' @param ... Currently unused.
#'
#' @return A one-row \code{data.frame} with columns
#'     \code{nobs}, \code{nresamples}, and
#'     \code{estimator_name}.
#'
#' @method glance ral_rep
#' @export
glance.ral_rep <- function(x, ...) {
  data.frame(
    nobs = x$nobs,
    nresamples = x$nresamples,
    estimator_name = if (is.null(x$estimator_name)) {
      class(x)[1]
    } else {
      x$estimator_name
    },
    stringsAsFactors = FALSE
  )
}#GLANCE.RAL_REP

# Plot methods ================================================================

#' @rdname plot.ral
#' @method plot ral_rep
#' @export
plot.ral_rep <- function(x, parm = NULL, level = 0.95,
                         uniform = TRUE,
                         type = "HC1",
                         xlab = NULL, ylab = NULL,
                         main = NULL,
                         col = "black", pch = 19, lwd = 1.5,
                         ...) {
  # Compute CIs using ral_rep method (aggregates across reps)
  ci <- confint(x, parm = parm, level = level, type = type, uniform = uniform)
  cf <- coef(x)
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
}#PLOT.RAL_REP
