# Internal helpers -----------------------------------------------

# Build inference results array from aggregated coef/SE.
build_inf_from_agg <- function(agg, coef_names,
                               ensemble_type) {
  p <- nrow(agg$coefficients)
  nensb <- ncol(agg$coefficients)
  inf_results <- array(0, dim = c(p, 4, nensb))
  for (j in seq_len(nensb)) {
    theta_j <- agg$coefficients[, j]
    se_j <- agg$se[, j]
    t_val <- theta_j / se_j
    p_val <- 2 * stats::pnorm(abs(t_val),
                               lower.tail = FALSE)
    inf_results[, 1, j] <- theta_j
    inf_results[, 2, j] <- se_j
    inf_results[, 3, j] <- t_val
    inf_results[, 4, j] <- p_val
  }#FOR
  dimnames(inf_results) <- list(
    coef_names,
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
    ensemble_type)
  inf_results
}#BUILD_INF_FROM_AGG

# Spectral-norm median of PSD matrices via SDP (CVXR).
#
# Finds V* = argmin_{V >= 0} sum_r ||V - V_r||_op where
# ||.||_op is the spectral norm (largest singular value).
# For p = 1 this reduces to the standard scalar median.
spectral_median_psd <- function(matrices) {
  p <- nrow(matrices[[1]])
  R <- length(matrices)

  if (p == 1) {
    vals <- vapply(matrices, function(m) m[1, 1],
                   numeric(1))
    return(matrix(stats::median(vals), 1, 1))
  }#IF

  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' is required for spectral ",
         "aggregation. Install it with:\n",
         "  install.packages('CVXR')",
         call. = FALSE)
  }#IF

  V <- CVXR::Variable(c(p, p), PSD = TRUE)
  obj <- 0
  for (r in seq_len(R)) {
    obj <- obj + CVXR::norm(V - matrices[[r]], "2")
  }#FOR

  result <- CVXR::solve(CVXR::Problem(CVXR::Minimize(obj)))

  if (result$status != "optimal") {
    warning("SDP solver returned status '",
            result$status,
            "'; falling back to element-wise median.",
            call. = FALSE)
    arr <- array(
      unlist(lapply(matrices, as.vector)),
      dim = c(p, p, R))
    return(apply(arr, c(1, 2), stats::median))
  }#IF

  V_sol <- result$getValue(V)
  (V_sol + t(V_sol)) / 2
}#SPECTRAL_MEDIAN_PSD

# Core aggregation workhorse.
#
# For each replication r the inflated covariance is
#   V_r = Sigma_r + (theta_r - theta_tilde)(theta_r - theta_tilde)'
# where theta_tilde is the aggregated coefficient vector.
# The three aggregation rules then differ only in how they
# summarise {V_1, ..., V_R} into a single matrix.
aggregate_reps <- function(object, aggregation = "median",
                           type = "HC1") {
  aggregation <- match.arg(aggregation,
                           c("median", "mean",
                             "spectral"))
  R <- object$nresamples
  nensb <- length(object$ensemble_type)
  p <- length(object$coef_names)

  # == Collect per-replication estimates ==
  coef_array <- array(0, dim = c(p, nensb, R))
  vcov_array <- array(0, dim = c(p, p, nensb, R))
  for (r in seq_len(R)) {
    fit <- object$fits[[r]]
    coef_array[, , r] <- fit$coefficients
    for (j in seq_len(nensb)) {
      vcov_array[, , j, r] <- stats::vcov(fit,
        ensemble_idx = j, type = type)
    }#FOR
  }#FOR

  # == Aggregate coefficients ==
  if (aggregation == "mean") {
    agg_coef <- apply(coef_array, c(1, 2), mean)
  } else {
    agg_coef <- apply(coef_array, c(1, 2), stats::median)
  }#IFELSE

  # == Inflate & aggregate covariance ==
  agg_vcov <- array(0, dim = c(p, p, nensb))

  for (j in seq_len(nensb)) {
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
      V_arr <- array(
        unlist(lapply(V_list, as.vector)),
        dim = c(p, p, R))
      agg_vcov[, , j] <- apply(V_arr, c(1, 2),
                                stats::median)
    }#IFELSE
  }#FOR

  # == Standard errors ==
  agg_se <- matrix(0, nrow = p, ncol = nensb)
  for (j in seq_len(nensb)) {
    V_j <- matrix(agg_vcov[, , j, drop = FALSE],
                  nrow = p, ncol = p)
    agg_se[, j] <- sqrt(diag(V_j))
  }#FOR

  list(coefficients = agg_coef, se = agg_se,
       vcov = agg_vcov, coef_array = coef_array,
       vcov_array = vcov_array)
}#AGGREGATE_REPS

# Exported functions ---------------------------------------------

#' Construct a Multi-Resample DDML Object
#'
#' Validates a list of \code{ddml} fits and stamps class
#' \code{"ddml_rep"} for multi-resample aggregation.
#'
#' @param fits A list of at least 2 objects inheriting from
#'     class \code{"ddml"}. All fits must share the same
#'     primary class, coefficient names, ensemble type, and
#'     number of observations.
#'
#' @return An object of class \code{"ddml_rep"} with fields:
#' \describe{
#'   \item{fits}{List of \code{ddml} objects.}
#'   \item{nresamples}{Number of resamples.}
#'   \item{model_type}{Primary class of the fits.}
#'   \item{coef_names}{Coefficient names.}
#'   \item{ensemble_type}{Ensemble types.}
#'   \item{nobs}{Number of observations.}
#'   \item{sample_folds}{Number of cross-fitting folds.}
#'   \item{shortstack}{Logical, whether short-stacking was
#'       used.}
#' }
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' fits = lapply(1:3, function(r) {
#'   ddml_plm(y, D, X,
#'            learners = list(what = ols),
#'            sample_folds = 2, silent = TRUE)
#' })
#' reps = ddml_rep(fits)
#' summary(reps)
#' }
#'
#' @family ddml replication
#' @seealso [ddml_replicate()]
#' @export
ddml_rep <- function(fits) {
  if (!is.list(fits) || length(fits) < 2) {
    stop("'fits' must be a list of at least 2 ddml objects.",
         call. = FALSE)
  }#IF

  for (i in seq_along(fits)) {
    if (!inherits(fits[[i]], "ddml")) {
      stop("Element ", i,
           " does not inherit from class 'ddml'.",
           call. = FALSE)
    }#IF
  }#FOR

  primary <- vapply(fits, function(f) class(f)[1],
                    character(1))
  if (length(unique(primary)) != 1) {
    stop("All fits must have the same primary class. ",
         "Found: ",
         paste(unique(primary), collapse = ", "),
         call. = FALSE)
  }#IF

  ref <- fits[[1]]
  for (i in seq_along(fits)[-1]) {
    if (!identical(fits[[i]]$coef_names,
                   ref$coef_names)) {
      stop("Fit ", i,
           " has different 'coef_names' than fit 1.",
           call. = FALSE)
    }#IF
    if (!identical(fits[[i]]$ensemble_type,
                   ref$ensemble_type)) {
      stop("Fit ", i,
           " has different 'ensemble_type' than fit 1.",
           call. = FALSE)
    }#IF
    if (!identical(fits[[i]]$nobs, ref$nobs)) {
      stop("Fit ", i,
           " has different 'nobs' than fit 1.",
           call. = FALSE)
    }#IF
  }#FOR

  ens_type <- ref$ensemble_type
  if (is.null(ens_type)) {
    ens_type <- "single base learner"
  }#IF

  structure(
    list(
      fits          = fits,
      nresamples    = length(fits),
      model_type    = primary[1],
      coef_names    = ref$coef_names,
      estimator_name= if (!is.null(ref$estimator_name)) ref$estimator_name else primary[1],
      ensemble_type = ens_type,
      nobs          = ref$nobs,
      sample_folds  = ref$sample_folds,
      shortstack    = ref$shortstack
    ),
    class = "ddml_rep"
  )
}#DDML_REP

#' Replicate a DDML Estimator Across Multiple Resamples
#'
#' Convenience wrapper that calls a \code{ddml_*} estimator
#' function multiple times with independent sample splits
#' and returns a \code{ddml_rep} object for aggregated
#' inference.
#'
#' @param fn A \code{ddml_*} estimator function
#'     (e.g., \code{ddml_plm}).
#' @param ... Arguments passed to \code{fn}.
#' @param resamples Integer number of independent resamples.
#'     Must be >= 2. Default 5.
#' @param silent Logical. If \code{TRUE}, suppresses all
#'     output at both the resample level and within each
#'     estimator call. Default \code{FALSE}.
#'
#' @return An object of class \code{"ddml_rep"}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
#'                       learners = list(what = ols),
#'                       sample_folds = 2,
#'                       resamples = 3, silent = TRUE)
#' summary(reps)
#' }
#'
#' @family ddml replication
#' @seealso [ddml_rep()]
#' @export
ddml_replicate <- function(fn, ..., resamples = 5,
                           silent = FALSE) {
  dots <- list(...)
  dots$silent <- silent
  # Suppress inner start/finish messages to avoid console spam
  if (is.null(dots$messages)) {
    dots$messages <- list(start = "", finish = "")
  } else {
    dots$messages$start <- ""
    dots$messages$finish <- ""
  }
  fits <- vector("list", resamples)
  for (r in seq_len(resamples)) {
    if (!silent) {
      message("[Resample ", r, "/", resamples, "]")
    }#IF
    fits[[r]] <- do.call(fn, dots)
  }#FOR
  ddml_rep(fits)
}#DDML_REPLICATE

# S3 methods -----------------------------------------------------

#' Extract a Single Fit from a ddml_rep Object
#'
#' @param x A \code{ddml_rep} object.
#' @param i Integer index of the fit to extract.
#'
#' @return The \code{i}-th \code{ddml} fit.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
#'                       learners = list(what = ols),
#'                       sample_folds = 2,
#'                       resamples = 3, silent = TRUE)
#' coef(reps[[1]])
#' }
#'
#' @export
#' @method [[ ddml_rep
`[[.ddml_rep` <- function(x, i) {
  x$fits[[i]]
}#[[.DDML_REP

#' Number of Resamples in a ddml_rep Object
#'
#' @param x A \code{ddml_rep} object.
#'
#' @return Integer number of resamples.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
#'                       learners = list(what = ols),
#'                       sample_folds = 2,
#'                       resamples = 3, silent = TRUE)
#' length(reps)
#' }
#'
#' @export
#' @method length ddml_rep
length.ddml_rep <- function(x) {
  x$nresamples
}#LENGTH.DDML_REP

#' Extract Number of Observations from a ddml_rep Object
#'
#' @param object A \code{ddml_rep} object.
#' @param ... Currently unused.
#'
#' @return An integer specifying the number of observations.
#'
#' @importFrom stats nobs
#' @export
#' @method nobs ddml_rep
nobs.ddml_rep <- function(object, ...) {
  object$nobs
}#NOBS.DDML_REP

#' @rdname ddml_rep
#'
#' @param x A \code{ddml_rep} object.
#' @param ... Currently unused.
#'
#' @export
#' @method print ddml_rep
print.ddml_rep <- function(x, ...) {
  cat("DDML replicated fits:", x$estimator_name, "\n")
  cat("  Resamples:", x$nresamples,
      "  Obs:", x$nobs,
      "  Folds:", x$sample_folds, "\n\n")
  cat("Use summary() for aggregated inference.\n")
  cat("Use x[[i]] to access individual fits.\n")
  invisible(x)
}#PRINT.DDML_REP

#' Extract Aggregated Coefficients from a ddml_rep Object
#'
#' @param object A \code{ddml_rep} object.
#' @param aggregation Character string: \code{"median"}
#'     (default), \code{"mean"}, or \code{"spectral"}.
#'     See \code{\link{summary.ddml_rep}} for details.
#' @param ... Additional arguments. 
#'
#' @return Named vector (single ensemble) or matrix
#'     (multiple ensembles).
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
#'                       learners = list(what = ols),
#'                       sample_folds = 2,
#'                       resamples = 3, silent = TRUE)
#' coef(reps)
#' coef(reps, aggregation = "mean")
#' }
#'
#' @seealso \code{\link{summary.ddml_rep}} for the
#'     aggregation equations.
#'
#' @export
#' @method coef ddml_rep
coef.ddml_rep <- function(object,
                          aggregation = c("median", "mean", "spectral"),
                          ...) {
  aggregation <- match.arg(aggregation)
  agg <- aggregate_reps(object,
                        aggregation = aggregation)
  cf <- agg$coefficients
  rownames(cf) <- object$coef_names
  colnames(cf) <- object$ensemble_type
  if (ncol(cf) == 1) cf <- drop(cf)
  cf
}#COEF.DDML_REP

#' Variance-Covariance Matrix for ddml_rep Objects
#'
#' Returns a variance-covariance matrix from cross-resample
#' aggregation.
#'
#' @param object A \code{ddml_rep} object.
#' @param ensemble_idx Integer index of the ensemble type.
#'     Defaults to 1.
#' @inheritParams coef.ddml_rep
#' @inheritParams vcov.ddml
#' @param ... Currently unused.
#'
#' @return A p x p variance-covariance matrix.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
#'                       learners = list(what = ols),
#'                       sample_folds = 2,
#'                       resamples = 3, silent = TRUE)
#' vcov(reps)
#' }
#'
#' @seealso \code{\link{summary.ddml_rep}} for the
#'     aggregation equations.
#'
#' @export
#' @method vcov ddml_rep
vcov.ddml_rep <- function(object, ensemble_idx = 1,
                          aggregation = c("median", "mean", "spectral"),
                          type = "HC1", ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  aggregation <- match.arg(aggregation)
  agg <- aggregate_reps(object,
                        aggregation = aggregation,
                        type = type)
  V <- agg$vcov[, , ensemble_idx]
  if (!is.matrix(V)) {
    V <- matrix(V, nrow = length(object$coef_names),
                ncol = length(object$coef_names))
  }#IF
  rownames(V) <- colnames(V) <- object$coef_names
  V
}#VCOV.DDML_REP

#' Confidence Intervals for ddml_rep Objects
#'
#' @param object A \code{ddml_rep} object.
#' @inheritParams confint.ddml
#' @param ensemble_idx Integer index of the ensemble type.
#'     Defaults to 1.
#' @inheritParams coef.ddml_rep
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
#' reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
#'                       learners = list(what = ols),
#'                       sample_folds = 2,
#'                       resamples = 3, silent = TRUE)
#' confint(reps)
#' confint(reps, level = 0.90)
#' }
#'
#' @seealso \code{\link{summary.ddml_rep}} for the
#'     aggregation equations.
#'
#' @export
#' @method confint ddml_rep
confint.ddml_rep <- function(object, parm,
                             level = 0.95,
                             ensemble_idx = 1,
                             aggregation = c("median", "mean", "spectral"),
                             type = "HC1", ...) {
  aggregation <- match.arg(aggregation)
  agg <- aggregate_reps(object,
                        aggregation = aggregation,
                        type = type)
  cf <- agg$coefficients[, ensemble_idx]
  se <- agg$se[, ensemble_idx]
  cf_names <- object$coef_names
  names(cf) <- cf_names
  names(se) <- cf_names

  if (missing(parm)) {
    parm <- cf_names
  } else if (is.numeric(parm)) {
    parm <- cf_names[parm]
  } else {
    parm <- intersect(parm, cf_names)
    if (length(parm) == 0) {
      stop("None of the specified 'parm' were found ",
           "in the model coefficients.",
           call. = FALSE)
    }#IF
  }#IFELSE

  cf <- cf[parm]
  se <- se[parm]
  z <- stats::qnorm((1 + level) / 2)
  ci <- cbind(cf - z * se, cf + z * se)
  pct <- c((1 - level) / 2, (1 + level) / 2) * 100
  colnames(ci) <- paste0(format(pct, digits = 3), " %")
  rownames(ci) <- parm
  ci
}#CONFINT.DDML_REP

#' Summary for ddml_rep Objects
#'
#' Aggregates coefficient estimates and covariance matrices
#' across \eqn{R} independent sample-splitting replications.
#'
#' @details
#' Let \eqn{\hat\theta_r} and \eqn{\hat\Sigma_r} denote
#' the coefficient vector and sandwich covariance matrix
#' from replication \eqn{r}.
#'
#' \strong{Coefficient aggregation.}
#' For \code{"mean"}:
#' \eqn{\tilde\theta = R^{-1} \sum_{r=1}^{R} \hat\theta_r}.
#' For \code{"median"} and \code{"spectral"}:
#' \eqn{\tilde\theta_j = \mathrm{median}_{r}(\hat\theta_{r,j})}.
#'
#' \strong{Covariance aggregation.}
#' Define the inflated per-replication covariance as
#' \deqn{V_r = \hat\Sigma_r + (\hat\theta_r - \tilde\theta)
#'     (\hat\theta_r - \tilde\theta)^\top}.
#' For \code{"mean"}:
#' \eqn{\tilde\Sigma = R^{-1} \sum_{r=1}^{R} V_r}.
#' For \code{"median"}:
#' \eqn{\tilde\Sigma_{ij} = \mathrm{median}_{r}(V_{r,ij})}.
#' For \code{"spectral"}:
#' \eqn{\tilde\Sigma = \arg\min_{\Sigma \succeq 0} \sum_{r} \|\Sigma - V_r\|_2},
#' solved via \pkg{CVXR}. Guarantees positive semi-definiteness of \eqn{\tilde\Sigma}.
#'
#' @references
#' Chernozhukov V, Chetverikov D, Demirer M, Duflo E,
#'     Hansen C B, Newey W, Robins J (2018).
#'     "Double/debiased machine learning for treatment
#'     and structural parameters." The Econometrics
#'     Journal, 21(1), C1-C68.
#'
#' @param object A \code{ddml_rep} object.
#' @inheritParams coef.ddml_rep
#' @param type Character string specifying the
#'     variance-covariance estimator. One of \code{"HC1"}
#'     (default), \code{"HC0"}, or \code{"HC3"}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"summary.ddml_rep"}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
#'                       learners = list(what = ols),
#'                       sample_folds = 2,
#'                       resamples = 3, silent = TRUE)
#' summary(reps)
#' summary(reps, aggregation = "mean")
#' }
#'
#' @export
#' @method summary ddml_rep
summary.ddml_rep <- function(object,
                             aggregation = c("median", "mean", "spectral"),
                             type = "HC1", ...) {
  aggregation <- match.arg(aggregation)
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  agg <- aggregate_reps(object,
                        aggregation = aggregation,
                        type = type)
  inf_results <- build_inf_from_agg(
    agg, object$coef_names, object$ensemble_type)
  result <- list(
    coefficients  = inf_results,
    type          = type,
    model_type    = object$model_type,
    estimator_name= object$estimator_name,
    nobs          = object$nobs,
    sample_folds  = object$sample_folds,
    shortstack    = object$shortstack,
    ensemble_type = object$ensemble_type,
    nresamples    = object$nresamples,
    aggregation   = aggregation
  )
  class(result) <- "summary.ddml_rep"
  result
}#SUMMARY.DDML_REP

#' @rdname summary.ddml_rep
#'
#' @param x An object of class \code{summary.ddml_rep}.
#' @param digits Number of significant digits. Default 3.
#'
#' @export
#' @method print summary.ddml_rep
print.summary.ddml_rep <- function(x, digits = 3, ...) {
  cat("DDML estimation:", x$estimator_name, "\n")
  cat("Obs:", x$nobs,
      "  Folds:", x$sample_folds,
      "  Resamples:", x$nresamples,
      "  Aggregation:", x$aggregation)
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
}#PRINT.SUMMARY.DDML_REP

#' Tidy a ddml_rep Object
#'
#' Extracts aggregated coefficient estimates, standard
#' errors, test statistics, and p-values from a
#' \code{ddml_rep} object in a format compatible with
#' \pkg{modelsummary} and the \pkg{broom} ecosystem.
#'
#' @param x A \code{ddml_rep} object.
#' @param ensemble_idx Integer index of the ensemble type
#'     to report. Defaults to 1. Set to \code{NULL} to
#'     return results for all ensemble types.
#' @inheritParams coef.ddml_rep
#' @param type Character string specifying the
#'     variance-covariance estimator. One of \code{"HC1"}
#'     (default), \code{"HC0"}, or \code{"HC3"}.
#' @param conf.int Logical. Include confidence interval
#'     columns? Default \code{FALSE}.
#' @param conf.level Confidence level for intervals.
#'     Default 0.95.
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns \code{term},
#'     \code{estimate}, \code{std.error}, \code{statistic},
#'     \code{p.value}, \code{ensemble_type}, and
#'     \code{aggregation}. If \code{conf.int = TRUE},
#'     also \code{conf.low} and \code{conf.high}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
#'                       learners = list(what = ols),
#'                       sample_folds = 2,
#'                       resamples = 3, silent = TRUE)
#' tidy(reps)
#' tidy(reps, conf.int = TRUE)
#' }
#'
#' @seealso \code{\link{summary.ddml_rep}} for the
#'     aggregation equations.
#'
#' @export
#' @method tidy ddml_rep
tidy.ddml_rep <- function(x, ensemble_idx = 1,
                          aggregation = c("median", "mean", "spectral"),
                          type = "HC1",
                          conf.int = FALSE,
                          conf.level = 0.95, ...) {
  aggregation <- match.arg(aggregation)
  s <- summary(x, aggregation = aggregation,
               type = type)
  inf <- s$coefficients
  nensb <- dim(inf)[3]
  p <- dim(inf)[1]

  if (is.null(ensemble_idx)) {
    j_seq <- seq_len(nensb)
  } else {
    j_seq <- ensemble_idx
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
        ensemble_type = dimnames(inf)[[3]][j],
        aggregation = aggregation,
        stringsAsFactors = FALSE
      )
      if (conf.int) {
        z <- stats::qnorm((1 + conf.level) / 2)
        row$conf.low <-
          inf[k, 1, j] - z * inf[k, 2, j]
        row$conf.high <-
          inf[k, 1, j] + z * inf[k, 2, j]
      }#IF
      rows[[length(rows) + 1]] <- row
    }#FOR
  }#FOR
  do.call(rbind, rows)
}#TIDY.DDML_REP

#' Glance at a ddml_rep Object
#'
#' Returns a one-row summary of model-level statistics,
#' compatible with \pkg{modelsummary} and the \pkg{broom}
#' ecosystem.
#'
#' @param x A \code{ddml_rep} object.
#' @param ... Currently unused.
#'
#' @return A one-row \code{data.frame} with columns
#'     \code{nobs}, \code{sample_folds}, \code{shortstack},
#'     \code{ensemble_type}, \code{model_type}, and
#'     \code{nresamples}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
#'                       learners = list(what = ols),
#'                       sample_folds = 2,
#'                       resamples = 3, silent = TRUE)
#' glance(reps)
#' }
#'
#' @export
#' @method glance ddml_rep
glance.ddml_rep <- function(x, ...) {
  data.frame(
    nobs = x$nobs,
    sample_folds = x$sample_folds,
    shortstack = if (is.null(x$shortstack)) {
      FALSE
    } else {
      x$shortstack
    },
    ensemble_type = paste(x$ensemble_type,
                          collapse = ", "),
    model_type = x$model_type,
    estimator_name = x$estimator_name,
    nresamples = x$nresamples,
    stringsAsFactors = FALSE
  )
}#GLANCE.DDML_REP
