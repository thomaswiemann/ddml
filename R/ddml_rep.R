# ddml_rep: DML-specific replicated inference ===================================
#
# Class hierarchy:
#   ddml_rep > ral_rep — adds DML-specific fields and display
#
# Inference methods (coef, vcov, confint, tidy, glance) are
# inherited from ral_rep in ral_rep.R.

# Exported functions ===========================================================

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
#' @return An object of class \code{c("ddml_rep", "ral_rep")}
#'     with fields:
#' \describe{
#'   \item{fits}{List of \code{ddml} objects.}
#'   \item{nresamples}{Number of resamples.}
#'   \item{model_type}{Primary class of the fits.}
#'   \item{coef_names}{Coefficient names.}
#'   \item{ensemble_type}{Ensemble types.}
#'   \item{nobs}{Number of observations.}
#'   \item{sample_folds}{Number of cross-fitting folds.}
#'   \item{shortstack}{Logical, whether short-stacking was used.}
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
  # Input validation
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
    if (!identical(fits[[i]]$coef_names, ref$coef_names)) {
      stop("Fit ", i, " has different 'coef_names' ",
           "than fit 1.", call. = FALSE)
    }#IF
    if (!identical(fits[[i]]$ensemble_type,
                   ref$ensemble_type)) {
      stop("Fit ", i, " has different 'ensemble_type' ",
           "than fit 1.", call. = FALSE)
    }#IF
    if (!identical(fits[[i]]$nobs, ref$nobs)) {
      stop("Fit ", i, " has different 'nobs' than fit 1.",
           call. = FALSE)
    }#IF
  }#FOR

  # Build via ral_rep, then layer DML fields
  ens_type <- ref$ensemble_type
  if (is.null(ens_type)) ens_type <- "single base learner"

  obj <- ral_rep(fits, subclass = "ddml_rep")
  obj$model_type <- primary[1]
  obj$ensemble_type <- ens_type
  obj$fit_labels <- ens_type
  obj$sample_folds <- ref$sample_folds
  obj$shortstack <- ref$shortstack
  obj
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
  # Suppress inner start/finish messages
  if (is.null(dots$messages)) {
    dots$messages <- list(start = "", finish = "")
  } else {
    dots$messages$start <- ""
    dots$messages$finish <- ""
  }#IFELSE
  fits <- vector("list", resamples)
  for (r in seq_len(resamples)) {
    if (!silent) message("[Resample ", r, "/", resamples, "]")
    fits[[r]] <- do.call(fn, dots)
  }#FOR
  ddml_rep(fits)
}#DDML_REPLICATE

# DML-specific S3 methods =====================================================
# Inference methods (coef, nobs, vcov, confint, [[, length)
# are inherited from ral_rep in ral_rep.R.

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

#' Summary for ddml_rep Objects
#'
#' DML-specific summary override. Adds ensemble type labels,
#' folds, shortstack status to the base \code{ral_rep}
#' summary.
#'
#' @details
#' See \code{\link{summary.ral_rep}} for the aggregation
#' formulas.
#'
#' @references
#' Chernozhukov V, Chetverikov D, Demirer M, Duflo E,
#'     Hansen C B, Newey W, Robins J (2018).
#'     "Double/debiased machine learning for treatment
#'     and structural parameters." The Econometrics
#'     Journal, 21(1), C1-C68.
#'
#' @param object A \code{ddml_rep} object.
#' @param aggregation Character string: \code{"median"}
#'     (default), \code{"mean"}, or \code{"spectral"}.
#' @param type Character. HC type. Default \code{"HC1"}.
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
                             aggregation = c("median", "mean",
                                             "spectral"),
                             type = "HC1", ...) {
  aggregation <- match.arg(aggregation)
  type <- match.arg(type, c("HC0", "HC1", "HC3"))

  # DML-specific fit labels
  object$fit_labels <- object$ensemble_type

  # Delegate table computation to ral_rep
  result <- summary.ral_rep(object,
                            aggregation = aggregation,
                            type = type, ...)

  # Attach DML-specific fields
  result$model_type    <- object$model_type
  result$sample_folds  <- object$sample_folds
  result$shortstack    <- object$shortstack
  result$ensemble_type <- object$ensemble_type

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

  print_coef_tables(x$coefficients,
                    fit_label = "Ensemble type",
                    digits = digits)

  invisible(x)
}#PRINT.SUMMARY.DDML_REP

#' Tidy a ddml_rep Object
#'
#' DML-specific tidy method. Adds \code{ensemble_type} and
#' \code{aggregation} columns. Delegates to
#' \code{tidy.ral_rep} for the base table computation.
#'
#' @param x A \code{ddml_rep} object.
#' @param ensemble_idx Integer index of the ensemble type
#'     to report. Defaults to 1. Set to \code{NULL} for
#'     all ensemble types.
#' @param aggregation Character string. Aggregation method.
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
#'     \code{p.value}, \code{ensemble_type}, and
#'     \code{aggregation}.
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
                          aggregation = c("median", "mean",
                                          "spectral"),
                          type = "HC1",
                          conf.int = FALSE,
                          conf.level = 0.95,
                          uniform = FALSE,
                          bootstraps = 999L, ...) {
  aggregation <- match.arg(aggregation)
  res <- tidy.ral_rep(x, fit_idx = ensemble_idx,
                      aggregation = aggregation,
                      type = type,
                      conf.int = conf.int,
                      conf.level = conf.level,
                      uniform = uniform,
                      bootstraps = bootstraps)
  # Rename fit_label -> ensemble_type for DML compatibility
  names(res)[names(res) == "fit_label"] <- "ensemble_type"
  res
}#TIDY.DDML_REP

#' Glance at a ddml_rep Object
#'
#' DML-specific glance method. Includes DML fields.
#'
#' @param x A \code{ddml_rep} object.
#' @param ... Currently unused.
#'
#' @return A one-row \code{data.frame}.
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
    shortstack = x$shortstack,
    ensemble_type = paste(x$ensemble_type,
                          collapse = ", "),
    model_type = x$model_type,
    estimator_name = x$estimator_name,
    nresamples = x$nresamples,
    stringsAsFactors = FALSE
  )
}#GLANCE.DDML_REP

# List conversion =============================================================

#' Split a ddml_rep Object by Ensemble Type
#'
#' Returns a named list of single-ensemble
#'     \code{ddml_rep} objects.
#'
#' @param x A \code{ddml_rep} object.
#' @param ... Currently unused.
#'
#' @return A named list of \code{ddml_rep} objects.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' reps = ddml_replicate(ddml_plm, y = y, D = D, X = X,
#'                       learners = list(what = ols),
#'                       resamples = 3,
#'                       sample_folds = 2,
#'                       silent = TRUE)
#' as.list(reps)
#' }
#'
#' @method as.list ddml_rep
#' @export
as.list.ddml_rep <- function(x, ...) {
  nfit <- x$nfit
  labels <- x$ensemble_type
  if (is.null(labels)) labels <- x$fit_labels
  if (is.null(labels)) labels <- paste0("fit", seq_len(nfit))

  out <- vector("list", nfit)
  names(out) <- labels
  for (j in seq_len(nfit)) {
    fits_j <- lapply(x$fits, function(f) as.list(f)[[j]])
    obj <- ral_rep(fits_j, subclass = "ddml_rep")
    obj$model_type <- x$model_type
    obj$ensemble_type <- labels[j]
    obj$fit_labels <- labels[j]
    obj$sample_folds <- x$sample_folds
    obj$shortstack <- x$shortstack
    out[[j]] <- obj
  }#FOR
  out
}#AS.LIST.DDML_REP
