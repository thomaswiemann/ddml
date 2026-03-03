# Internal helpers -----------------------------------------------

# Extract SEs from a single ddml fit.
extract_se <- function(object, type = "HC1") {
  nensb <- length(object$scores)
  p <- length(object$coef_names)
  se_mat <- matrix(0, p, nensb)
  for (j in seq_len(nensb)) {
    V <- compute_ddml_variance(
      object$scores[[j]], object$J[[j]],
      object$cluster_variable, type = type)
    se_mat[, j] <- sqrt(diag(V))
  }#FOR
  se_mat
}#EXTRACT_SE

# Ensure coefficients are always a p x nensb matrix.
normalize_coef_matrix <- function(coefficients, p, nensb) {
  if (is.matrix(coefficients)) {
    return(coefficients)
  }#IF
  if (p == 1) {
    matrix(coefficients, nrow = 1, ncol = nensb)
  } else {
    matrix(coefficients, nrow = p, ncol = 1)
  }#IFELSE
}#NORMALIZE_COEF_MATRIX

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
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    ensemble_type)
  inf_results
}#BUILD_INF_FROM_AGG

# Core aggregation workhorse (Ahrens et al., 2024, Remark 2).
aggregate_reps <- function(object, aggregation = "median",
                           type = "HC1") {
  aggregation <- match.arg(aggregation,
                           c("median", "mean"))
  R <- object$nresamples
  nensb <- length(object$ensemble_type)
  p <- length(object$coef_names)

  coef_array <- array(0, dim = c(p, nensb, R))
  se_array   <- array(0, dim = c(p, nensb, R))
  for (r in seq_len(R)) {
    fit <- object$fits[[r]]
    coef_array[, , r] <- normalize_coef_matrix(
      fit$coefficients, p, nensb)
    se_array[, , r] <- extract_se(fit, type = type)
  }#FOR

  if (aggregation == "median") {
    agg_coef <- apply(coef_array, c(1, 2),
                      stats::median)
    var_total <- se_array^2 +
      sweep(coef_array, c(1, 2), agg_coef)^2
    agg_se <- sqrt(apply(var_total, c(1, 2),
                         stats::median))
  } else {
    agg_coef <- apply(coef_array, c(1, 2), mean)
    var_total <- se_array^2 +
      sweep(coef_array, c(1, 2), agg_coef)^2
    agg_se <- sqrt(apply(var_total, c(1, 2),
      function(x) length(x) / sum(1 / x)))
  }#IFELSE

  list(coefficients = agg_coef, se = agg_se,
       coef_array = coef_array, se_array = se_array)
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
#' @family ddml
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
#' @family ddml
#' @seealso [ddml_rep()]
#' @export
ddml_replicate <- function(fn, ..., resamples = 5,
                           silent = FALSE) {
  dots <- list(...)
  dots$silent <- silent
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

#' Print a ddml_rep Object
#'
#' Displays a brief overview of a \code{ddml_rep} object.
#'
#' @param x A \code{ddml_rep} object.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
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
#' reps
#' }
#'
#' @export
#' @method print ddml_rep
print.ddml_rep <- function(x, ...) {
  type_labels <- c(
    ddml_plm  = "Partially Linear Model",
    ddml_pliv = "Partially Linear IV Model",
    ddml_fpliv =
      "Flexible Partially Linear IV Model",
    ddml_ate  = "Average Treatment Effect",
    ddml_att  =
      "Average Treatment Effect on the Treated",
    ddml_late = "Local Average Treatment Effect")
  model_name <- type_labels[x$model_type]
  if (is.na(model_name)) model_name <- x$model_type

  cat("DDML replicated fits:", model_name, "\n")
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
#' @param ... Additional arguments. Supports
#'     \code{aggregation} (\code{"median"} or \code{"mean"}).
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
#' @export
#' @method coef ddml_rep
coef.ddml_rep <- function(object,
                          aggregation = c("median",
                                          "mean"),
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
#' Returns a diagonal variance-covariance matrix from
#' cross-resample aggregation.
#'
#' @param object A \code{ddml_rep} object.
#' @param ensemble_idx Integer index of the ensemble type.
#'     Defaults to 1.
#' @param aggregation Character string, either
#'     \code{"median"} (default) or \code{"mean"}.
#' @param type Character string specifying the
#'     variance-covariance estimator. One of \code{"HC1"}
#'     (default), \code{"HC0"}, or \code{"HC3"}.
#' @param ... Currently unused.
#'
#' @return A diagonal p x p variance-covariance matrix.
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
#' @export
#' @method vcov ddml_rep
vcov.ddml_rep <- function(object, ensemble_idx = 1,
                          aggregation = c("median",
                                          "mean"),
                          type = "HC1", ...) {
  aggregation <- match.arg(aggregation)
  agg <- aggregate_reps(object,
                        aggregation = aggregation,
                        type = type)
  se_j <- agg$se[, ensemble_idx]
  V <- diag(se_j^2, nrow = length(se_j))
  rownames(V) <- colnames(V) <- object$coef_names
  V
}#VCOV.DDML_REP

#' Confidence Intervals for ddml_rep Objects
#'
#' @param object A \code{ddml_rep} object.
#' @param parm Not used (included for generic compatibility).
#' @param level Confidence level. Default 0.95.
#' @param ensemble_idx Integer index of the ensemble type.
#'     Defaults to 1.
#' @param aggregation Character string, either
#'     \code{"median"} (default) or \code{"mean"}.
#' @param type Character string specifying the
#'     variance-covariance estimator. One of \code{"HC1"}
#'     (default), \code{"HC0"}, or \code{"HC3"}.
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
#' @export
#' @method confint ddml_rep
confint.ddml_rep <- function(object, parm,
                             level = 0.95,
                             ensemble_idx = 1,
                             aggregation = c("median",
                                             "mean"),
                             type = "HC1", ...) {
  aggregation <- match.arg(aggregation)
  agg <- aggregate_reps(object,
                        aggregation = aggregation,
                        type = type)
  cf <- agg$coefficients[, ensemble_idx]
  se <- agg$se[, ensemble_idx]
  z <- stats::qnorm((1 + level) / 2)
  ci <- cbind(cf - z * se, cf + z * se)
  pct <- c((1 - level) / 2, (1 + level) / 2) * 100
  colnames(ci) <- paste0(format(pct, digits = 3), " %")
  rownames(ci) <- object$coef_names
  ci
}#CONFINT.DDML_REP

#' Summary for ddml_rep Objects
#'
#' @param object A \code{ddml_rep} object.
#' @param aggregation Character string, either
#'     \code{"median"} (default) or \code{"mean"}.
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
                             aggregation = c("median",
                                             "mean"),
                             type = "HC1", ...) {
  aggregation <- match.arg(aggregation)
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  agg <- aggregate_reps(object,
                        aggregation = aggregation,
                        type = type)
  inf_results <- build_inf_from_agg(
    agg, object$coef_names, object$ensemble_type)
  result <- list(
    inf_results   = inf_results,
    type          = type,
    model_type    = object$model_type,
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

#' Print Summary for ddml_rep Objects
#'
#' @param x An object of class \code{summary.ddml_rep}.
#' @param digits Number of significant digits. Default 3.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#'
#' @export
#' @method print summary.ddml_rep
print.summary.ddml_rep <- function(x, digits = 3, ...) {
  type_labels <- c(
    ddml_plm  = "Partially Linear Model",
    ddml_pliv = "Partially Linear IV Model",
    ddml_fpliv =
      "Flexible Partially Linear IV Model",
    ddml_ate  = "Average Treatment Effect",
    ddml_att  =
      "Average Treatment Effect on the Treated",
    ddml_late = "Local Average Treatment Effect")
  model_name <- type_labels[x$model_type]
  if (is.na(model_name)) model_name <- x$model_type

  cat("DDML estimation:", model_name, "\n")
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
#' @param aggregation Character string, either
#'     \code{"median"} (default) or \code{"mean"}.
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
#'     \code{p.value}, and \code{ensemble_type}. If
#'     \code{conf.int = TRUE}, also \code{conf.low} and
#'     \code{conf.high}.
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
#' @family ddml
#' @export
#' @method tidy ddml_rep
tidy.ddml_rep <- function(x, ensemble_idx = 1,
                          aggregation = c("median",
                                          "mean"),
                          type = "HC1",
                          conf.int = FALSE,
                          conf.level = 0.95, ...) {
  s <- summary(x, aggregation = aggregation,
               type = type)
  inf <- s$inf_results
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
#' @family ddml
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
    nresamples = x$nresamples,
    stringsAsFactors = FALSE
  )
}#GLANCE.DDML_REP
