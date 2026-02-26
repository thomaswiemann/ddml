#' Stacking Diagnostics for DDML Estimators
#'
#' Computes per-learner diagnostics including MSPE, R-squared,
#' ensemble weights, and optionally cross-validation comparison
#' (CVC) p-values for each nuisance equation.
#'
#' @param object An object of class \code{ddml}.
#' @param cvc Logical. Compute CVC p-values via multiplier
#'     bootstrap? Default \code{FALSE}. CVC tests whether each
#'     learner is significantly outperformed by the others.
#' @param bootnum Number of bootstrap replications for CVC.
#'     Default 500. Ignored when \code{cvc = FALSE}.
#' @param alpha Significance level for the model confidence
#'     set. Default 0.05.
#' @param ... Currently unused.
#'
#' @return An object of class \code{ddml_diagnostics}
#'     containing per-equation learner diagnostics. Use
#'     \code{print()} for formatted output or \code{tidy()}
#'     for a flat data.frame.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = list(what = ols),
#'                     sample_folds = 2, silent = TRUE)
#' diagnostics(plm_fit)
#' tidy(diagnostics(plm_fit))
#' }
#'
#' @export
diagnostics <- function(object, cvc = FALSE,
                        bootnum = 500, alpha = 0.05,
                        ...) {
  if (!inherits(object, "ddml")) {
    stop("object must be of class 'ddml'.")
  }#IF

  eq_names <- names(object$weights)
  single_learner <- is_single_learner(object$learners)
  tables <- list()

  for (eq in eq_names) {
    w <- object$weights[[eq]]
    m <- object$mspe[[eq]]
    r <- object$r2[[eq]]

    # Build per-learner table
    if (single_learner) {
      m_val <- if (is.null(m) || length(m) == 0) {
        NA_real_
      } else {
        as.numeric(m[1])
      }#IFELSE
      r_val <- if (is.null(r) || length(r) == 0) {
        NA_real_
      } else {
        as.numeric(r[1])
      }#IFELSE
      tbl <- data.frame(
        learner = "single",
        mspe = m_val,
        r2 = r_val,
        weight = 1,
        stringsAsFactors = FALSE)
    } else {
      nlearners <- nrow(w)
      learner_names <- paste0("learner_", seq_len(nlearners))

      # Weights: average across folds if 3D array
      if (length(dim(w)) == 3) {
        w_avg <- apply(w, c(1, 2), mean)
      } else {
        w_avg <- as.matrix(w)
      }#IFELSE

      # MSPE: average across folds if matrix
      if (is.matrix(m) && ncol(m) > 1) {
        m_avg <- rowMeans(m)
      } else {
        m_avg <- as.numeric(m)
      }#IFELSE

      # R-squared: average across folds if matrix
      if (is.matrix(r) && ncol(r) > 1) {
        r_avg <- rowMeans(r)
      } else {
        r_avg <- as.numeric(r)
      }#IFELSE

      # Use first ensemble type weights for display
      w_display <- if (ncol(w_avg) >= 1) {
        w_avg[, 1]
      } else {
        rep(NA_real_, nlearners)
      }#IFELSE

      tbl <- data.frame(
        learner = learner_names,
        mspe = m_avg,
        r2 = r_avg,
        weight = w_display,
        stringsAsFactors = FALSE,
        row.names = NULL)
    }#IFELSE

    # CVC p-values (opt-in)
    if (cvc && !single_learner) {
      resid <- object$oos_resid_bylearner[[eq]]
      subs <- get_diag_subsamples(object, eq)
      if (!is.null(resid) && !is.null(subs) &&
          ncol(resid) > 1) {
        cvc_res <- cvc_test(resid, subs, bootnum, alpha)
        tbl$cvc_pval <- cvc_res$pvalues
        tbl$in_conf_set <- cvc_res$confidence_set
      } else {
        tbl$cvc_pval <- NA_real_
        tbl$in_conf_set <- NA
      }#IFELSE
    }#IF

    tables[[eq]] <- tbl
  }#FOR

  result <- list(
    tables = tables,
    model_type = class(object)[1],
    nobs = object$nobs,
    shortstack = object$shortstack,
    cvc = cvc,
    alpha = alpha)
  class(result) <- "ddml_diagnostics"
  result
}#DIAGNOSTICS

# Resolve the correct subsamples for a given equation.
get_diag_subsamples <- function(object, eq) {
  model <- class(object)[1]
  if (model %in% c("ddml_plm", "ddml_pliv",
                    "ddml_fpliv")) {
    return(object$subsamples)
  }#IF
  if (model %in% c("ddml_ate", "ddml_att")) {
    if (grepl("D0", eq)) return(object$subsamples_byD[[1]])
    if (grepl("D1", eq)) return(object$subsamples_byD[[2]])
    return(merge_subsamples(object$subsamples_byD))
  }#IF
  if (model == "ddml_late") {
    if (grepl("Z0", eq)) return(object$subsamples_byZ[[1]])
    if (grepl("Z1", eq)) return(object$subsamples_byZ[[2]])
    return(merge_subsamples(object$subsamples_byZ))
  }#IF
  object$subsamples
}#GET_DIAG_SUBSAMPLES

# Reconstruct full-sample folds from stratified folds.
merge_subsamples <- function(subsamples_by) {
  K <- length(subsamples_by[[1]])
  lapply(seq_len(K), function(k) {
    sort(unlist(lapply(subsamples_by,
                       function(s) s[[k]])))
  })
}#MERGE_SUBSAMPLES

#' Print Stacking Diagnostics
#'
#' @param x An object of class \code{ddml_diagnostics}.
#' @param digits Number of significant digits. Default 4.
#' @param ... Currently unused.
#'
#' @return \code{x}, invisibly.
#'
#' @export
print.ddml_diagnostics <- function(x, digits = 4, ...) {
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

  cat("Stacking diagnostics:", model_name, "\n")
  cat("Obs:", x$nobs, "\n\n")

  for (eq in names(x$tables)) {
    tbl <- x$tables[[eq]]
    cat("  ", eq, ":\n", sep = "")

    # Format numeric columns
    display <- tbl
    for (col in c("mspe", "r2", "weight", "cvc_pval")) {
      if (col %in% names(display)) {
        display[[col]] <- round(display[[col]], digits)
      }#IF
    }#FOR
    if ("in_conf_set" %in% names(display)) {
      display$in_conf_set <- NULL
    }#IF

    print(display, row.names = FALSE, right = TRUE)
    cat("\n")
  }#FOR

  if (x$cvc && !is.null(x$shortstack) && x$shortstack) {
    cat("Note: CVC compares individual base learners.",
        "\nShortstacked ensemble CVC is not available",
        "(weights use all data).\n")
  }#IF

  invisible(x)
}#PRINT.DDML_DIAGNOSTICS

#' Tidy Stacking Diagnostics
#'
#' Returns a flat data.frame of per-learner stacking
#' diagnostics for all nuisance equations. Suitable for
#' table creation with \code{kable()}, \code{gt()}, or
#' \code{modelsummary::datasummary()}.
#'
#' @param x An object of class \code{ddml_diagnostics}.
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns \code{equation},
#'     \code{learner}, \code{mspe}, \code{r2}, \code{weight},
#'     and optionally \code{cvc_pval} and
#'     \code{in_conf_set}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = list(what = ols),
#'                     sample_folds = 2, silent = TRUE)
#' tidy(diagnostics(plm_fit))
#' }
#'
#' @export
#' @method tidy ddml_diagnostics
tidy.ddml_diagnostics <- function(x, ...) {
  rows <- list()
  for (eq in names(x$tables)) {
    tbl <- x$tables[[eq]]
    tbl$equation <- eq
    rows[[length(rows) + 1]] <- tbl
  }#FOR
  out <- do.call(rbind, rows)
  # Reorder: equation first
  eq_col <- which(names(out) == "equation")
  out <- out[, c(eq_col, seq_along(out)[-eq_col]),
             drop = FALSE]
  rownames(out) <- NULL
  out
}#TIDY.DDML_DIAGNOSTICS
