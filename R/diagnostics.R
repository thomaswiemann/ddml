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
#' @param ... Currently unused.
#'
#' @return An object of class \code{ddml_diagnostics} containing per-equation
#'     learner diagnostics. Use \code{print()} for formatted output or
#'     \code{tidy()} for a flat data.frame.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' learners = list(list(what = ols),
#'                list(what = mdl_glmnet))
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = learners,
#'                     sample_folds = 2, silent = TRUE)
#' diagnostics(plm_fit, cvc = TRUE)
#' tidy(diagnostics(plm_fit, cvc = TRUE))
#' }
#'
#' @references
#' Lei J (2020). "Cross-Validation With Confidence." Journal of the American
#'     Statistical Association, 115(532), 1978-1997.
#'
#' @family utilities
#' @export
diagnostics <- function(object, cvc = FALSE,
                        bootnum = 500,
                        ...) {
  if (!inherits(object, "ddml")) {
    stop("object must be of class 'ddml'.", call. = FALSE)
  }#IF

  eq_names <- names(object$ensemble_weights)
  single_learner <- is_single_learner(object$learners)
  tables <- list()

  for (eq in eq_names) {
    w <- object$ensemble_weights[[eq]]
    m <- object$mspe[[eq]]
    r <- object$r2[[eq]]

    # Build per-learner table
    nlearners <- if (single_learner) 1L else nrow(w)
    
    # Per-learner (and per-ensemble) OOS mspe and r2
    m_vec <- if (is.null(m) || length(m) == 0) rep(NA_real_, nlearners) else as.numeric(m)
    r_vec <- if (is.null(r) || length(r) == 0) rep(NA_real_, nlearners) else as.numeric(r)
    n_total <- length(m_vec)

    # Use names assigned upstream, or fallback to default
    if (!is.null(names(m))) {
      learner_names <- names(m)
    } else {
      learner_names <- if (single_learner) "single" else paste0("learner_", seq_len(n_total))
    }#IFELSE

    tbl <- data.frame(
      learner = learner_names,
      mspe = m_vec,
      r2 = r_vec,
      stringsAsFactors = FALSE,
      row.names = NULL)

    if (single_learner) {
      tbl$weight <- 1
    } else {
      # Weights: average across folds if 3D array
      if (length(dim(w)) == 3) {
        w_avg <- apply(w, c(1, 2), mean)
      } else {
        w_avg <- as.matrix(w)
      }#IFELSE

      ens_names <- colnames(w_avg)
      if (is.null(ens_names)) {
        ens_names <- paste0("weight_", seq_len(ncol(w_avg)))
      }#IF

      for (j in seq_len(ncol(w_avg))) {
        col_name <- paste0("weight_", ens_names[j])
        # Pad with NA for ensemble rows
        w_col <- c(w_avg[, j], rep(NA_real_, n_total - nlearners))
        tbl[[col_name]] <- w_col
      }#FOR
    }#IFELSE

    # CVC p-values (opt-in)
    if (cvc && !single_learner) {
      cvc_pvals <- cvc_pvalues(object$fitted, object$splits, eq, bootnum)
      tbl$cvc_pval <- c(cvc_pvals, rep(NA_real_, n_total - length(cvc_pvals)))
    }#IF

    tables[[eq]] <- tbl
  }#FOR

  result <- list(
    tables = tables,
    model_type = class(object)[1],
    estimator_name = object$estimator_name,
    nobs = object$nobs,
    shortstack = object$shortstack,
    cvc = cvc)
  class(result) <- "ddml_diagnostics"
  result
}#DIAGNOSTICS

# Cross-validation comparison p-values (Lei, 2020).
#
# For each learner, tests whether it is significantly outperformed
# by any other learner using a sup-type multiplier bootstrap on
# the squared-residual differences.
#
# @param fitted  Named list of per-equation cross-fitted objects.
# @param splits  The splits structure from a ddml object.
# @param eq      Character name of the nuisance equation.
# @param bootnum Number of bootstrap replications.
#
# @return Numeric vector of p-values (length = nlearners),
#     or NA_real_ vector if residuals are unavailable.
cvc_pvalues <- function(fitted, splits, eq, bootnum = 500) {
  entry <- fitted[[eq]]
  resid <- if (!is.null(entry)) entry$cf_resid_bylearner else NULL
  subs <- splits[[eq]]$subsamples
  if (is.null(resid) || is.null(subs) || ncol(resid) < 2) {
    nL <- if (!is.null(resid)) ncol(resid) else 1L
    return(rep(NA_real_, nL))
  }#IF

  # Derive fold IDs from subsamples
  n <- nrow(resid)
  fid <- integer(n)
  for (k in seq_along(subs)) {
    fid[subs[[k]]] <- k
  }#FOR

  nlearners <- ncol(resid)
  pvalues <- numeric(nlearners)
  for (i in seq_len(nlearners)) {
    pvalues[i] <- cvc_one_vs_many(
      resid[, i], resid[, -i, drop = FALSE], fid, bootnum)
  }#FOR
  pvalues
}#CVC_PVALUES

# One-vs-many cross-validation comparison test.
#
# Tests whether learner i is dominated by any learner in the
# comparison set. Uses a sup-type statistic with multiplier
# bootstrap.
#
# Reference: Lei J (2020). "Cross-Validation With Confidence."
#   Journal of the American Statistical Association,
#   115(532), 1978-1997.
#
# @param resid_base   Numeric vector of OOS residuals from the
#     base learner being tested.
# @param resid_others Matrix (n x K) of OOS residuals from
#     the comparison learners.
# @param fid          Integer vector of fold IDs.
# @param bootnum      Number of bootstrap replications.
#
# @return P-value. Large p-value means the base learner is
#     not significantly worse than the best alternative.
cvc_one_vs_many <- function(resid_base, resid_others, fid,
                            bootnum = 500) {
  n <- length(resid_base)
  K <- ncol(resid_others)
  resid_base_mat <- matrix(resid_base, n, K)

  zeta <- resid_base_mat^2 - resid_others^2

  # Demean by fold
  fid_unique <- unique(fid)
  zeta_til <- zeta
  for (k in fid_unique) {
    sel <- which(fid == k)
    fold_means <- colMeans(zeta[sel, , drop = FALSE])
    zeta_til[sel, ] <- sweep(
      zeta[sel, , drop = FALSE], 2, fold_means)
  }#FOR

  zeta_m <- colMeans(zeta)
  zeta_sd <- apply(zeta_til, 2, stats::sd)
  zeta_sd[zeta_sd < .Machine$double.eps] <- Inf

  # Sup-type test statistic
  Tx <- max(sqrt(n) * zeta_m / zeta_sd)

  # Multiplier bootstrap
  zeta_scaled <- sweep(zeta_til, 2, zeta_sd, FUN = "/")
  W <- matrix(stats::rnorm(n * bootnum), n, bootnum)
  boot_stats <- crossprod(zeta_scaled, W) / sqrt(n)  # K × bootnum
  Txb <- apply(boot_stats, 2, max)

  mean(Txb > Tx)
}#CVC_ONE_VS_MANY

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
  model_name <- x$estimator_name
  if (is.null(model_name)) model_name <- x$model_type

  cat("Stacking diagnostics:", model_name, "\n")
  cat("Obs:", x$nobs, "\n\n")

  for (eq in names(x$tables)) {
    tbl <- x$tables[[eq]]
    cat("  ", eq, ":\n", sep = "")

    # Format numeric columns
    display <- tbl
    num_cols <- c("mspe", "r2", "cvc_pval",
      grep("^weight_", names(display), value = TRUE))
    for (col in num_cols) {
      if (col %in% names(display)) {
        display[[col]] <- round(display[[col]], digits)
      }#IF
    }#FOR

    print(display, row.names = FALSE, right = TRUE)
    cat("\n")
  }#FOR

  if (!is.null(x$shortstack) && x$shortstack) {
    if (x$cvc) cat("Note: CVC compares individual base learners.\n",
                   "      Shortstacked ensemble CVC is not available.\n")
    cat("Note: Ensemble MSPE and R2 for short-stacking rely on full-sample weights\n",
        "      and represent in-sample fit over cross-fitted base predictions.\n")
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
#' @return A \code{data.frame} with columns \code{equation}, \code{learner},
#'     \code{mspe}, \code{r2}, \code{weight}, and optionally \code{cvc_pval}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' learners = list(list(what = ols),
#'                list(what = mdl_glmnet))
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = learners,
#'                     sample_folds = 2, silent = TRUE)
#' tidy(diagnostics(plm_fit, cvc = TRUE))
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
