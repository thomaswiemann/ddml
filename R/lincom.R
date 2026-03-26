#' Linear Combinations of DDML Coefficients
#'
#' @family ddml inference
#' @importFrom stats coef confint
#'
#' @description Computes linear combinations \eqn{R'\hat\theta}
#'     of DDML coefficient estimates.
#'
#' @details
#' For a \eqn{p}-dimensional coefficient vector
#' \eqn{\hat\theta} and a \eqn{(p \times q)}{p x q} contrast
#' matrix \eqn{R}, the linear combination is
#' \eqn{\gamma = R'\hat\theta}.
#'
#' The influence function for \eqn{\gamma} is
#'
#' \deqn{\phi_\gamma(W_i; \theta, R)
#'   = R'\, \phi_\theta(W_i; \theta)
#'   + \Phi_{R}(W_i)\, \theta,}
#'
#' where \eqn{\phi_\theta} is the influence function of
#' \eqn{\hat\theta} (see \code{\link{ral}} and
#' \code{\link{ddml-intro}}),
#' \eqn{\Phi_R(W_i)} is the \eqn{(p \times q)}{p x q}
#' matrix of influence functions for the contrast matrix
#' \eqn{R} (the \eqn{i}-th slice of \code{inf_func_R}),
#' and the second term vanishes when \eqn{R} is fixed.
#' The estimated influence function is
#'
#' \deqn{\hat\phi_{\gamma,i}
#'   = R'\, \hat\phi_{\theta,i}
#'   + \hat\Phi_{R,i}\, \hat\theta.}
#'
#' The leverage for \eqn{\gamma} is
#'
#' \deqn{h_\gamma(W_i; \theta, R)
#'   = \mathrm{tr}\! \left(
#'   R'\,h_\theta(W_i; \theta)\,R
#'   + h_R(W_i)\right),}
#'
#' where \eqn{h_\theta} is the structural leverage from
#' the parent estimator (see \code{\link{hatvalues.ral}})
#' and \eqn{h_R} is the weighting leverage. The sample
#' analog is
#'
#' \deqn{\hat{h}_{\gamma,i}
#'   = \mathrm{tr}\! \left(
#'   R'\,\hat{h}_{\theta,i}\,R
#'   + \hat{h}_{R,i}\right),}
#'
#' where \eqn{\hat{h}_{\theta,i}} is mapped from the
#' parent's \code{dinf_dtheta} and \eqn{\hat{h}_{R,i}}
#' from the optional \code{dinf_dR} argument.
#'
#' The resulting \code{lincom} object inherits from
#' \code{ral} and supports all standard inference methods:
#' \code{vcov}, \code{confint}, \code{summary}, \code{tidy},
#' and \code{hatvalues}. For \code{ddml_rep} objects,
#' \code{lincom} returns a \code{lincom_rep} inheriting
#' from \code{ral_rep}.
#'
#' Note that \code{inf_func_R} is needed for inference when
#' \eqn{R} is estimated. Leverage computation further requires
#' \code{dinf_dR}. See \code{\link{vcov.ral}} and
#' \code{\link{hatvalues.ral}} for more details.
#'

#' @param fit A \code{ddml} or \code{ddml_rep} object.
#' @param R A \eqn{(p \times q)}{p x q} contrast matrix.
#'     Each column defines one linear combination.
#' @param fit_idx Integer index of the fit to use, or
#'     \code{NULL} (default) for all ensemble types.
#'     When \code{NULL}, the output carries all ensembles
#'     from the parent fit.
#' @param labels Optional character vector of length \eqn{q}
#'     naming the linear combinations. Defaults to column
#'     names of \code{R}, or \code{"lc1"}, \code{"lc2"}, etc.
#' @param inf_func_R An optional \eqn{(n \times p \times q)}
#'     {n x p x q} array of influence functions
#'     \eqn{\Phi_{R,i}} for the contrast matrix \eqn{R}.
#'     Slice \code{[,,k]} contains the IFs for column
#'     \eqn{k} of \eqn{R}. When supplied, a delta-method
#'     correction is applied to the variance. When
#'     \code{NULL} (default), \eqn{R} is treated as fixed
#'     and only the first term contributes.
#' @param dinf_dR An optional \eqn{(n \times q \times q)}{n x q x q}
#'     array of observation-level derivatives
#'     \eqn{-n^{-1} \partial \phi_{R,i} / \partial R},
#'     representing the Weighting Leverage. When supplied,
#'     it is added to the Structural Leverage to form the
#'     total leverage used by HC3.
#' @param ... Currently unused.
#'
#' @return An object of class \code{"lincom"} (inheriting
#'     from \code{"ral"}) for \code{ddml} input, or
#'     \code{"lincom_rep"} (inheriting from \code{"ral_rep"})
#'     for \code{ddml_rep} input.
#'
#' @seealso \code{\link{lincom_weights_did}} for constructing
#'     DiD aggregation weights.
#'
#' @references
#' Chernozhukov V, Chetverikov D, Demirer M, Duflo E,
#'     Hansen C B, Newey W, Robins J (2018).
#'     "Double/debiased machine learning for treatment
#'     and structural parameters." The Econometrics
#'     Journal, 21(1), C1-C68.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 200; T_ <- 4
#' X <- matrix(rnorm(n * 2), n, 2)
#' G <- sample(c(3, 4, Inf), n, replace = TRUE,
#'             prob = c(0.3, 0.3, 0.4))
#' y <- matrix(rnorm(n * T_), n, T_)
#' for (i in seq_len(n)) {
#'   if (is.finite(G[i])) {
#'     for (j in seq_len(T_)) {
#'       if (j >= G[i]) y[i, j] <- y[i, j] + 1
#'     }
#'   }
#' }
#' fit <- ddml_attgt(y, X, t = 1:T_, G = G,
#'                 learners = list(what = ols),
#'                 sample_folds = 2,
#'                 silent = TRUE)
#' # Simple contrast: first cell minus second
#' p <- nrow(fit$coefficients)
#' R <- matrix(0, p, 1)
#' R[1, 1] <- 1; R[2, 1] <- -1
#' lc <- lincom(fit, R = R, labels = "ATT1-ATT2")
#' summary(lc)
#' }
#'
#' @export
lincom <- function(fit, R, ...) {
  UseMethod("lincom")
}#LINCOM

#' @rdname lincom
#' @export
#' @method lincom ddml
lincom.ddml <- function(fit, R,
                        fit_idx = NULL,
                        labels = NULL,
                        inf_func_R = NULL,
                        dinf_dR = NULL,
                        ...) {
  R <- as.matrix(R)
  p <- nrow(fit$coefficients)
  q <- ncol(R)
  n <- fit$nobs
  stopifnot(nrow(R) == p)

  # Resolve ensemble indices
  nensb <- ncol(fit$coefficients)
  if (is.null(fit_idx)) {
    j_seq <- seq_len(nensb)
  } else {
    j_seq <- fit_idx
  }#IFELSE
  nfit <- length(j_seq)

  # Validate inf_func_R dimensions (ensemble-invariant)
  if (!is.null(inf_func_R)) {
    stopifnot(length(dim(inf_func_R)) == 3,
              dim(inf_func_R)[1] == n,
              dim(inf_func_R)[2] == p,
              dim(inf_func_R)[3] == q)
  }#IF

  # Allocate output arrays
  coef_mat <- matrix(NA_real_, q, nfit)
  lincom_if <- array(0, dim = c(n, q, nfit))
  dinf_dtheta <- NULL

  for (jj in seq_along(j_seq)) {
    j <- j_seq[jj]

    # Point estimate: R'theta_j
    theta <- fit$coefficients[, j]
    coef_mat[, jj] <- as.numeric(crossprod(R, theta))

    # Unit-level IF: R' phi_i^theta_j
    theta_if <- fit$inf_func[, , j, drop = FALSE]
    dim(theta_if) <- dim(theta_if)[1:2]
    lincom_if[, , jj] <- theta_if %*% R

    # Delta-method correction (when R estimated)
    if (!is.null(inf_func_R)) {
      for (k in seq_len(q)) {
        lincom_if[, k, jj] <- lincom_if[, k, jj] + inf_func_R[, , k] %*% theta
      }#FOR
    }#IF

    # Map Structural Leverage
    if (!is.null(fit$dinf_dtheta)) {
      if (is.null(dinf_dtheta)) {
        dinf_dtheta <- array(NA_real_, dim = c(n, q, q, nfit))
      }#IF
      dinf_j <- fit$dinf_dtheta[, , , j, drop = FALSE]
      for (i in seq_len(n)) {
        dinf_i <- matrix(dinf_j[i, , , 1], p, p)
        dinf_dtheta[i, , , jj] <- t(R) %*% dinf_i %*% R
      }#FOR

      # Add weighting leverage
      if (!is.null(dinf_dR)) {
        if (length(dim(dinf_dR)) == 3) {
          # 3D dinf_dR: broadcast to current ensemble
          dinf_dtheta[, , , jj] <- dinf_dtheta[, , , jj] +
            dinf_dR
        } else {
          # 4D dinf_dR: per-ensemble slice
          dinf_dtheta[, , , jj] <- dinf_dtheta[, , , jj] +
            dinf_dR[, , , jj]
        }#IFELSE
      }#IF
    }#IF
  }#FOR

  # Labels
  if (is.null(labels)) labels <- colnames(R)
  if (is.null(labels)) labels <- paste0("lc", seq_len(q))

  # Fit labels (ensemble type names)
  fit_labels <- colnames(fit$coefficients)[j_seq]
  if (is.null(fit_labels)) {
    fit_labels <- if (!is.null(fit$fit_labels)) {
      fit$fit_labels[j_seq]
    } else {
      paste0("lincom", seq_len(nfit))
    }#IFELSE
  }#IF

  rownames(coef_mat) <- labels
  colnames(coef_mat) <- fit_labels

  ral(coefficients = coef_mat,
      inf_func = lincom_if,
      dinf_dtheta = dinf_dtheta,
      nobs = n,
      coef_names = labels,
      cluster_variable = fit$cluster_variable,
      estimator_name = "Linear Combination",
      subclass = "lincom",
      fixed_R = is.null(inf_func_R))
}#LINCOM.DDML

# lincom.ddml_rep ==============================================================

#' @rdname lincom
#' @export
#' @method lincom ddml_rep
lincom.ddml_rep <- function(fit, R, inf_func_R = NULL,
                            dinf_dR = NULL,
                            fit_idx = NULL,
                            labels = NULL, ...) {
  # Apply lincom.ddml to each rep
  lc_fits <- lapply(fit$fits, function(f) {
    lincom(f, R = R, inf_func_R = inf_func_R,
           dinf_dR = dinf_dR, fit_idx = fit_idx,
           labels = labels, ...)
  })

  ral_rep(lc_fits, subclass = "lincom_rep",
          fixed_R = is.null(inf_func_R))
}#LINCOM.DDML_REP

# Print methods ================================================================

#' @rdname lincom
#' @param x A \code{lincom} or \code{lincom_rep} object.
#' @export
#' @method print lincom
print.lincom <- function(x, ...) {
  cat("Linear Combination\n")
  cat("Obs:", x$nobs)
  nfit <- ncol(x$coefficients)
  if (nfit > 1) cat("  Ensembles:", nfit)
  if (!x$fixed_R) cat("  (delta-method)")
  cat("\n\n")
  cat("Use summary() for inference.\n")
  invisible(x)
}#PRINT.LINCOM

#' @rdname lincom
#' @export
#' @method print lincom_rep
print.lincom_rep <- function(x, ...) {
  cat("Linear Combination (replicated)\n")
  cat("Obs:", x$nobs,
      "  Resamples:", x$nresamples)
  nfit <- x$nfit
  if (!is.null(nfit) && nfit > 1) cat("  Ensembles:", nfit)
  if (!is.null(x$fixed_R) && !x$fixed_R) {
    cat("  (delta-method)")
  }#IF
  cat("\n\n")
  cat("Use summary() for aggregated inference.\n")
  invisible(x)
}#PRINT.LINCOM_REP
