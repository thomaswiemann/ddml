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
#' The influence function for \eqn{\gamma} is obtained via
#' the product rule:
#'
#' \deqn{\phi_\gamma(W_i; \theta, \eta, J, R)
#'   = R' \phi_\theta(W_i; \theta, \eta, J)
#'   + \phi_R(W_i)' \theta}
#'
#' where \eqn{\phi_\theta} is the influence function of
#' \eqn{\hat\theta} (see \code{\link{ddml-intro}}) and
#' \eqn{\phi_R(W_i)} is the \eqn{p}-vector of influence
#' functions for the weight vector \eqn{R}, i.e., the
#' \eqn{i}-th row of \code{inf_func_R}. The second term
#' vanishes when \eqn{R} is fixed.
#'
#' The leverage for \eqn{\gamma} is defined in
#' \code{\link{hatvalues.ddml_lincom}}.
#' The resulting \code{ddml_lincom} object supports all standard
#' inference methods: \code{vcov}, \code{confint}, \code{summary},
#' and \code{tidy}. For \code{ddml_rep} objects, \code{lincom}
#' can be called directly to obtain a \code{ddml_lincom_rep}.
#'
#' Note that \code{inf_func_R} is needed for inference when
#' \eqn{R} is estimated. Leverage computation further requires
#' \code{dinf_dR}. See \code{\link{vcov.ddml_lincom}} and
#' \code{\link{hatvalues.ddml_lincom}} for more details.
#'
#' @inheritParams ddml
#' @param fit A \code{ddml} or \code{ddml_rep} object.
#' @param R A \eqn{(p \times q)}{p x q} contrast matrix.
#'     Each column defines one linear combination.
#' @param ensemble_idx Integer index of the ensemble type to
#'     use. Default 1.
#' @param labels Optional character vector of length \eqn{q}
#'     naming the linear combinations. Defaults to column
#'     names of \code{R}, or \code{"lc1"}, \code{"lc2"}, etc.
#' @param inf_func_R An optional \eqn{(n \times p)}{n x p} matrix
#'     of influence functions \eqn{\phi_{R,i}} for the
#'     weight vector \eqn{R}. When supplied, a delta-method
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
#' @return An object of class \code{"ddml_lincom"} (for
#'     \code{ddml} input) or \code{"ddml_lincom_rep"} (for
#'     \code{ddml_rep} input).
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
                        ensemble_idx = 1,
                        labels = NULL, 
                        inf_func_R = NULL, 
                        dinf_dR = NULL,
                        ...) {
  R <- as.matrix(R)
  p <- nrow(fit$coefficients)
  q <- ncol(R)
  n <- fit$nobs
  stopifnot(nrow(R) == p)

  # Point estimate: R'theta
  att <- fit$coefficients[, ensemble_idx]
  lincom_coef <- as.numeric(crossprod(R, att))

  # Unit-level IF for theta: phi_i
  theta_if <- fit$inf_func[, , ensemble_idx, drop = FALSE]
  dim(theta_if) <- dim(theta_if)[1:2]

  # Combined IF for R'theta
  # Term 1: R' phi_i^theta
  lincom_if <- theta_if %*% R  # n x q
  # Term 2: delta-method correction (when R estimated)
  if (!is.null(inf_func_R)) {
    stopifnot(nrow(inf_func_R) == n, ncol(inf_func_R) == p)
    lincom_if <- lincom_if + inf_func_R %*% (att * (R != 0))
  }#IF

  # Map Structural Leverage 
  dinf_dtheta <- NULL
  if (!is.null(fit$dinf_dtheta)) {
    dinf_j <- fit$dinf_dtheta[, , , ensemble_idx, drop = FALSE]
    dinf_lc <- array(NA_real_, dim = c(n, q, q, 1))
    for (i in seq_len(n)) {
      dinf_i <- matrix(dinf_j[i, , , 1], p, p)
      # Structural Leverage: R^T * dinf_i * R
      dinf_lc[i, , , 1] <- t(R) %*% dinf_i %*% R
    }#FOR
    dinf_dtheta <- dinf_lc

    if (!is.null(dinf_dR)) {
      stopifnot(dim(dinf_dR)[1:3] == c(n, q, q))
      if (length(dim(dinf_dR)) == 3) dim(dinf_dR) <- c(n, q, q, 1)
      dinf_dtheta <- dinf_dtheta + dinf_dR
    }#IF
  }#IF

  # Labels
  if (is.null(labels)) labels <- colnames(R)
  if (is.null(labels)) labels <- paste0("lc", seq_len(q))

  # Package as standalone ddml_lincom
  coef_mat <- matrix(lincom_coef, nrow = q, ncol = 1)
  rownames(coef_mat) <- labels
  colnames(coef_mat) <- "lincom"

  obj <- list(
    coefficients     = coef_mat,
    scores           = array(-lincom_if, dim = c(n, q, 1)),
    J                = array(diag(-1, q), dim = c(q, q, 1)),
    inf_func         = array(lincom_if, dim = c(n, q, 1)),
    dinf_dtheta      = dinf_dtheta,
    coef_names       = labels,
    nobs             = n,
    estimator_name   = "Linear Combination",
    ensemble_type    = "lincom",
    cluster_variable = fit$cluster_variable,
    fixed_R          = is.null(inf_func_R)
  )
  class(obj) <- "ddml_lincom"
  obj
}#LINCOM.DDML

# Internal helper ==============================================================

# Build a temporary ddml object for delegation.
lincom_to_ddml <- function(lc) {
  # dinf_dtheta evaluates to NULL by default
  q <- nrow(lc$coefficients)
  ddml(
    coefficients     = lc$coefficients,
    scores           = lc$scores,
    J                = lc$J,
    inf_func         = lc$inf_func,
    dinf_dtheta      = lc$dinf_dtheta,
    nobs             = lc$nobs,
    coef_names       = lc$coef_names,
    estimator_name   = lc$estimator_name,
    ensemble_type    = lc$ensemble_type,
    cluster_variable = lc$cluster_variable,
    call             = NULL
  )
}#LINCOM_TO_DDML

# Build a temporary ddml_rep from per-rep ddml_lincom fits.
lincom_to_ddml_rep <- function(lc_rep) {
  ddml_rep(lapply(lc_rep$fits, lincom_to_ddml))
}#LINCOM_TO_DDML_REP

# ddml_lincom S3 methods =======================================================

#' @rdname lincom
#' @param object A \code{ddml_lincom} object.
#' @export
#' @method coef ddml_lincom
coef.ddml_lincom <- function(object, ...) {
  cf <- object$coefficients[, 1]
  names(cf) <- object$coef_names
  cf
}#COEF.DDML_LINCOM

#' @rdname lincom
#' @export
#' @method nobs ddml_lincom
nobs.ddml_lincom <- function(object, ...) object$nobs

#' @rdname lincom
#' @details
#' \code{vcov.ddml_lincom} computes the sandwich variance of
#' \eqn{\lambda = R'\hat\theta} using the influence function
#' \eqn{\phi_{\lambda,i} = R' \phi_{\theta,i}
#'   + \sum_c \theta_c \phi_{R_c,i}},
#' where the second term is included only when
#' \code{inf_func_R} was supplied to \code{lincom}.
#' @param type Character. Variance estimator type
#'     (\code{"HC0"} or \code{"HC1"}). \code{"HC3"} is not
#'     supported for linear combinations.
#' @export
#' @method vcov ddml_lincom
vcov.ddml_lincom <- function(object, type = "HC1", ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  if (type == "HC3") {
    stop("HC3 standard errors (generalized leverage) are not ",
         "currently supported for linear combinations.", call. = FALSE)
  }#IF
  vcov(lincom_to_ddml(object), ensemble_idx = 1, type = type)
}#VCOV.DDML_LINCOM

#' @rdname lincom
#' @param parm Parameters for confidence intervals.
#' @param level Confidence level. Default 0.95.
#' @param uniform Logical. If \code{TRUE}, computes uniform
#'     confidence bands via the multiplier bootstrap.
#' @param bootstraps Integer bootstrap draws. Default 999.
#' @export
#' @method confint ddml_lincom
confint.ddml_lincom <- function(object, parm, level = 0.95,
                                 type = "HC1",
                                 uniform = FALSE,
                                 bootstraps = 999L, ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  if (type == "HC3") {
    stop("HC3 standard errors (generalized leverage) are not ",
         "currently supported for linear combinations.", call. = FALSE)
  }#IF
  confint(lincom_to_ddml(object), parm = parm, level = level,
          ensemble_idx = 1, type = type, uniform = uniform,
          bootstraps = bootstraps)
}#CONFINT.DDML_LINCOM

#' @rdname lincom
#' @export
#' @method summary ddml_lincom
summary.ddml_lincom <- function(object, type = "HC1", ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  if (type == "HC3") {
    stop("HC3 standard errors (generalized leverage) are not ",
         "currently supported for linear combinations.", call. = FALSE)
  }#IF
  summary(lincom_to_ddml(object), type = type)
}#SUMMARY.DDML_LINCOM

#' @rdname lincom
#' @param x A \code{ddml_lincom} or \code{ddml_lincom_rep}
#'     object.
#' @param conf.int Logical. Include confidence intervals?
#' @param conf.level Confidence level for intervals.
#' @export
#' @method tidy ddml_lincom
tidy.ddml_lincom <- function(x, conf.int = FALSE,
                              conf.level = 0.95,
                              type = "HC1",
                              uniform = FALSE,
                              bootstraps = 999L, ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  if (type == "HC3") {
    stop("HC3 standard errors (generalized leverage) are not ",
         "currently supported for linear combinations.", call. = FALSE)
  }#IF
  tidy(lincom_to_ddml(x), ensemble_idx = 1,
       conf.int = conf.int, conf.level = conf.level,
       type = type, uniform = uniform, bootstraps = bootstraps)
}#TIDY.DDML_LINCOM

#' @rdname lincom
#' @details
#' \code{hatvalues.ddml_lincom} computes the observation-level
#' leverage for \eqn{\gamma}. The leverage decomposes as:
#'
#' \deqn{h_\gamma(W_i; \theta, \eta, J, R)
#'   = \underbrace{R'
#'   \left(-\frac{1}{n}
#'   \frac{\partial \phi_\theta(W_i; \theta,
#'   \eta, J)}{\partial \theta}\right)
#'   R}_{\text{structural leverage}}
#'   + \underbrace{\theta'
#'   \left(-\frac{1}{n}
#'   \frac{\partial \phi_R(W_i)}
#'   {\partial R}\right)
#'   \theta}_{\text{weighting leverage}}}
#'
#' The structural leverage is computed from the parent
#' model's \code{dinf_dtheta}; the weighting leverage
#' requires \code{dinf_dR} to be supplied to \code{lincom}.
#' @param model A \code{ddml_lincom} or \code{ddml_lincom_rep} object.
#' @export
#' @method hatvalues ddml_lincom
hatvalues.ddml_lincom <- function(model, ...) {
  stop("HC3 standard errors (generalized leverage) are not currently ",
       "supported for linear combinations.", call. = FALSE)
}#HATVALUES.DDML_LINCOM

# lincom.ddml_rep ==============================================================

#' @rdname lincom
#' @export
#' @method lincom ddml_rep
lincom.ddml_rep <- function(fit, R, inf_func_R = NULL, dinf_dR = NULL,
                              ensemble_idx = 1,
                              labels = NULL, ...) {
  # Apply lincom.ddml to each rep
  lc_fits <- lapply(fit$fits, function(f) {
    lincom(f, R = R, inf_func_R = inf_func_R, dinf_dR = dinf_dR,
           ensemble_idx = ensemble_idx, labels = labels, ...)
  })

  obj <- list(
    fits           = lc_fits,
    nresamples     = length(lc_fits),
    nobs           = lc_fits[[1]]$nobs,
    coef_names     = lc_fits[[1]]$coef_names,
    estimator_name = "Linear Combination",
    ensemble_type  = "lincom",
    fixed_R        = is.null(inf_func_R)
  )
  class(obj) <- "ddml_lincom_rep"
  obj
}#LINCOM.DDML_REP

# ddml_lincom_rep S3 methods ===================================================

#' @rdname lincom
#' @param aggregation Aggregation method: \code{"median"},
#'     \code{"mean"}, or \code{"spectral"}.
#' @export
#' @method coef ddml_lincom_rep
coef.ddml_lincom_rep <- function(object, aggregation = c("median", "mean"), 
                                 ...) {
  aggregation <- match.arg(aggregation)
  coef(lincom_to_ddml_rep(object), aggregation = aggregation)
}#COEF.DDML_LINCOM_REP

#' @rdname lincom
#' @export
#' @method nobs ddml_lincom_rep
nobs.ddml_lincom_rep <- function(object, ...) object$nobs

#' @rdname lincom
#' @export
#' @method vcov ddml_lincom_rep
vcov.ddml_lincom_rep <- function(object,
                                 aggregation = c("median", "mean",
                                                 "spectral"),
                                 type = "HC1", ...) {
  aggregation <- match.arg(aggregation)
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  if (type == "HC3") {
    stop("HC3 standard errors (generalized leverage) are not ",
         "currently supported for linear combinations.", call. = FALSE)
  }#IF
  vcov(lincom_to_ddml_rep(object), aggregation = aggregation, type = type)
}#VCOV.DDML_LINCOM_REP

#' @rdname lincom
#' @export
#' @method confint ddml_lincom_rep
confint.ddml_lincom_rep <- function(object, parm, level = 0.95,
                                    aggregation = c("median", "mean",
                                                    "spectral"),
                                    type = "HC1",
                                    uniform = FALSE,
                                    bootstraps = 999L, ...) {
  aggregation <- match.arg(aggregation)
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  if (type == "HC3") {
    stop("HC3 standard errors (generalized leverage) are not ",
         "currently supported for linear combinations.", call. = FALSE)
  }#IF
  confint(lincom_to_ddml_rep(object), parm = parm, level = level,
          aggregation = aggregation, type = type,
          uniform = uniform, bootstraps = bootstraps)
}#CONFINT.DDML_LINCOM_REP

#' @rdname lincom
#' @export
#' @method summary ddml_lincom_rep
summary.ddml_lincom_rep <- function(object,
                                    aggregation = c("median", "mean",
                                                    "spectral"),
                                    type = "HC1", ...) {
  aggregation <- match.arg(aggregation)
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  if (type == "HC3") {
    stop("HC3 standard errors (generalized leverage) are not ",
         "currently supported for linear combinations.", call. = FALSE)
  }#IF
  summary(lincom_to_ddml_rep(object), aggregation = aggregation, type = type)
}#SUMMARY.DDML_LINCOM_REP

#' @rdname lincom
#' @export
#' @method tidy ddml_lincom_rep
tidy.ddml_lincom_rep <- function(x,
                                 aggregation = c("median", "mean"),
                                 type = "HC1",
                                 conf.int = FALSE,
                                 conf.level = 0.95,
                                 uniform = FALSE,
                                 bootstraps = 999L, ...) {
  aggregation <- match.arg(aggregation)
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  if (type == "HC3") {
    stop("HC3 standard errors (generalized leverage) are not ",
         "currently supported for linear combinations.", call. = FALSE)
  }#IF
  tidy(lincom_to_ddml_rep(x), aggregation = aggregation,
       type = type, conf.int = conf.int, conf.level = conf.level,
       uniform = uniform, bootstraps = bootstraps)
}#TIDY.DDML_LINCOM_REP

#' @rdname lincom
#' @export
#' @method hatvalues ddml_lincom_rep
hatvalues.ddml_lincom_rep <- function(model, ...) {
  stop("HC3 standard errors (generalized leverage) are not currently ",
       "supported for linear combinations.", call. = FALSE)
}#HATVALUES.DDML_LINCOM_REP

# Print methods ================================================================

#' @rdname lincom
#' @export
#' @method print ddml_lincom
print.ddml_lincom <- function(x, ...) {
  cat("DDML Linear Combination\n")
  cat("Obs:", x$nobs)
  if (!x$fixed_R) cat("  (delta-method)")
  cat("\n\n")
  cat("Use summary() for inference.\n")
  invisible(x)
}#PRINT.DDML_LINCOM

#' @rdname lincom
#' @export
#' @method print ddml_lincom_rep
print.ddml_lincom_rep <- function(x, ...) {
  cat("DDML Linear Combination (replicated)\n")
  cat("Obs:", x$nobs, "  Resamples:", x$nresamples)
  if (!is.null(x$fixed_R) && !x$fixed_R) cat("  (delta-method)")
  cat("\n\n")
  cat("Use summary() for aggregated inference.\n")
  invisible(x)
}#PRINT.DDML_LINCOM_REP
