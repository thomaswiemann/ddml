#' Estimator for the Partially Linear IV Model.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml()], [ddml::coef.ddml()],
#'     [ddml::confint.ddml()], [ddml::tidy.ddml()],
#'     [ddml::glance.ddml()], [ddml::diagnostics()],
#'     [AER::ivreg()]
#'
#' @description Estimator for the partially linear IV model.
#'
#' @details \code{ddml_pliv} provides a Double/Debiased Machine Learning
#'     estimator for the target parameter \eqn{\theta_0} in the partially
#'     linear IV model given by
#'
#' \eqn{Y = \theta_0D + g_0(X) + U,}
#'
#' where \eqn{(Y, D, X, Z, U)} is a random vector such that
#'     \eqn{E[Cov(U, Z\vert X)] = 0} and \eqn{E[Cov(D, Z\vert X)] \neq 0}, and
#'     \eqn{g_0} is an unknown nuisance function.
#'
#' In this model, the target parameter \eqn{\theta_0} is identified by the 
#'     estimating equation 
#'     \eqn{E[m(W; \theta_0, \eta_0)] = 0}, where \eqn{W = (Y, D, Z, X)} and
#'     \eqn{m(W; \theta, \eta)} is the Neyman orthogonal score
#'
#' \eqn{m(W; \theta, \eta) = (Y - \ell(X) - \theta(D - r_D(X)))(Z - r_Z(X)),}
#'
#'     with nuisance parameters \eqn{\eta = (\ell, r_D, r_Z)} taking true values
#'     \eqn{\ell_0(X) = E[Y|X]}, \eqn{r_{D,0}(X) = E[D|X]}, and
#'     \eqn{r_{Z,0}(X) = E[Z|X]}.
#'
#' @inheritParams ddml_plm
#' @param Z A matrix of instruments.
#' @param learners May take one of two forms, depending on whether a single
#'     learner or stacking with multiple learners is used for estimation of the
#'     conditional expectation functions.
#'     If a single learner is used, \code{learners} is a list with two named
#'     elements:
#'     \itemize{
#'         \item{\code{what} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}.}
#'     }
#'     If stacking with multiple learners is used, \code{learners} is a list of
#'     lists, each containing four named elements:
#'     \itemize{
#'         \item{\code{fun} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}.}
#'         \item{\code{assign_X} An optional vector of column indices
#'             corresponding to control variables in \code{X} that are passed to
#'             the base learner.}
#'         \item{\code{assign_Z} An optional vector of column indices
#'             corresponding to instruments in \code{Z} that are passed to the
#'             base learner.}
#'     }
#'     Omission of the \code{args} element results in default arguments being
#'     used in \code{fun}. Omission of \code{assign_X} (and/or \code{assign_Z})
#'     results in inclusion of all variables in \code{X} (and/or \code{Z}).
#' @param learners_DX,learners_ZX Optional arguments to allow for different
#'     base learners for estimation of \eqn{E[D|X]}, \eqn{E[Z|X]}. Setup is
#'     identical to \code{learners}.
#' @param custom_ensemble_weights_DX,custom_ensemble_weights_ZX Optional
#'     arguments to allow for different
#'     custom ensemble weights for \code{learners_DX},\code{learners_ZX}. Setup
#'     is identical to \code{custom_ensemble_weights}. Note:
#'     \code{custom_ensemble_weights} and
#'     \code{custom_ensemble_weights_DX},\code{custom_ensemble_weights_ZX} must
#'     have the same number of columns.
#'
#' @return \code{ddml_pliv} returns an object of S3 class
#'     \code{ddml_pliv}. An object of class \code{ddml_pliv} is a list
#'     containing the following components:
#'     \describe{
#'         \item{\code{coef}}{A vector with the \eqn{\theta_0} estimates and
#'             the second-stage intercept (last element).}
#'         \item{\code{ensemble_weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{iv_fit}}{Object of class \code{ivreg} from the IV
#'             regression of \eqn{Y - \hat{E}[Y\vert X]} on
#'             \eqn{D - \hat{E}[D\vert X]} using \eqn{Z - \hat{E}[Z\vert X]} as
#'             the instrument. See also [AER::ivreg()] for details.}
#'         \item{\code{learners},\code{learners_DX},\code{learners_ZX},
#'             \code{cluster_variable}, \code{subsamples},
#'             \code{cv_subsamples},\code{ensemble_type}}{Pass-through of
#'             selected user-provided arguments. See above.}
#'     }
#' @export
#'
#' @references
#' Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). "Model Averaging and 
#'     Double Machine Learning." Journal of Applied Econometrics, 40(3): 249-269.
#'
#' Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
#'     Robins J (2018). "Double/debiased machine learning for treatment and
#'     structural parameters." The Econometrics Journal, 21(1), C1-C68.
#'
#' Kleiber C, Zeileis A (2008). Applied Econometrics with R. Springer-Verlag,
#'     New York.
#'
#' Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2), 241-259.
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' Z = AE98[, "samesex"]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the partially linear IV model using a single base learner, ridge.
#' pliv_fit <- ddml_pliv(y, D, Z, X,
#'                       learners = list(what = mdl_glmnet,
#'                                       args = list(alpha = 0)),
#'                       sample_folds = 2,
#'                       silent = TRUE)
#' summary(pliv_fit)
ddml_pliv <- function(y, D, Z, X,
                      learners,
                      learners_DX = learners,
                      learners_ZX = learners,
                      sample_folds = 10,
                      ensemble_type = "nnls",
                      shortstack = FALSE,
                      cv_folds = 10,
                      custom_ensemble_weights = NULL,
                      custom_ensemble_weights_DX = custom_ensemble_weights,
                      custom_ensemble_weights_ZX = custom_ensemble_weights,
                      cluster_variable = seq_along(y),
                      silent = FALSE,
                      parallel = NULL,
                      fitted = NULL,
                      splits = NULL,
                      save_crossval = TRUE,
                      ...) {
  cl <- match.call()

  # == Preliminaries ================================================

  dots <- list(...)
  messages <- resolve_messages(dots, "ddml_pliv", list(
    y_X = "E[Y|X]"))

  validate_inputs(y = y, D = D, X = X, Z = Z,
                  learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX,
                          learners_DX)
  validate_custom_weights(custom_ensemble_weights_ZX,
                          learners_ZX)

  nobs <- length(y)
  D <- as.matrix(D)
  nD <- ncol(D)
  Z <- as.matrix(Z)
  nZ <- ncol(Z)

  splits <- normalize_splits(splits = splits, ...)
  validate_fitted_splits_pair(fitted, splits, !shortstack)

  indxs <- get_sample_splits(
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    subsamples = splits$subsamples,
    cv_subsamples = splits$cv_subsamples)
  check_subsamples(indxs$subsamples, NULL, stratify = FALSE)

  t0 <- proc.time()[3]
  mode_str <- if (!is.null(parallel)) {
    p <- parse_parallel(parallel)
    paste0("parallel, ", p$num_cores, " cores")
  } else {
    "sequential"
  }#IFELSE
  if (!is.null(messages$start) && messages$start != "") {
    info_msg(sprintf(messages$start, mode_str),
             silent = silent)
  }#IF

  # == Reduced-form estimation ======================================

  # E[Y|X]
  y_X_res <- get_CEF(y, X,
                     learners = learners,
                     ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights =
                       custom_ensemble_weights,
                     subsamples = indxs$subsamples,
                     cv_subsamples = indxs$cv_subsamples,
                     silent = silent, label = messages$y_X,
                     parallel = parallel,
                     fitted = fitted$y_X)

  # E[Z|X]
  Z_X_res_list <- compute_CEF_list(
    Z, X, learners = learners_ZX,
    ensemble_type = ensemble_type,
    shortstack = shortstack,
    custom_ensemble_weights = custom_ensemble_weights_ZX,
    subsamples = indxs$subsamples,
    cv_subsamples = indxs$cv_subsamples,
    silent = silent,
    label_prefix = "E[Z", label_suffix = "|X]",
    parallel = parallel,
    fitted = fitted$Z_X)

  # E[D|X]
  D_X_res_list <- compute_CEF_list(
    D, X, learners = learners_DX,
    ensemble_type = ensemble_type,
    shortstack = shortstack,
    custom_ensemble_weights = custom_ensemble_weights_DX,
    subsamples = indxs$subsamples,
    cv_subsamples = indxs$cv_subsamples,
    silent = silent,
    label_prefix = "E[D", label_suffix = "|X]",
    parallel = parallel,
    fitted = fitted$D_X)

  ensb_info <- update_ensemble_info(y_X_res$weights)
  ensemble_type <- ensb_info$ensemble_type
  nensb <- ensb_info$nensb

  # == Score construction ===========================================

  coef <- matrix(0, nD + 1, nensb)
  iv_fit <- rep(list(1), nensb)
  scores <- vector("list", nensb)
  J_list <- vector("list", nensb)
  psi_a <- vector("list", nensb)
  psi_b <- vector("list", nensb)

  for (j in seq_len(nensb)) {
    y_r <- y - cbind(y_X_res$oos_fitted)[, j]
    D_r <- D - get_oosfitted(D_X_res_list, j)
    V_r <- Z - get_oosfitted(Z_X_res_list, j)

    iv_fit_j <- AER::ivreg(y_r ~ D_r | V_r, x = TRUE)

    coef_iv_j <- stats::coef(iv_fit_j)
    coef[, j] <- c(coef_iv_j[-1], coef_iv_j[1])
    iv_fit[[j]] <- iv_fit_j

    D_r_mat <- as.matrix(D_r)
    D_hat <- as.matrix(
      iv_fit_j$x$projected[, -1, drop = FALSE])
    X_hat <- cbind(D_hat, 1)
    X_full <- cbind(D_r_mat, 1)
    e_j <- as.vector(stats::residuals(iv_fit_j))
    scores[[j]] <- X_hat * e_j
    J_list[[j]] <- -crossprod(X_hat, X_full) / nobs

    psi_b[[j]] <- X_hat * as.vector(y_r)
    psi_a[[j]] <- -sapply(seq_len(nD + 1),
                          function(k) X_hat * X_full[, k],
                          simplify = "array")
  }#FOR

  # == Target parameter =============================================

  cn_weights <- dimnames(y_X_res$weights)[[2]]
  colnames(coef) <- names(iv_fit) <-
    if (is.null(cn_weights)) ensemble_type else cn_weights
  cn_j <- names(iv_fit_j$coefficients)
  rownames(coef) <- c(cn_j[-1], cn_j[1])
  coef_names <- rownames(coef)

  ensemble_weights <- list(y_X = y_X_res$weights)
  mspe <- list(y_X = y_X_res$mspe)
  r2 <- list(y_X = y_X_res$r2)
  for (k in seq_len(nD)) {
    eq <- paste0("D", k, "_X")
    ensemble_weights[[eq]] <- D_X_res_list[[k]]$weights
    mspe[[eq]] <- D_X_res_list[[k]]$mspe
    r2[[eq]] <- D_X_res_list[[k]]$r2
  }#FOR
  for (k in seq_len(nZ)) {
    eq <- paste0("Z", k, "_X")
    ensemble_weights[[eq]] <- Z_X_res_list[[k]]$weights
    mspe[[eq]] <- Z_X_res_list[[k]]$mspe
    r2[[eq]] <- Z_X_res_list[[k]]$r2
  }#FOR

  # == Output =======================================================

  ddml_fit <- list(
    coefficients = coef,
    iv_fit = iv_fit,
    ensemble_weights = ensemble_weights,
    mspe = mspe,
    r2 = r2,
    psi_a = psi_a, psi_b = psi_b,
    scores = scores, J = J_list,
    coef_names = coef_names,
    estimator_name = "Partially Linear IV Model",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    learners = learners,
    learners_DX = learners_DX,
    learners_ZX = learners_ZX,
    cluster_variable = cluster_variable,
    fitted = list(
      y_X = build_fitted_entry(y_X_res, save_crossval),
      D_X = build_fitted_from_list(D_X_res_list,
                                   save_crossval),
      Z_X = build_fitted_from_list(Z_X_res_list,
                                   save_crossval)),
    splits = list(subsamples = indxs$subsamples,
                  cv_subsamples = indxs$cv_subsamples),
    call = cl)

  elapsed <- round(proc.time()[3] - t0, 1)
  if (!is.null(messages$finish) && messages$finish != "") {
    info_msg(sprintf(messages$finish, elapsed),
             silent = silent)
  }#IF

  class(ddml_fit) <- c("ddml_pliv", "ddml")
  return(ddml_fit)
}#DDML_PLIV
