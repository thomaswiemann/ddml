#' Estimator for the Partially Linear Model.
#'
#' @family ddml estimators
#'
#' @description Estimator for the partially linear model.
#'
#' @details \code{ddml_plm} provides a Double/Debiased Machine Learning
#'     estimator for the target parameter \eqn{\theta_0} in the partially
#'     linear model given by
#'
#' \eqn{Y = \theta_0D + g_0(X) + U,}
#'
#' where \eqn{(Y, D, X, U)} is a random vector such that
#'     \eqn{E[Cov(U, D\vert X)] = 0} and \eqn{E[Var(D\vert X)] \neq 0}, and
#'     \eqn{g_0} is an unknown nuisance function.
#'
#' In this model, the target parameter \eqn{\theta_0} is identified by the 
#'     estimating equation 
#'     \eqn{E[m(W; \theta_0, \eta_0)] = 0}, where \eqn{W = (Y, D, X)} and
#'     \eqn{m(W; \theta, \eta)} is the Neyman orthogonal score
#'
#' \eqn{m(W; \theta, \eta) = (Y - \ell(X) - \theta(D - r(X)))(D - r(X)),}
#'
#'     with nuisance parameters \eqn{\eta = (\ell, r)} taking true values
#'     \eqn{\ell_0(X) = E[Y|X]} and \eqn{r_0(X) = E[D|X]}.
#'
#' @inheritParams ddml-class
#' @param learners_DX Optional argument to allow for different estimators of
#'     \eqn{E[D|X]}. Setup is identical to \code{learners}.
#' @param custom_ensemble_weights_DX Optional argument to allow for different
#'     custom ensemble weights for \code{learners_DX}. Setup is identical to
#'     \code{custom_ensemble_weights}. Note: \code{custom_ensemble_weights} and
#'     \code{custom_ensemble_weights_DX} must have the same number of columns.
#'
#' @return \code{ddml_plm} returns an object of S3 class
#'     \code{ddml_plm} and \code{ddml}. See \code{\link{ddml-class}}
#'     for the common output structure. Additional pass-through
#'     fields: \code{learners}, \code{learners_DX}.
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
#' Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2), 241-259.
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the partially linear model using a single base learner, ridge.
#' plm_fit <- ddml_plm(y, D, X,
#'                     learners = list(what = mdl_glmnet,
#'                                     args = list(alpha = 0)),
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(plm_fit)
#'
#' # Estimate the partially linear model using short-stacking with base learners
#' #     ols, lasso, and ridge. We can also use custom_ensemble_weights
#' #     to estimate the ATE using every individual base learner.
#' weights_everylearner <- diag(1, 3)
#' colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
#' plm_fit <- ddml_plm(y, D, X,
#'                     learners = list(list(what = ols),
#'                                     list(what = mdl_glmnet),
#'                                     list(what = mdl_glmnet,
#'                                          args = list(alpha = 0))),
#'                     ensemble_type = 'nnls',
#'                     custom_ensemble_weights = weights_everylearner,
#'                     shortstack = TRUE,
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(plm_fit)
#'
#' \donttest{
#' # Re-estimate with a different ensemble type using pass-through
#' #     (skips cross-fitting, only recomputes ensemble weights).
#' plm_fit2 <- ddml_plm(y, D, X,
#'                      learners = list(list(what = ols),
#'                                      list(what = mdl_glmnet),
#'                                      list(what = mdl_glmnet,
#'                                           args = list(alpha = 0))),
#'                      ensemble_type = 'average',
#'                      shortstack = TRUE,
#'                      sample_folds = 2,
#'                      silent = TRUE,
#'                      fitted = plm_fit$fitted,
#'                      splits = plm_fit$splits)
#' summary(plm_fit2)
#' }
ddml_plm <- function(y, D, X,
                     learners,
                     learners_DX = learners,
                     sample_folds = 10,
                     ensemble_type = "nnls",
                     shortstack = FALSE,
                     cv_folds = 10,
                     custom_ensemble_weights = NULL,
                     custom_ensemble_weights_DX = custom_ensemble_weights,
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
  messages <- resolve_messages(dots, "ddml_plm", list(
    y_X = "E[Y|X]"))

  validate_inputs(y = y, D = D, X = X, learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type,
                  cluster_variable = cluster_variable)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX,
                          learners_DX)

  nobs <- length(y)
  D <- as.matrix(D)
  nD <- ncol(D)

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
                     custom_ensemble_weights = custom_ensemble_weights,
                     subsamples = indxs$subsamples,
                     cv_subsamples = indxs$cv_subsamples,
                     silent = silent, label = "E[Y|X]",
                     parallel = parallel,
                     fitted = fitted$y_X)

  # Compute estimates of E[D|X], loop through endogenous variables
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
  scores <- array(NA_real_, dim = c(nobs, nD + 1, nensb))
  J <- array(NA_real_, dim = c(nD + 1, nD + 1, nensb))
  psi_a <- vector("list", nensb)
  psi_b <- vector("list", nensb)
  for (j in seq_len(nensb)) {
    D_r <- D - get_cf_fitted(D_X_res_list, j)
    y_r <- y - cbind(y_X_res$cf_fitted)[, j]

    D_r_mat <- as.matrix(D_r)
    X_full <- cbind(D_r_mat, 1)
    
    coef_ols_j <- as.vector(qr.solve(X_full, y_r))
    coef[, j] <- coef_ols_j
    e_j <- as.vector(y_r - X_full %*% coef_ols_j)
    
    scores[, , j] <- X_full * e_j
    J[, , j] <- -crossprod(X_full) / nobs

    psi_b[[j]] <- X_full * as.vector(y_r)
    psi_a[[j]] <- -sapply(seq_len(nD + 1),
                          function(k) X_full * X_full[, k],
                          simplify = "array")
  }#FOR

  # == Target parameter =============================================

  cn_weights <- dimnames(y_X_res$weights)[[2]]
  colnames(coef) <- if (is.null(cn_weights)) ensemble_type else cn_weights
  cn_j <- colnames(D)
  if (is.null(cn_j)) cn_j <- paste0("D", seq_len(nD))
  rownames(coef) <- c(cn_j, "(Intercept)")
  coef_names <- rownames(coef)

  # == Output =======================================================

  # Ensemble metrics
  ensemble_weights <- list(y_X = y_X_res$weights)
  mspe <- list(y_X = y_X_res$mspe)
  r2 <- list(y_X = y_X_res$r2)
  for (k in seq_len(nD)) {
    eq <- paste0("D", k, "_X")
    ensemble_weights[[eq]] <- D_X_res_list[[k]]$weights
    mspe[[eq]] <- D_X_res_list[[k]]$mspe
    r2[[eq]] <- D_X_res_list[[k]]$r2
  }#FOR

  ddml_fit <- list(
    coefficients = coef,
    ensemble_weights = ensemble_weights,
    mspe = mspe,
    r2 = r2,
    psi_a = psi_a, psi_b = psi_b,
    scores = scores, J = J,
    coef_names = coef_names,
    estimator_name = "Partially Linear Model",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    learners = learners,
    learners_DX = learners_DX,
    cluster_variable = cluster_variable,
    fitted = list(
      y_X = build_fitted_entry(y_X_res, save_crossval),
      D_X = build_fitted_from_list(D_X_res_list,
                                   save_crossval)),
    splits = list(subsamples = indxs$subsamples,
                  cv_subsamples = indxs$cv_subsamples),
    call = cl)

  elapsed <- round(proc.time()[3] - t0, 1)
  if (!is.null(messages$finish) && messages$finish != "") {
    info_msg(sprintf(messages$finish, elapsed),
             silent = silent)
  }#IF

  class(ddml_fit) <- c("ddml_plm", "ddml")
  return(ddml_fit)
}#DDML_PLM
