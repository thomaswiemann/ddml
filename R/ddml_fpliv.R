#' Estimator for the Flexible Partially Linear IV Coefficient
#'
#' @family ddml estimators
#'
#' @description Estimator for the flexible partially linear IV coefficient.
#'
#' @details
#' \strong{Parameter of Interest:} \code{ddml_fpliv} provides a Double/Debiased Machine Learning
#'     estimator for the target parameter \eqn{\theta_0} in the partially
#'     linear IV model given by:
#'
#' \deqn{Y = \theta_0 D + g_0(X) + U,}
#'
#' where \eqn{(Y, D, X, Z, U)} is a random vector such that
#'     \eqn{E[U\vert X, Z] = 0} and \eqn{E[Var(E[D\vert X, Z]\vert X)] \neq 0},
#'     and \eqn{g_0} is an unknown nuisance function.
#'
#' \strong{Neyman Orthogonal Score:} The Neyman orthogonal score is:
#'
#' \deqn{m(W; \theta, \eta) = [(Y - \ell(X)) - \theta(D - r(X))](v(X, Z) - r(X))}
#'
#' where the nuisance parameters are \eqn{\eta = (\ell, r, v)} taking
#'     true values \eqn{\ell_0(X) = E[Y|X]}, \eqn{r_0(X) = E[D|X]}, and \eqn{v_0(X, Z) = E[D|X, Z]}.
#'
#' \strong{Jacobian:}
#'
#' \deqn{J = -E[(D - r(X))(v(X, Z) - r(X))^\top]}
#'
#' See \code{\link{ddml-intro}} for how the influence function
#' and inference are derived from these components.
#'
#' @inheritParams ddml-intro
#' @param Z A (sparse) matrix of instruments.
#' @param learners_DXZ,learners_DX Optional arguments to allow for different
#'     base learners for estimation of \eqn{E[D \vert X, Z]}, \eqn{E[D \vert X]}. Setup is
#'     identical to \code{learners}.
#' @param custom_ensemble_weights_DXZ,custom_ensemble_weights_DX Optional
#'     arguments to allow for different
#'     custom ensemble weights for \code{learners_DXZ},\code{learners_DX}. Setup
#'     is identical to \code{custom_ensemble_weights}. Note:
#'     \code{custom_ensemble_weights} and
#'     \code{custom_ensemble_weights_DXZ},\code{custom_ensemble_weights_DX} must
#'     have the same number of columns.
#'
#' @return \code{ddml_fpliv} returns an object of S3 class
#'     \code{ddml_fpliv} and \code{ddml}. See
#'     \code{\link{ddml-intro}} for the common output structure.
#'     Additional pass-through fields: \code{learners},
#'     \code{learners_DXZ}, \code{learners_DX}.
#' @export
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' Z = AE98[, "samesex", drop = FALSE]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the partially linear IV model using a single base learner, ridge.
#' fpliv_fit <- ddml_fpliv(y, D, Z, X,
#'                         learners = list(what = mdl_glmnet,
#'                                         args = list(alpha = 0)),
#'                         sample_folds = 2,
#'                         silent = TRUE)
#' summary(fpliv_fit)
ddml_fpliv <- function(y, D, Z, X,
                       learners,
                       learners_DXZ = learners,
                       learners_DX = learners,
                       sample_folds = 10,
                       ensemble_type = "nnls",
                       shortstack = FALSE,
                       cv_folds = 10,
                       custom_ensemble_weights = NULL,
                       custom_ensemble_weights_DXZ = custom_ensemble_weights,
                       custom_ensemble_weights_DX = custom_ensemble_weights,
                       cluster_variable = seq_along(y),
                       silent = FALSE,
                       parallel = NULL,
                       fitted = NULL,
                       splits = NULL,
                       save_crossval = TRUE,
                       ...) {
  cl <- match.call()

  # Preliminaries --------------------------------------------------------------

  dots <- list(...)
  messages <- resolve_messages(dots, "ddml_fpliv", list(
    y_X = "E[Y|X]"))

  validate_inputs(y = y, D = D, X = X, Z = Z,
                  learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type,
                  cluster_variable = cluster_variable)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DXZ,
                          learners_DXZ)
  validate_custom_weights(custom_ensemble_weights_DX,
                          learners_DX)

  nobs <- length(y)
  D <- as.matrix(D)
  nD <- ncol(D)

  validate_fitted_splits_pair(fitted, splits)

  indxs <- get_sample_splits(
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    subsamples = splits[[1]]$subsamples,
    cv_subsamples = splits[[1]]$cv_subsamples)
  check_subsamples(indxs$subsamples, NULL, stratify = FALSE)

  t0 <- proc.time()[3]
  announce_start(messages, parallel, silent)

  # Reduced-form estimation ----------------------------------------------------

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

  # E[D|X,Z] — merge X and Z, translate assign_Z to assign_X offsets
  nX <- ncol(X)
  nZ <- ncol(Z)
  XZ <- cbind(X, Z)
  if (!is_single_learner(learners_DXZ)) {
    learners_DXZ <- lapply(learners_DXZ, function(l) {
      az <- l$assign_Z
      if (is.null(az)) az <- seq_len(nZ)
      ax <- l$assign_X
      if (is.null(ax)) ax <- seq_len(nX)
      l$assign_X <- c(ax, nX + az)
      l$assign_Z <- NULL
      l
    })#LAPPLY
  }#IF
  # E[D|X,Z]
  D_XZ_res_list <- vector("list", nD)
  for (k in seq_len(nD)) {
    D_XZ_res_list[[k]] <- get_CEF(
      D[, k, drop = FALSE], XZ,
      learners = learners_DXZ,
      ensemble_type = ensemble_type,
      shortstack = shortstack,
      custom_ensemble_weights = custom_ensemble_weights_DXZ,
      subsamples = indxs$subsamples,
      cv_subsamples = indxs$cv_subsamples,
      silent = silent,
      label = paste0("E[D", k, "|X,Z]"),
      parallel = parallel,
      fitted = fitted[[paste0("D", k, "_XZ")]])
  }#FOR

  # E[D|X]
  D_X_res_list <- vector("list", nD)
  for (k in seq_len(nD)) {
    D_X_res_list[[k]] <- get_CEF(
      D[, k, drop = FALSE], X,
      learners = learners_DX,
      ensemble_type = ensemble_type,
      shortstack = shortstack,
      custom_ensemble_weights = custom_ensemble_weights_DX,
      subsamples = indxs$subsamples,
      cv_subsamples = indxs$cv_subsamples,
      silent = silent,
      label = paste0("E[D", k, "|X]"),
      parallel = parallel,
      fitted = fitted[[paste0("D", k, "_X")]])
  }#FOR

  ensemble_type <- y_X_res$ensemble_type
  nensb <- if (is.null(ensemble_type)) 1L
    else length(ensemble_type)

  # Target parameter & influence function --------------------------------------

  coef <- matrix(0, nD + 1, nensb)
  scores <- array(NA_real_, dim = c(nobs, nD + 1, nensb))
  J <- array(NA_real_, dim = c(nD + 1, nD + 1, nensb))
  inf_func <- array(NA_real_, dim = c(nobs, nD + 1, nensb))
  dinf_dtheta <- array(NA_real_, dim = c(nobs, nD + 1, nD + 1, nensb))

  cn_weights <- dimnames(y_X_res$weights)[[2]]
  colnames(coef) <- if (is.null(cn_weights)) ensemble_type else cn_weights

  for (j in seq_len(nensb)) {
    D_r <- D - get_cf_fitted(D_X_res_list, j)
    V_r <- get_cf_fitted(D_XZ_res_list, j) -
      get_cf_fitted(D_X_res_list, j)
    y_r <- y - cbind(y_X_res$cf_fitted)[, j]

    D_r_mat <- as.matrix(D_r)
    V_r_mat <- as.matrix(V_r)
    
    D_fit <- cbind(D_r_mat, 1)
    V_fit <- cbind(V_r_mat, 1)
    
    # 2SLS: First stage projection
    pi_hat <- qr.solve(V_fit, D_fit)
    D_hat <- V_fit %*% pi_hat

    # 2SLS: Second stage fit
    coef_iv_j <- as.vector(qr.solve(D_hat, y_r))
    coef[, j] <- coef_iv_j

    X_hat <- D_hat
    X_full <- D_fit
    e_j <- as.vector(y_r - X_full %*% coef_iv_j)

    scores[, , j] <- X_hat * e_j
    J[, , j] <- -crossprod(X_hat, X_full) / nobs

    J_inv <- csolve(matrix(J[, , j], nD + 1, nD + 1))
    inf_func[, , j] <- matrix(scores[, , j], nobs, nD + 1) %*% t(J_inv)
    
    U <- X_hat %*% t(J_inv)
    dinf_dtheta[, , , j] <- sapply(seq_len(nD + 1), function(k) {
      -X_full[, k] * U
    }, simplify = "array")
  }#FOR


  cn_j <- colnames(D)
  if (is.null(cn_j)) cn_j <- paste0("D", seq_len(nD))
  rownames(coef) <- c(cn_j, "(Intercept)")
  coef_names <- rownames(coef)

  # Output ---------------------------------------------------------------------

  ensemble_weights <- list(y_X = y_X_res$weights)
  mspe <- list(y_X = y_X_res$mspe)
  r2 <- list(y_X = y_X_res$r2)
  for (k in seq_len(nD)) {
    eq <- paste0("D", k, "_X")
    ensemble_weights[[eq]] <- D_X_res_list[[k]]$weights
    mspe[[eq]] <- D_X_res_list[[k]]$mspe
    r2[[eq]] <- D_X_res_list[[k]]$r2
  }#FOR
  for (k in seq_len(nD)) {
    eq <- paste0("D", k, "_XZ")
    ensemble_weights[[eq]] <- D_XZ_res_list[[k]]$weights
    mspe[[eq]] <- D_XZ_res_list[[k]]$mspe
    r2[[eq]] <- D_XZ_res_list[[k]]$r2
  }#FOR

  announce_finish(t0, messages, silent)

  ddml(
    coefficients = coef,
    ensemble_weights = ensemble_weights,
    mspe = mspe,
    r2 = r2,
    inf_func = inf_func, dinf_dtheta = dinf_dtheta,
    scores = scores, J = J,
    coef_names = coef_names,
    estimator_name = "Flexible Partially Linear IV Model",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    cluster_variable = cluster_variable,
    fitted = c(
      list(y_X = build_fitted_entry(y_X_res, save_crossval)),
      build_fitted_flat(D_XZ_res_list, save_crossval,
                        "D", "_XZ"),
      build_fitted_flat(D_X_res_list, save_crossval,
                        "D", "_X")),
    splits = stats::setNames(
      rep(list(list(subsamples = indxs$subsamples,
                    cv_subsamples = indxs$cv_subsamples)),
          length(ensemble_weights)),
      names(ensemble_weights)),
    call = cl,
    subclass = "ddml_fpliv",
    # ddml_fpliv-specific fields
    learners = learners,
    learners_DXZ = learners_DXZ,
    learners_DX = learners_DX
  )
}#DDML_FPLIV
