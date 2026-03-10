#' @rdname ddml_ate
#'
#' @export
ddml_att <- function(y, D, X,
                     learners,
                     learners_DX = learners,
                     sample_folds = 10,
                     ensemble_type = "nnls",
                     shortstack = FALSE,
                     cv_folds = 10,
                     custom_ensemble_weights = NULL,
                     custom_ensemble_weights_DX = custom_ensemble_weights,
                     cluster_variable = seq_along(y),
                     stratify = TRUE,
                     trim = 0.01,
                     silent = FALSE,
                     parallel = NULL,
                     fitted = NULL,
                     splits = NULL,
                     save_crossval = TRUE,
                     ...) {
  cl <- match.call()

  # == Preliminaries ================================================

  dots <- list(...)
  messages <- resolve_messages(dots, "ddml_att", list(
    y_D0 = "E[Y|D=0,X]",
    D_X = "E[D|X]",
    D = "E[D]"))

  validate_inputs(y = y, D = D, X = X, learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type, trim = trim,
                  cluster_variable = cluster_variable,
                  require_binary_D = TRUE)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX,
                          learners_DX)

  nobs <- length(y)

  validate_fitted_splits_pair(fitted, splits)

  indxs <- get_sample_splits(
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    D = D, stratify = stratify,
    subsamples = splits$D_X$subsamples,
    subsamples_byD = list(
      splits$y_X_D0$subsamples,
      splits$y_X_D1$subsamples),
    cv_subsamples = splits$D_X$cv_subsamples,
    cv_subsamples_byD = list(
      splits$y_X_D0$cv_subsamples,
      splits$y_X_D1$cv_subsamples))
  check_subsamples(indxs$subsamples, indxs$subsamples_byD,
                   stratify, D)

  t0 <- proc.time()[3]
  announce_start(messages, parallel, silent)

  # == Reduced-form estimation ======================================

  # E[D|X]
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX,
                     ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights =
                       custom_ensemble_weights_DX,
                     subsamples = indxs$subsamples,
                     cv_subsamples = indxs$cv_subsamples,
                     silent = silent, label = messages$D_X,
                     parallel = parallel,
                     fitted = fitted$D_X)

  # E[D] (unconditional treatment probability, cross-fitted)
  D_res <- get_CEF(D, matrix(1, nobs, 1),
                   learners = list(what = ols,
                                   args = list(const = FALSE)),
                   ensemble_type = "average",
                   shortstack = FALSE,
                   subsamples = indxs$subsamples,
                   cv_subsamples = indxs$cv_subsamples,
                   silent = silent, label = messages$D,
                   fitted = fitted$D)

  # E[Y|D=0,X]
  is_d0 <- which(D == 0)
  y_X_D0_res <- get_CEF(y[is_d0], X[is_d0, , drop = FALSE],
                         learners = learners,
                         ensemble_type = ensemble_type,
                         shortstack = shortstack,
                         custom_ensemble_weights =
                           custom_ensemble_weights,
                         subsamples = indxs$subsamples_byD[[1]],
                         cv_subsamples =
                           indxs$cv_subsamples_byD[[1]],
                         silent = silent,
                         label = messages$y_D0,
                         auxiliary_X = get_auxiliary_X(
                           indxs$aux_indx[[1]], X),
                         parallel = parallel,
                         fitted = fitted$y_X_D0)

  ensemble_type <- y_X_D0_res$ensemble_type
  nensb <- if (is.null(ensemble_type)) 1L
    else length(ensemble_type)

  # == Score construction ===========================================

  # Extrapolate E[Y|D=0,X] to full sample. aux_indx is indexed by
  # sorted D levels {0, 1}: [[1]] holds positions of {D=1}
  # observations per fold, where the D=0 model must extrapolate.
  g_X_D0 <- extrapolate_CEF(
    D = D,
    CEF_res_byD = list(list(fit = y_X_D0_res, d = 0)),
    aux_indx = indxs$aux_indx)[, , 1]

  m_X <- D_X_res$cf_fitted
  m_X_tr <- trim_propensity_scores(m_X, trim, ensemble_type)
  p <- D_res$cf_fitted[, 1]

  D_mat <- matrix(D, nobs, nensb)
  y_mat <- matrix(y, nobs, nensb)
  p_mat <- matrix(p, nobs, nensb)

  psi_b_mat <- D_mat * (y_mat - g_X_D0) / p_mat -
    m_X_tr * (1 - D_mat) * (y_mat - g_X_D0) / (p_mat * (1 - m_X_tr))
  psi_a_vec <- -D / p
  psi_b <- lapply(seq_len(nensb), function(j) psi_b_mat[, j, drop = FALSE])
  psi_a <- lapply(seq_len(nensb), function(j) array(psi_a_vec, dim = c(nobs, 1, 1)))

  # == Target parameter =============================================

  att <- -colMeans(psi_b_mat) / mean(psi_a_vec)

  scores <- array(NA_real_, dim = c(nobs, 1, nensb))
  J <- array(NA_real_, dim = c(1, 1, nensb))
  for (j in seq_len(nensb)) {
    scores[, 1, j] <- psi_a_vec * att[j] + psi_b_mat[, j]
    J[1, 1, j] <- mean(psi_a_vec)
  }#FOR

  coef_names <- "ATT"
  coef <- matrix(att, nrow = 1, ncol = nensb)
  rownames(coef) <- coef_names
  colnames(coef) <- ensemble_type

  # == Output =======================================================

  announce_finish(t0, messages, silent)

  ddml(
    coefficients = coef,
    ensemble_weights = list(
      y_X_D0 = y_X_D0_res$weights,
      D_X = D_X_res$weights,
      D = D_res$weights),
    mspe = list(y_X_D0 = y_X_D0_res$mspe,
                D_X = D_X_res$mspe,
                D = D_res$mspe),
    r2 = list(y_X_D0 = y_X_D0_res$r2,
              D_X = D_X_res$r2,
              D = D_res$r2),
    psi_a = psi_a, psi_b = psi_b,
    scores = scores, J = J,
    coef_names = coef_names,
    estimator_name =
      "Average Treatment Effect on the Treated",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    cluster_variable = cluster_variable,
    fitted = list(
      y_X_D0 = build_fitted_entry(y_X_D0_res,
                                   save_crossval,
                                   include_auxiliary = TRUE),
      D_X = build_fitted_entry(D_X_res, save_crossval),
      D = build_fitted_entry(D_res, save_crossval)),
    splits = list(
      y_X_D0 = list(
        subsamples = indxs$subsamples_byD[[1]],
        cv_subsamples = indxs$cv_subsamples_byD[[1]]),
      y_X_D1 = list(
        subsamples = indxs$subsamples_byD[[2]],
        cv_subsamples = indxs$cv_subsamples_byD[[2]]),
      D_X = list(
        subsamples = indxs$subsamples,
        cv_subsamples = indxs$cv_subsamples),
      D = list(
        subsamples = indxs$subsamples,
        cv_subsamples = indxs$cv_subsamples)),
    call = cl,
    subclass = "ddml_att",
    learners = learners,
    learners_DX = learners_DX
  )
}#DDML_ATT
