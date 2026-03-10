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
    D_X = "E[D|X]"))

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

  validate_fitted_splits_pair(fitted, splits, !shortstack)

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
  mode_str <- if (!is.null(parallel)) {
    p <- parse_parallel(parallel)
    paste0("parallel, ", p$num_cores, " cores")
  } else {
    "sequential"
  }#IFELSE
  if (!is.null(messages$start) && messages$start != "") {
    info_msg(sprintf(messages$start, mode_str), silent = silent)
  }#IF

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

  # E[Y|D=0,X] via ddml_apo(d=0)
  # Swap byD indices: apo uses D_ind = 1*(D==0), so byD
  # levels are reversed relative to the ATT's D.
  shared_splits <- list(subsamples = indxs$subsamples,
                        cv_subsamples = indxs$cv_subsamples)
  apo_splits_0 <- list(
    y_X = list(subsamples = indxs$subsamples_byD[[1]],
               cv_subsamples = indxs$cv_subsamples_byD[[1]]),
    D_X = shared_splits)

  # Pre-ensembled propensity: P(D=0|X) = 1 - E[D|X]
  fitted_D_X_0 <- list(
    cf_fitted = 1 - D_X_res$cf_fitted)

  apo_0 <- ddml_apo(
    y = y, D = D, X = X, d = 0, weights = NULL,
    learners = learners, learners_DX = learners_DX,
    sample_folds = sample_folds, cv_folds = cv_folds,
    custom_ensemble_weights = custom_ensemble_weights,
    custom_ensemble_weights_DX = custom_ensemble_weights_DX,
    cluster_variable = cluster_variable,
    ensemble_type = ensemble_type, shortstack = shortstack,
    trim = trim, parallel = parallel, silent = silent,
    splits = apo_splits_0,
    fitted = list(y_X = fitted$y_X_D0,
                  D_X = fitted_D_X_0),
    messages = list(start = "", finish = "",
                    y_X = messages$y_D0, D_X = ""))

  ensemble_type <- apo_0$ensemble_type
  nensb <- ncol(apo_0$coefficients)

  # == Score construction ===========================================

  # Extrapolate E[Y|D=0,X] to full sample. aux_indx is indexed by
  # sorted D levels {0, 1}: [[1]] holds positions of {D=1}
  # observations per fold, where the D=0 model must extrapolate.
  g_X_D0 <- extrapolate_CEF(
    D = D,
    CEF_res_byD = list(list(
      fit = apo_0$fitted$y_X, d = 0)),
    aux_indx = indxs$aux_indx)[, , 1]

  m_X <- D_X_res$cf_fitted
  m_X_tr <- trim_propensity_scores(m_X, trim, ensemble_type)
  p <- mean(D)

  D_mat <- matrix(D, nobs, nensb)
  y_mat <- matrix(y, nobs, nensb)

  psi_b_mat <- D_mat * (y_mat - g_X_D0) / p -
    m_X_tr * (1 - D_mat) * (y_mat - g_X_D0) /
    (p * (1 - m_X_tr))
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

  elapsed <- round(proc.time()[3] - t0, 1)
  if (!is.null(messages$finish) && messages$finish != "") {
    info_msg(sprintf(messages$finish, elapsed),
             silent = silent)
  }#IF

  ddml(
    coefficients = coef,
    ensemble_weights = list(
      y_X_D0 = apo_0$ensemble_weights$y_X,
      D_X = D_X_res$weights),
    mspe = list(y_X_D0 = apo_0$mspe$y_X,
                D_X = D_X_res$mspe),
    r2 = list(y_X_D0 = apo_0$r2$y_X,
              D_X = D_X_res$r2),
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
      y_X_D0 = apo_0$fitted$y_X,
      D_X = build_fitted_entry(D_X_res, save_crossval)),
    splits = list(
      y_X_D0 = list(
        subsamples = indxs$subsamples_byD[[1]],
        cv_subsamples = indxs$cv_subsamples_byD[[1]]),
      y_X_D1 = list(
        subsamples = indxs$subsamples_byD[[2]],
        cv_subsamples = indxs$cv_subsamples_byD[[2]]),
      D_X = list(
        subsamples = indxs$subsamples,
        cv_subsamples = indxs$cv_subsamples)),
    call = cl,
    subclass = "ddml_att",
    learners = learners,
    learners_DX = learners_DX
  )
}#DDML_ATT
