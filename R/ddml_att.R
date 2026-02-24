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
                     subsamples = NULL,
                     subsamples_byD = NULL,
                     cv_subsamples = NULL,
                     cv_subsamples_byD = NULL,
                     trim = 0.01,
                     silent = FALSE,
                     parallel = NULL) {
  # Validate inputs
  validate_inputs(y = y, D = D, X = X, learners = learners,
                  sample_folds = sample_folds, cv_folds = cv_folds,
                  ensemble_type = ensemble_type, trim = trim,
                  require_binary_D = TRUE)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX, learners_DX)

  # Data parameters
  nobs <- length(y)
  is_D0 <- which(D == 0)

  # Check whether ddml uses conventional stacking w/ data driven weights
  w_cv <- !shortstack &
    any(ensemble_type %in% c("nnls", "nnls1", "singlebest", "ols")) &
    (class(learners[[1]]) != "function")

  # Create crossfitting and cv tuples
  indxs <- get_sample_splits(cluster_variable = cluster_variable,
                             sample_folds = sample_folds,
                             cv_folds = if (w_cv) cv_folds,
                             D = D, stratify = stratify,
                             subsamples = subsamples,
                             subsamples_byD = subsamples_byD,
                             cv_subsamples = cv_subsamples,
                             cv_subsamples_byD = cv_subsamples_byD)
  check_subsamples(indxs$subsamples, indxs$subsamples_byD,
                   stratify, D)

  # Estimation start
  t0 <- proc.time()[3]
  mode_str <- if (!is.null(parallel)) {
    p <- parse_parallel(parallel)
    paste0("parallel, ", p$num_cores, " cores")
  } else {
    "sequential"
  }
  info_msg("ddml_att: estimating (", mode_str, ")",
           silent = silent)

  # Compute estimates of E[y|D=0,X]
  y_X_D0_res <- get_CEF(y[is_D0], X[is_D0, , drop = FALSE],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        subsamples = indxs$subsamples_byD[[1]],
                        cv_subsamples = indxs$cv_subsamples_byD[[1]],
                        silent = silent, label = "E[Y|D=0,X]",
                        auxiliary_X = get_auxiliary_X(indxs$aux_indx[[1]], X),
                        parallel = parallel)

  # Compute estimates of E[D|X]
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights_DX,
                     subsamples = indxs$subsamples,
                     cv_subsamples = indxs$cv_subsamples,
                     silent = silent, label = "E[D|X]",
                     parallel = parallel)

  # Compute estimates of E[D] -- simple computation of averages here
  D_res <- get_CEF(D, matrix(1, nobs, 1),
                   learners = list(what = ols),
                   ensemble_type = "average",
                   shortstack = FALSE,
                   cv_subsamples = NULL,
                   subsamples = indxs$subsamples,
                   silent = TRUE,
                   label = "E[D]",
                   parallel = parallel)

  # Update ensemble type to account for (optional) custom weights
  ensb_info <- update_ensemble_info(y_X_D0_res$weights)
  ensemble_type <- ensb_info$ensemble_type
  nensb <- ensb_info$nensb
  multiple_ensembles <- ensb_info$multiple_ensembles

  # Construct reduced form variables
  g_X_D0<- extrapolate_CEF(D = D,
                             CEF_res_byD = list(list(fit = y_X_D0_res, d = 0)),
                             aux_indx = indxs$aux_indx)[, , 1]
  m_X <- D_X_res$oos_fitted

  # Trim propensity scores, return warnings
  m_X_tr <- trim_propensity_scores(m_X, trim, ensemble_type)

  # Compute the ATT using the constructed variables
  if (!multiple_ensembles) {
    p <- as.vector(D_res$oos_fitted)
    g0 <- as.vector(g_X_D0)
    m <- as.vector(m_X_tr)
    psi_b <- matrix(
      D * (y - g0) / p - m * (1 - D) * (y - g0) / (p * (1 - m)),
      nobs, 1)
    psi_a <- matrix(-D / p, nobs, 1)
    att <- -mean(psi_b) / mean(psi_a)
    names(att) <- ensemble_type
  } else {
    y_copy <- matrix(rep(y, nensb), nobs, nensb)
    D_copy <- matrix(rep(D, nensb), nobs, nensb)
    p_copy <- matrix(rep(D_res$oos_fitted, nensb), nobs, nensb)
    psi_b <- D_copy * (y_copy - g_X_D0) / p_copy -
      m_X_tr * (1 - D_copy) * (y_copy - g_X_D0) /
      (p_copy * (1 - m_X_tr))
    psi_a <- -D_copy / p_copy
    att <- -colMeans(psi_b) / colMeans(psi_a)
    names(att) <- ensemble_type
  }#IFELSE

  # Compute scores and Jacobian from psi_a/psi_b
  scores <- lapply(seq_len(nensb), function(j) {
    as.matrix(psi_a[, j] * att[j] + psi_b[, j])
  })
  J_list <- lapply(seq_len(nensb), function(j) {
    as.matrix(mean(psi_a[, j]))
  })
  coef_names <- "ATT"

  # Ensemble metrics
  weights <- list(y_X_D0 = y_X_D0_res$weights,
                  D_X = D_X_res$weights)
  mspe <- list(y_X_D0 = y_X_D0_res$mspe,
               D_X = D_X_res$mspe)
  r2 <- list(y_X_D0 = y_X_D0_res$r2,
             D_X = D_X_res$r2)

  # Predictions and per-learner residuals
  oos_pred <- list(EY_D0_X = g_X_D0,
                   ED_X = m_X,
                   ED = D_res$oos_fitted)
  oos_resid_bylearner <- list(
    y_X_D0 = y_X_D0_res$oos_resid_bylearner,
    D_X = D_X_res$oos_resid_bylearner)

  # Organize output
  ddml_fit <- list(att = att, weights = weights, mspe = mspe,
                   psi_a = psi_a, psi_b = psi_b,
                   oos_pred = oos_pred,
                   learners = learners,
                   learners_DX = learners_DX,
                   cluster_variable = cluster_variable,
                   subsamples_byD = indxs$subsamples_byD,
                   cv_subsamples_byD =
                     indxs$cv_subsamples_byD,
                   ensemble_type = ensemble_type,
                   coefficients = att,
                   scores = scores,
                   J = J_list,
                   coef_names = coef_names,
                   nobs = nobs,
                   sample_folds = sample_folds,
                   cv_folds = if (shortstack) NULL
                     else cv_folds,
                   shortstack = shortstack,
                   oos_resid_bylearner =
                     oos_resid_bylearner,
                   r2 = r2)

  # Print estimation completion
  elapsed <- round(proc.time()[3] - t0, 1)
  info_msg("ddml_att: completed in ", elapsed, "s",
           silent = silent)

  # Amend class and return
  class(ddml_fit) <- c("ddml_att", "ddml")
  return(ddml_fit)
}#DDML_ATT
