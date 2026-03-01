# Estimate a conditional expectation function via cross-fitting.
#
# Dispatches to crosspred() (standard stacking) or shortstacking()
# depending on the shortstack flag. Returns out-of-sample fitted
# values, ensemble weights, and MSPE.
#
# @param y Outcome vector.
# @param X Feature matrix (may be sparse).
# @param Z Optional instrument matrix.
# @param learners List of base learner specifications.
# @param ensemble_type Character vector of ensemble types.
# @param shortstack Logical; use short-stacking if TRUE.
# @param subsamples List of sample fold indices.
# @param cv_subsamples List of cross-validation fold indices.
# @param parallel Optional list with parallel config.
get_CEF <- function(y, X, Z = NULL,
                    learners,
                    ensemble_type,
                    shortstack,
                    custom_ensemble_weights = NULL,
                    compute_insample_predictions = FALSE,
                    subsamples,
                    cv_subsamples,
                    silent = FALSE,
                    label = NULL,
                    auxiliary_X = NULL,
                    shortstack_y = y,
                    parallel = NULL,
                    fitted = NULL) {
  t0 <- proc.time()[3]

  if (!is.null(label)) {
    info_msg("  Estimating ", label, "...", silent = silent)
  }#IF

  # Use pre-computed predictions if supplied
  if (!is.null(fitted)) {
    if (!is.null(auxiliary_X) &&
        is.null(fitted$auxiliary_fitted_bylearner)) {
      stop(paste("When auxiliary predictions are required, fitted objects",
                 "must contain 'auxiliary_fitted_bylearner'."))
    }#IF
    res <- build_CEF_from_crossfit(
      y, fitted$crossfit_fitted,
      ensemble_type, custom_ensemble_weights,
      crossval_resid = fitted$crossval_resid,
      subsamples = subsamples,
      auxiliary_fitted_bylearner = fitted$auxiliary_fitted_bylearner)
    return(res)
  }#IF

  # Compute CEF via cross-fitting
  if (shortstack) {
    res <- shortstacking(y, X, Z,
                         learners = learners,
                         ensemble_type = ensemble_type,
                         custom_ensemble_weights =
                           custom_ensemble_weights,
                         compute_insample_predictions =
                           compute_insample_predictions,
                         subsamples = subsamples,
                         silent = silent,
                         auxiliary_X = auxiliary_X,
                         shortstack_y = shortstack_y,
                         parallel = parallel)
  } else {
    res <- crosspred(y, X, Z,
                     learners = learners,
                     ensemble_type = ensemble_type,
                     custom_ensemble_weights =
                       custom_ensemble_weights,
                     compute_insample_predictions =
                       compute_insample_predictions,
                     subsamples = subsamples,
                     cv_subsamples = cv_subsamples,
                     silent = silent,
                     auxiliary_X = auxiliary_X,
                     parallel = parallel)
  }#IFELSE

  # Return estimates
  return(res)
}#GET_CEF

# Extrapolate CEF predictions across treatment levels.
#
# For each level d of D, populates out-of-sample and auxiliary
# predictions into an (nobs x nensb x nlevels) array.
#
# @param D Treatment vector.
# @param CEF_res_byD List of get_CEF results, one per D level.
# @param aux_indx Auxiliary sample indices for extrapolation.
extrapolate_CEF <- function(D, CEF_res_byD, aux_indx) {
  # Data parameters
  nCEF <- length(CEF_res_byD)
  nobs <- length(D)
  D_levels <- lapply(CEF_res_byD, function(x) x$d)
  is_D <- rep(list(NULL), nCEF)
  for (d in seq_len(nCEF)) is_D[[d]] <- which(D == D_levels[d])
  nensb <- ncol(as.matrix(CEF_res_byD[[1]][[1]]$oos_fitted))
  sample_folds <- length(CEF_res_byD[[1]][[1]]$auxiliary_fitted)

  # Populate CEF
  g_X_byD <- array(0, dim = c(nobs, nensb, nCEF))
  for (d in seq_len(nCEF)) {
    g_X_byD[is_D[[d]], , d] <- CEF_res_byD[[d]][[1]]$oos_fitted
    for (k in seq_len(sample_folds)) {
      g_X_byD[aux_indx[[d]][[k]], , d] <-
        CEF_res_byD[[d]][[1]]$auxiliary_fitted[[k]]
    }#FOR
  }#FOR

  # return as array, third dimension is different levels of d
  g_X_byD
}#EXTRAPOLATE_CEF
