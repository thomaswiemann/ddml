# Estimate a conditional expectation function via cross-fitting.
#
# Dispatches to crosspred() (standard stacking) or shortstacking()
# depending on the shortstack flag. Returns out-of-sample fitted
# values, ensemble weights, and MSPE.
#
# @param y Outcome vector.
# @param X Feature matrix (may be sparse).
# @param learners List of base learner specifications.
# @param ensemble_type Character vector of ensemble types.
# @param shortstack Logical; use short-stacking if TRUE.
# @param subsamples List of sample fold indices.
# @param cv_subsamples List of cross-validation fold indices.
# @param parallel Optional list with parallel config.
get_CEF <- function(y, X,
                    learners,
                    ensemble_type,
                    shortstack,
                    custom_ensemble_weights = NULL,
                    subsamples,
                    cv_subsamples,
                    silent = FALSE,
                    label = NULL,
                    auxiliary_X = NULL,
                    shortstack_y = y,
                    parallel = NULL,
                    fitted = NULL) {
  t0 <- proc.time()[3]

  # Use pre-computed predictions if supplied
  if (!is.null(fitted)) {
    if (!is.null(label)) {
      info_msg("  Estimating ", label, "...", silent = silent)
    }#IF
    if (!is.null(fitted$cf_fitted_bylearner)) {
      # Rule 2: recompute ensemble from per-learner predictions
      if (!is.null(auxiliary_X) &&
          is.null(fitted$auxiliary_fitted_bylearner)) {
        stop(paste("When auxiliary predictions are required,",
                   "fitted objects must contain",
                   "'auxiliary_fitted_bylearner'."))
      }#IF
      res <- build_CEF_from_crossfit(
        y, fitted$cf_fitted_bylearner,
        ensemble_type, custom_ensemble_weights,
        cv_resid_byfold = fitted$cv_resid_byfold,
        subsamples = subsamples,
        auxiliary_fitted_bylearner =
          fitted$auxiliary_fitted_bylearner)
      return(res)
    }#IF
    if (!is.null(fitted$cf_fitted)) {
      # Rule 1: pre-ensembled predictions, use directly
      res <- list(cf_fitted = fitted$cf_fitted,
                  weights = NULL, mspe = NULL, r2 = NULL,
                  auxiliary_fitted = NULL,
                  cf_fitted_bylearner = NULL,
                  cf_resid_bylearner = NULL,
                  cv_resid_byfold = NULL)
      return(res)
    }#IF
  }#IF

  # Constant or empty y: return trivial predictions
  if (length(unique(y)) <= 1) {
    constant_val <- if (length(y) > 0) y[1] else 0
    nlearners <- if (is_single_learner(learners)) 1L
      else length(learners)
    ncustom <- if (!is.null(custom_ensemble_weights))
      ncol(custom_ensemble_weights) else 0L
    nensb <- length(ensemble_type) + ncustom
    n <- length(y)

    aux_fitted <- aux_fitted_bl <- NULL
    if (!is.null(auxiliary_X)) {
      aux_fitted <- lapply(auxiliary_X, function(ax)
        matrix(constant_val, nrow(ax), nensb))
      aux_fitted_bl <- lapply(auxiliary_X, function(ax)
        matrix(constant_val, nrow(ax), nlearners))
    }#IF

    if (!is.null(label) && label != "") {
      info_msg("  ", label,
               " .......... skipped (constant outcome)",
               silent = silent)
    }#IF

    return(list(
      cf_fitted = matrix(constant_val, n, nensb),
      weights = NULL, mspe = NULL, r2 = NULL,
      cf_fitted_bylearner = matrix(constant_val, n, nlearners),
      cf_resid_bylearner = matrix(0, n, nlearners),
      cv_resid_byfold = NULL,
      auxiliary_fitted = aux_fitted,
      auxiliary_fitted_bylearner = aux_fitted_bl))
  }#IF

  # Compute CEF via cross-fitting
  if (!is.null(label)) {
    info_msg("  Estimating ", label, "...", silent = silent)
  }#IF
  if (shortstack) {
    res <- shortstacking(y, X,
                         learners = learners,
                         ensemble_type = ensemble_type,
                         custom_ensemble_weights =
                           custom_ensemble_weights,
                         subsamples = subsamples,
                         silent = silent,
                         auxiliary_X = auxiliary_X,
                         shortstack_y = shortstack_y,
                         parallel = parallel)
  } else {
    res <- crosspred(y, X,
                     learners = learners,
                     ensemble_type = ensemble_type,
                     custom_ensemble_weights =
                       custom_ensemble_weights,
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
  nensb <- ncol(as.matrix(
    CEF_res_byD[[1]][[1]]$cf_fitted))
  sample_folds <- length(
    CEF_res_byD[[1]][[1]]$auxiliary_fitted)

  # Populate CEF
  g_X_byD <- array(0, dim = c(nobs, nensb, nCEF))
  for (d in seq_len(nCEF)) {
    g_X_byD[is_D[[d]], , d] <-
      CEF_res_byD[[d]][[1]]$cf_fitted
    for (k in seq_len(sample_folds)) {
      g_X_byD[aux_indx[[d]][[k]], , d] <-
        CEF_res_byD[[d]][[1]]$auxiliary_fitted[[k]]
    }#FOR
  }#FOR

  # return as array, third dimension is different levels of d
  g_X_byD
}#EXTRAPOLATE_CEF
