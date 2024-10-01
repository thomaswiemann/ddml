# Wrapper for [ddml::crosspred()] and [ddml::shortstacking()].
get_CEF <- function(y, X, Z = NULL,
                    learners,
                    ensemble_type,
                    shortstack,
                    custom_ensemble_weights = NULL,
                    compute_insample_predictions = FALSE,
                    compute_predictions_bylearner = FALSE,
                    subsamples,
                    cv_subsamples_list,
                    silent = FALSE,
                    progress = NULL,
                    auxiliary_X = NULL,
                    shortstack_y = y) {
  # Compute CEF
  if (shortstack) {
    res <- shortstacking(y, X, Z,
                         learners = learners,
                         ensemble_type = ensemble_type,
                         custom_ensemble_weights = custom_ensemble_weights,
                         compute_insample_predictions =
                           compute_insample_predictions,
                         subsamples = subsamples,
                         silent = silent, progress = progress,
                         auxiliary_X = auxiliary_X,
                         shortstack_y = shortstack_y)
  } else {
    res <- crosspred(y, X, Z,
                     learners = learners,
                     ensemble_type = ensemble_type,
                     custom_ensemble_weights = custom_ensemble_weights,
                     compute_insample_predictions =
                       compute_insample_predictions,
                     compute_predictions_bylearner =
                       compute_predictions_bylearner,
                     subsamples = subsamples,
                     cv_subsamples_list = cv_subsamples_list,
                     silent = silent, progress = progress,
                     auxiliary_X = auxiliary_X)
  }#IFELSE
  update_progress(silent)

  # Return estimates
  return(res)
}#GET_CEF

# Utility to print progress to console
update_progress <- function(silent) {
  if (!silent) cat(" -- Done! \n")
}#UPDATE_PROGRESS

# Construct CEF from auxiliary_X
extrapolate_CEF <- function(D, CEF_res_byD, aux_indxs) {
  # Data parameters
  nCEF <- length(CEF_res_byD)
  nobs <- length(D)
  D_levels <- lapply(CEF_res_byD, function(x) x$d)
  is_D <- rep(list(NULL), nCEF)
  for (d in 1:nCEF) is_D[[d]] <- which(D == D_levels[d])
  nensb <- ncol(as.matrix(CEF_res_byD[[1]][[1]]$oos_fitted))
  sample_folds <- length(CEF_res_byD[[1]][[1]]$auxiliary_fitted)

  # Populate CEF
  g_X_byD <- array(0, dim = c(nobs, nensb, nCEF))
  for (d in 1:nCEF) {
    g_X_byD[is_D[[d]], , d] <- CEF_res_byD[[d]][[1]]$oos_fitted
    for (k in 1:sample_folds) {
      g_X_byD[aux_indxs[[d]][[k]], , d] <-
        CEF_res_byD[[d]][[1]]$auxiliary_fitted[[k]]
    }#FOR
  }#FOR

  # return as array, third dimension is different levels of d
  g_X_byD
}#EXTRAPOLATE_CEF
