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
                    auxilliary_X = NULL,
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
                         auxilliary_X = auxilliary_X,
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
                     auxilliary_X = auxilliary_X)
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
  nobs <- length(D)
  D_levels <- sort(unique(D))
  n_D_levels <- length(D_levels)
  is_D <- rep(list(NULL), n_D_levels)
  for (d in 1:n_D_levels) is_D[[d]] <- which(D == D_levels[d])
  nensb <- ncol(as.matrix(CEF_res_byD[[1]]$oos_fitted))
  sample_folds <- length(CEF_res_byD[[1]]$auxilliary_fitted)

  # Populate CEF
  g_X_byD <- array(0, dim = c(nobs, nensb, n_D_levels))
  for (d in 1:n_D_levels) {
    g_X_byD[is_D[[d]], , d] <- CEF_res_byD[[d]]$oos_fitted
    for (k in 1:sample_folds) {
      g_X_byD[aux_indxs[[d]][[k]], , d] <-
        CEF_res_byD[[d]]$auxilliary_fitted[[k]]
    }#FOR
  }#FOR

  # return as array, third dimension is different levels of d
  g_X_byD
}#EXTRAPOLATE_CEF
