# Wrapper for [ddml::crosspred()] and [ddml::shortstacking()].
get_CEF <- function(y, X, Z = NULL,
                    learners,
                    ensemble_type,
                    shortstack,
                    custom_ensemble_weights = NULL,
                    compute_insample_predictions = FALSE,
                    compute_predictions_bylearner = FALSE,
                    subsamples,
                    cv_subsamples,
                    silent = FALSE,
                    label = NULL,
                    auxiliary_X = NULL,
                    shortstack_y = y,
                    parallel = NULL) {
  t0 <- proc.time()[3]

  if (!is.null(label)) {
    info_msg("  Estimating ", label, "...", silent = silent)
  }#IF

  # Compute CEF
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
                     compute_predictions_bylearner =
                       compute_predictions_bylearner,
                     subsamples = subsamples,
                     cv_subsamples = cv_subsamples,
                     silent = silent,
                     auxiliary_X = auxiliary_X,
                     parallel = parallel)
  }#IFELSE

  # Return estimates
  return(res)
}#GET_CEF

# Construct CEF from auxiliary_X
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
