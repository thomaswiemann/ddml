#' Title
#'
#' @param y abc
#' @param X abc
#' @param Z abc
#' @param learners abc
#' @param sample_folds abc
#' @param ensemble_type abc
#' @param cv_folds abc
#' @param compute_insample_predictions abc
#' @param subsamples abc
#' @param cv_subsamples_list abc
#' @param silent abc
#' @param progress abc
#' @param auxilliary_X abc
#'
#' @return object
#' @export
#'
#' @examples
#' 1 + 1
crosspred <- function(y, X, Z = NULL,
                      learners,
                      sample_folds = 2,
                      ensemble_type = c("average"),
                      cv_folds = 5,
                      compute_insample_predictions = FALSE,
                      subsamples = NULL,
                      cv_subsamples_list = NULL,
                      silent = F,
                      progress = NULL,
                      auxilliary_X = NULL) {
  # Data parameters
  nobs <- nrow(X)
  nlearners <- length(learners)
  calc_ensemble <- !("what" %in% names(learners))
  # Create sample fold tuple
  if (is.null(subsamples)) {
    subsamples <- generate_subsamples(nobs, sample_folds)
  }#IF
  sample_folds <- length(subsamples)

  # Create cv-subsamples tuple
  if (is.null(cv_subsamples_list)) {
    cv_subsamples_list <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      nobs_k <- nobs - length(subsamples[[k]])
      cv_subsamples_list[[k]] <- generate_subsamples(nobs_k, cv_folds)
    }# FOR
  }#IF
  cv_folds <- length(cv_subsamples_list[[1]])

  # Initialize output matrices
  oos_fitted <- matrix(0, nobs, length(ensemble_type)^(calc_ensemble))
  is_fitted <- rep(list(NULL), sample_folds)
  auxilliary_fitted <- rep(list(NULL), sample_folds)
  mspe <- matrix(0, nlearners^(calc_ensemble), sample_folds)
  colnames(mspe) <- paste("sample fold ", c(1:sample_folds))
  weights <- array(0, dim = c(nlearners, length(ensemble_type), sample_folds))
  dimnames(weights) <- list(NULL, NULL,
                            paste("sample fold ", c(1:sample_folds)))
  # Loop over training samples
  for (k in 1:sample_folds) {
    # Compute fit on training data. Check whether a single model or an ensemble
    #     should be computed. Check whether the user-supplied response is
    #     training-sample specific.
    if (!calc_ensemble) {
      # When a single model should be fitted, call the constructor function.
      #     Begin with assigning features and response to model arguments.
      #     Note: this is effectively copying the data -- improvement needed.
      learners$args$X <- cbind(X[-subsamples[[k]], ],
                               Z[-subsamples[[k]], ])
      if ("list" %in% class(y)) {
        learners$args$y <- y[[k]]
      } else {
        learners$args$y <- y[-subsamples[[k]]]
      }#IFELSE
      # Compute learner
      mdl_fit <- do.call(do.call, learners)
      # Compute out-of-sample predictions
      oos_fitted[subsamples[[k]], ] <-
        as.numeric(stats::predict(mdl_fit, cbind(X[subsamples[[k]], ],
                                          Z[subsamples[[k]], ])))
      # Print progress
      if (!silent) {
        cat(paste0("\r", progress, " sample fold ", k, "/", sample_folds))
      }#IF

    } else if (calc_ensemble) {
      # When multiple learners are passed, fit an ensemble on the training data.
      if ("list" %in% class(y)) {
        y_ <- y[[k]]
      } else {
        y_ <- y[-subsamples[[k]]]
      }#IFELSE

      # Compile progress-preamble
      if (!silent) {
        progress_k = paste0(progress,
                            "sample fold ", k,
                            "/", sample_folds)
        # Print immediately if no cv is needed
        cv_stacking <- c("stacking", "stacking_nn",
                         "stacking_01", "stacking_best")
        if (!any(cv_stacking %in% ensemble_type)) cat(paste0("\r", progress_k))
      }#IF

      # Compute ensemble
      mdl_fit <- ensemble(y_, X[-subsamples[[k]], , drop = F],
                          Z[-subsamples[[k]], , drop = F],
                          ensemble_type, learners,
                          cv_folds, cv_subsamples_list[[k]],
                          silent = silent,
                          progress = paste0(progress_k, ", "))
      # Compute out-of-sample predictions
      oos_fitted[subsamples[[k]], ] <-
        as.numeric(predict.ensemble(mdl_fit,
                           newdata = X[subsamples[[k]], ,
                                    drop = F],
                           newZ = Z[subsamples[[k]], ,
                                    drop = F]))
      # Record ensemble weights
      weights[, , k] <- mdl_fit$weights
      # Record model MSPEs when weights were computed via cross validation
      if (!is.null(mdl_fit$cv_res)) {
        mspe[,k] <- colMeans(mdl_fit$cv_res$oos_resid^2)
      }#IF
    }#IFELSE
    # Compute in-sample predictions (optional)
    if (compute_insample_predictions) {
      if (!calc_ensemble) {
        is_fitted[[k]] <- stats::predict(mdl_fit, cbind(X[-subsamples[[k]], ],
                                                 Z[-subsamples[[k]], ]))
      } else if (calc_ensemble) {
        is_fitted[[k]] <- predict.ensemble(mdl_fit,
                                  newdata = X[-subsamples[[k]], ,drop = F],
                                  newZ = Z[-subsamples[[k]], , drop = F])
      }#IFELSE
    }#IF
    # Compute auxilliary predictions (optional)
    if (!is.null(auxilliary_X)) {
      auxilliary_fitted[[k]] <- stats::predict(mdl_fit,
                                               auxilliary_X[[k]])
    }#if
  }#FOR
  # When multiple ensembles are computed, need to reorganize is_fitted
  nensb <- length(ensemble_type)
  if (compute_insample_predictions & calc_ensemble & nensb > 1) {
    # Loop over each ensemble type to creat list of is_fitted's
    new_is_fitted <- rep(list(rep(list(1), sample_folds)), nensb)
    for (i in 1:nensb) {
      for (k in 1:sample_folds) {
        new_is_fitted[[i]][[k]] <- is_fitted[[k]][, i, drop = F]
      }#FOR
    }#FOR
    is_fitted <- new_is_fitted
  }#IF
  # Organize and return output
  if (!calc_ensemble) weights <- mspe <- NULL
  output <- list(oos_fitted = oos_fitted,
                 weights = weights, mspe = mspe,
                 is_fitted = is_fitted, auxilliary_fitted = auxilliary_fitted)
  return(output)
}#CROSSPRED
