#' Compute cross-sample predictions.
#'
#' Compute cross-sample predictions.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param Z An optional instrument matrix.
#' @param learners May take one of two forms, depending on whether a single model
#'     should be used for prediction, or whether an ensemble procedure should be
#'     employed.
#'     If a single model should be used, \code{learners} is a list with two named
#'     elements:
#'     \itemize{
#'         \item{\code{what} The function used to trained the model. The
#'             function must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}.}
#'     }
#'     If an ensemble should be used, \code{learners} is a list of lists, each
#'     containing four named elements:
#'     \itemize{
#'         \item{\code{fun} The function used to trained the model. The
#'             function must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{assign_X} A vector of indices corresponding to features
#'             in \code{X} that should be used for training.}
#'         \item{\code{assign_Z} An optional vector of indices corresponding to
#'             instruments in \code{Z} that should be used for training.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}}
#'     }
#' @param ensemble_type A string indicating the type of ensemble. Multiple types may
#'     be passed in form of a vector of strings.
#' @param cv_folds The number for cross-validation folds.
#' @param sample_folds The number of split sample folds used for calculation of
#'     the out-of-sample predictions.
#' @param subsamples An optional list of vectors, each containing indices of
#'     a test-sample. If not used-provided, the split sample folds are randomly
#'     drawn.
#' @param compute_insample_predictions A boolean for whether in-sample predictions
#'     should be calculated.
#' @param setup_parallel An list containing details on the parallelization of
#'    \code{crossval}.
#' @param silent A boolean indicating whether current learners and folds should be
#'     printed to the console.
#'
#' @return \code{crosspred} returns a list containig the following
#'     components:
#' \describe{
#' \item{\code{oos_fitted}}{A matrix of out-of-sample predictions, each column
#'     corresponding to an ensemble as passed via \code{ensemble_type}.}
#' \item{\code{is_fitted}}{A list of matrices of out-of-sample predictions, each
#'     list element corresponding to a training sample, each column
#'     corresponding to an ensemble as passed via \code{ensemble_type}.}
#' }
#'
#' @export crosspred
crosspred <- function(y, X, Z = NULL,
                      learners,
                      sample_folds = 2,
                      ensemble_type = c("average"),
                      cv_folds = 5,
                      compute_insample_predictions = FALSE,
                      subsamples = NULL,
                      cv_subsamples_list = NULL,
                      silent = F, progress = NULL) {
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
        as.numeric(predict(mdl_fit, cbind(X[subsamples[[k]], ],
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
        as.numeric(predict(mdl_fit,
                           newX = X[subsamples[[k]], ,
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
        is_fitted[[k]] <- predict(mdl_fit, cbind(X[-subsamples[[k]], ],
                                                 Z[-subsamples[[k]], ]))
      } else if (calc_ensemble) {
        is_fitted[[k]] <- predict(mdl_fit, newX = X[-subsamples[[k]],,drop = F],
                                  newZ = Z[-subsamples[[k]], , drop = F])
      }#IFELSE
    }#IF
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
  output <- list(oos_fitted = oos_fitted, is_fitted = is_fitted,
                 weights = weights, mspe = mspe)
  return(output)
}#CROSSPRED
