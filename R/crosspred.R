#' Cross-Predictions using Stacking.
#'
#' @family utilities
#'
#' @description Cross-predictions using stacking.
#'
#' @inheritParams crossval
#' @param learners May take one of two forms, depending on whether a single
#'     learner or stacking with multiple learners is used for estimation of the
#'     predictor.
#'     If a single learner is used, \code{learners} is a list with two named
#'     elements:
#'     \itemize{
#'         \item{\code{what} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}.}
#'     }
#'     If stacking with multiple learners is used, \code{learners} is a list of
#'     lists, each containing four named elements:
#'     \itemize{
#'         \item{\code{fun} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}.}
#'         \item{\code{assign_X} An optional vector of column indices
#'             corresponding to predictive variables in \code{X} that are passed to
#'             the base learner.}
#'         \item{\code{assign_Z} An optional vector of column indices
#'             corresponding to predictive in \code{Z} that are passed to the
#'             base learner.}
#'     }
#'     Omission of the \code{args} element results in default arguments being
#'     used in \code{fun}. Omission of \code{assign_X} (and/or \code{assign_Z})
#'     results in inclusion of all variables in \code{X} (and/or \code{Z}).
#' @param sample_folds Number of cross-fitting folds.
#' @param ensemble_type Ensemble method to combine base learners into final
#'     estimate of the conditional expectation functions. Possible values are:
#'     \itemize{
#'         \item{\code{"nnls"} Non-negative least squares.}
#'         \item{\code{"nnls1"} Non-negative least squares with the constraint
#'             that all weights sum to one.}
#'         \item{\code{"singlebest"} Select base learner with minimum MSPE.}
#'         \item{\code{"ols"} Ordinary least squares.}
#'         \item{\code{"average"} Simple average over base learners.}
#'     }
#'     Multiple ensemble types may be passed as a vector of strings.
#' @param cv_folds Number of folds used for cross-validation in ensemble
#'     construction.
#' @param custom_ensemble_weights A numerical matrix with user-specified
#'     ensemble weights. Each column corresponds to a custom ensemble
#'     specification, each row corresponds to a base learner in \code{learners}
#'     (in chronological order). Optional column names are used to name the
#'     estimation results corresponding the custom ensemble specification.
#' @param compute_insample_predictions Indicator equal to 1 if in-sample
#'     predictions should also be computed.
#' @param compute_predictions_bylearner Indicator equal to 1 if in-sample
#'     predictions should also be computed for each learner (rather than the
#'     entire ensemble).
#' @param subsamples List of vectors with sample indices for cross-fitting.
#' @param cv_subsamples_list List of lists, each corresponding to a subsample
#'     containing vectors with subsample indices for cross-validation.
#' @param auxiliary_X An optional list of matrices of length
#'     \code{sample_folds}, each containing additional observations to calculate
#'     predictions for.
#'
#' @return \code{crosspred} returns a list containing the following components:
#'     \describe{
#'         \item{\code{oos_fitted}}{A matrix of out-of-sample predictions,
#'             each column corresponding to an ensemble type (in chronological
#'             order).}
#'         \item{\code{weights}}{An array, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedures.}
#'         \item{\code{is_fitted}}{When \code{compute_insample_predictions = T}.
#'             a list of matrices with in-sample predictions by sample fold.}
#'         \item{\code{auxiliary_fitted}}{When \code{auxiliary_X} is not
#'             \code{NULL}, a list of matrices with additional predictions.}
#'         \item{\code{oos_fitted_bylearner}}{When
#'             \code{compute_predictions_bylearner = T}, a matrix of
#'             out-of-sample predictions, each column corresponding to a base
#'             learner (in chronological order).}
#'         \item{\code{is_fitted_bylearner}}{When
#'             \code{compute_insample_predictions = T} and
#'             \code{compute_predictions_bylearner = T}, a list of matrices with
#'             in-sample predictions by sample fold.}
#'         \item{\code{auxiliary_fitted_bylearner}}{When \code{auxiliary_X} is
#'             not \code{NULL} and \code{compute_predictions_bylearner = T}, a
#'             list of matrices with additional predictions for each learner.}
#'     }
#' @export
#'
#' @references
#' Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). "Model Averaging and 
#'    Double Machine Learning." Journal of Applied Econometrics, 40(3): 249-269.
#'
#' Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2), 241-259.
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' X = AE98[, c("morekids", "age","agefst","black","hisp","othrace","educ")]
#'
#' # Compute cross-predictions using stacking with base learners ols and lasso.
#' #     Two stacking approaches are simultaneously computed: Equally
#' #     weighted (ensemble_type = "average") and MSPE-minimizing with weights
#' #     in the unit simplex (ensemble_type = "nnls1"). Predictions for each
#' #     learner are also calculated.
#' crosspred_res <- crosspred(y, X,
#'                            learners = list(list(fun = ols),
#'                                            list(fun = mdl_glmnet)),
#'                            ensemble_type = c("average",
#'                                              "nnls1",
#'                                              "singlebest"),
#'                            compute_predictions_bylearner = TRUE,
#'                            sample_folds = 2,
#'                            cv_folds = 2,
#'                            silent = TRUE)
#' dim(crosspred_res$oos_fitted) # = length(y) by length(ensemble_type)
#' dim(crosspred_res$oos_fitted_bylearner) # = length(y) by length(learners)
crosspred <- function(y, X, Z = NULL,
                      learners,
                      sample_folds = 2,
                      ensemble_type = "average",
                      cv_folds = 5,
                      custom_ensemble_weights = NULL,
                      compute_insample_predictions = FALSE,
                      compute_predictions_bylearner = FALSE,
                      subsamples = NULL,
                      cv_subsamples_list = NULL,
                      silent = FALSE,
                      progress = NULL,
                      auxiliary_X = NULL) {
  # Data parameters
  nobs <- nrow(X)
  nlearners <- length(learners)
  calc_ensemble <- !("what" %in% names(learners))
  ncustom <- ncol(custom_ensemble_weights)
  ncustom <- ifelse(is.null(ncustom), 0, ncustom)
  nensb <- length(ensemble_type) + ncustom
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
  oos_fitted <- matrix(0, nobs, nensb^(calc_ensemble))
  oos_fitted_bylearner <- matrix(0, nobs, nlearners)
  is_fitted <- rep(list(NULL), sample_folds)
  is_fitted_bylearner <- rep(list(NULL), sample_folds)
  auxiliary_fitted <- rep(list(NULL), sample_folds)
  auxiliary_fitted_bylearner <- rep(list(NULL), sample_folds)
  mspe <- matrix(0, nlearners^(calc_ensemble), sample_folds)
  colnames(mspe) <- paste("sample fold ", 1:sample_folds)
  weights <- array(0, dim = c(nlearners, nensb, sample_folds))
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
                          custom_weights = custom_ensemble_weights,
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
        mspe[,k] <- mdl_fit$cv_res$mspe
      }#IF
    }#IFELSE
    # Assign names to weights
    dimnames(weights) <- list(NULL, colnames(mdl_fit$weights),
                              paste("sample fold ", 1:sample_folds))
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
    if (!is.null(auxiliary_X)) {
      auxiliary_fitted[[k]] <- stats::predict(mdl_fit,
                                               auxiliary_X[[k]])
    }#if
    # Compute out-of-sample predictions for each learner (optional)
    if (compute_predictions_bylearner) {
      # Adjust ensemble weights
      mdl_fit$weights <- diag(1, nlearners)
      oos_fitted_bylearner[subsamples[[k]], ] <-
        as.numeric(predict.ensemble(mdl_fit,
                                    newdata = X[subsamples[[k]], , drop = F],
                                    newZ = Z[subsamples[[k]], , drop = F]))
      # Compute in-sample predictions (optional)
      if (compute_insample_predictions) {
        is_fitted_bylearner[[k]] <-
          predict.ensemble(mdl_fit, newdata = X[-subsamples[[k]], ,drop = F],
                           newZ = Z[-subsamples[[k]], , drop = F])
      }#IF
      # Compute auxilliary predictions by learner (optional)
      if (!is.null(auxiliary_X)) {
        auxiliary_fitted_bylearner[[k]] <- stats::predict(mdl_fit,
                                                           auxiliary_X[[k]])
      }#if
    }#IF
  }#FOR
  # When multiple ensembles are computed, need to reorganize is_fitted
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
                 is_fitted = is_fitted,
                 auxiliary_fitted = auxiliary_fitted,
                 oos_fitted_bylearner = oos_fitted_bylearner,
                 is_fitted_bylearner = is_fitted_bylearner,
                 auxiliary_fitted_bylearner = auxiliary_fitted_bylearner)
  return(output)
}#CROSSPRED
