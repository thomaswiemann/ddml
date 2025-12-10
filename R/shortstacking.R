#' Predictions using Short-Stacking.
#'
#' @family utilities
#'
#' @description Predictions using short-stacking.
#'
#' @inheritParams crosspred
#' @param shortstack_y Optional vector of the outcome variable to form
#'     short-stacking predictions for. Base learners are always trained on
#'     \code{y}.
#'
#' @return \code{shortstack} returns a list containing the following components:
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
#'         \item{\code{oos_fitted_bylearner}}{A matrix of
#'             out-of-sample predictions, each column corresponding to a base
#'             learner (in chronological order).}
#'         \item{\code{is_fitted_bylearner}}{When
#'             \code{compute_insample_predictions = T}, a list of matrices with
#'             in-sample predictions by sample fold.}
#'         \item{\code{auxiliary_fitted_bylearner}}{When \code{auxiliary_X} is
#'             not \code{NULL}, a
#'             list of matrices with additional predictions for each learner.}
#'     }
#'     Note that unlike \code{crosspred}, \code{shortstack} always computes
#'        out-of-sample predictions for each base learner (at no additional
#'        computational cost).
#' @export
#'
#' @references
#' Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). "Model Averaging and 
#'     Double Machine Learning." Journal of Applied Econometrics, 40(3): 249-269.
#'
#' Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2), 241-259.
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' X = AE98[, c("morekids", "age","agefst","black","hisp","othrace","educ")]
#'
#' # Compute predictions using shortstacking with base learners ols and lasso.
#' #     Two stacking approaches are simultaneously computed: Equally
#' #     weighted (ensemble_type = "average") and MSPE-minimizing with weights
#' #     in the unit simplex (ensemble_type = "nnls1"). Predictions for each
#' #     learner are also calculated.
#' shortstack_res <- shortstacking(y, X,
#'                                 learners = list(list(fun = ols),
#'                                                 list(fun = mdl_glmnet)),
#'                                 ensemble_type = c("average",
#'                                                   "nnls1",
#'                                                   "singlebest"),
#'                                 sample_folds = 2,
#'                                 silent = TRUE)
#' dim(shortstack_res$oos_fitted) # = length(y) by length(ensemble_type)
#' dim(shortstack_res$oos_fitted_bylearner) # = length(y) by length(learners)
shortstacking <- function (y, X, Z = NULL,
                           learners,
                           sample_folds = 2,
                           ensemble_type = "average",
                           custom_ensemble_weights = NULL,
                           compute_insample_predictions = FALSE,
                           subsamples = NULL,
                           silent = FALSE,
                           progress = NULL,
                           auxiliary_X = NULL,
                           shortstack_y = y) {

  # Data parameters
  nobs <- nrow(X)
  nlearners <- length(learners)

  # Throw error if no ensemble is estimated
  calc_ensemble <- !("what" %in% names(learners))
  if (!calc_ensemble) {
    stop("shortstacking cannot be estimated with a single learner.")
  }#IF

  # Create sample fold tuple
  if (is.null(subsamples)) {
    subsamples <- generate_subsamples(nobs, sample_folds)
  }#IF
  sample_folds <- length(subsamples)
  ncustom <- ncol(custom_ensemble_weights)
  ncustom <- ifelse(is.null(ncustom), 0, ncustom)
  nensb <- length(ensemble_type) + ncustom

  # Compute out-of-sample predictions for each learner
  res <- crosspred(y, X, Z,
                   learners = learners,
                   ensemble_type = "average",
                   compute_insample_predictions = compute_insample_predictions,
                   compute_predictions_bylearner = TRUE,
                   subsamples = subsamples,
                   silent = silent, progress = progress,
                   auxiliary_X = auxiliary_X)

  # Compute ensemble weights via subsample cross-fitted residual
  fakecv <- list()
  fakecv$oos_resid <- kronecker(shortstack_y, t(rep(1, nlearners))) -
    res$oos_fitted_bylearner
  weights <- ensemble_weights(shortstack_y, X, learners = learners,
                              type = ensemble_type,
                              custom_weights = custom_ensemble_weights,
                              cv_results = fakecv)$weights

  # Compute predictions
  oos_fitted <- res$oos_fitted_bylearner %*% weights

  # Compute auxilliary predictions (optional)
  auxiliary_fitted <- rep(list(NULL), sample_folds)
  if (!is.null(auxiliary_X)) {
    for (k in 1:sample_folds) {
      auxiliary_fitted[[k]] <- res$auxiliary_fitted_bylearner[[k]] %*% weights
    }#FOR
  }#if

  # Compute in-sample predictions (optional)
  is_fitted <- rep(list(NULL), sample_folds)
  fakecv_k <- list()
  if (compute_insample_predictions) {
    for (k in 1:sample_folds) {
      # Compute shortstacking weights in-sample
      fakecv_k$oos_resid <- kronecker(y[-subsamples[[k]]], t(rep(1, nlearners))) -
        res$is_fitted_bylearner[[k]]
      weights_k <- ensemble_weights(y[-subsamples[[k]]], X[-subsamples[[k]], ],
                                    learners = learners,
                                    type = ensemble_type,
                                    custom_weights = custom_ensemble_weights,
                                    cv_results = fakecv_k)$weights
      # Combine base learners
      is_fitted[[k]] <- res$is_fitted_bylearner[[k]] %*% weights_k
    }#FOR

    # When multiple ensembles are computed, need to reorganize is_fitted
    if (nensb > 1) {
      # Loop over each ensemble type to creat list of is_fitted's
      new_is_fitted <- rep(list(rep(list(1), sample_folds)), nensb)
      for (i in 1:nensb) {
        for (k in 1:sample_folds) {
          new_is_fitted[[i]][[k]] <- is_fitted[[k]][, i, drop = F]
        }#FOR
      }#FOR
      is_fitted <- new_is_fitted
    }#IF
  }#IF

  # Compute mspe
  mspe <- colMeans((kronecker(shortstack_y, t(rep(1, nensb))) - oos_fitted)^2)

  # return shortstacking output
  output <- list(oos_fitted = oos_fitted,
                 weights = weights, mspe = mspe,
                 is_fitted = is_fitted,
                 auxiliary_fitted = auxiliary_fitted,
                 oos_fitted_bylearner = res$oos_fitted_bylearner,
                 is_fitted_bylearner = res$is_fitted_bylearner,
                 auxiliary_fitted_bylearner = res$auxiliary_fitted_bylearner)
  return(output)
}#SHORTSTACKING
