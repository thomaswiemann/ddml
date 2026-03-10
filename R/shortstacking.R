#' Predictions using Short-Stacking
#'
#' @family utilities
#'
#' @description Predictions using short-stacking.
#'
#' @inheritParams crosspred
#' @inheritParams ddml-intro
#'
#' @return \code{shortstack} returns a list containing the following components:
#'     \describe{
#'         \item{\code{cf_fitted}}{A matrix of out-of-sample predictions,
#'             each column corresponding to an ensemble type (in chronological
#'             order).}
#'         \item{\code{weights}}{An array, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedures.}
#'         \item{\code{mspe}}{A numeric vector of per-learner out-of-sample
#'             MSPEs, computed from cross-fitted residuals.}
#'         \item{\code{r2}}{A numeric vector of per-learner out-of-sample
#'             R-squared values.}
#'         \item{\code{auxiliary_fitted}}{When \code{auxiliary_X} is not
#'             \code{NULL}, a list of matrices with additional predictions.}
#'         \item{\code{cf_fitted_bylearner}}{A matrix of out-of-sample
#'             predictions, each column corresponding to a base learner (in
#'             chronological order).}
#'         \item{\code{cf_resid_bylearner}}{A matrix of per-learner
#'             out-of-sample residuals used for weight estimation.}
#'         \item{\code{auxiliary_fitted_bylearner}}{When \code{auxiliary_X} is
#'             not \code{NULL}, a list of matrices with additional predictions
#'             for each learner.}
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
#'                                 learners = list(list(what = ols),
#'                                                 list(what = mdl_glmnet)),
#'                                 ensemble_type = c("average",
#'                                                   "nnls1",
#'                                                   "singlebest"),
#'                                 sample_folds = 2,
#'                                 silent = TRUE)
#' dim(shortstack_res$cf_fitted) # = length(y) by length(ensemble_type)
#' dim(shortstack_res$cf_fitted_bylearner) # = length(y) by length(learners)
shortstacking <- function(y, X,
                          learners,
                          sample_folds = 2,
                          ensemble_type = "average",
                          custom_ensemble_weights = NULL,
                          cluster_variable = seq_along(y),
                          subsamples = NULL,
                          silent = FALSE,
                          auxiliary_X = NULL,
                          parallel = NULL) {
  # Data parameters
  nobs <- nrow(X)
  nlearners <- length(learners)

  if (is_single_learner(learners)) {
    stop("shortstacking cannot be estimated with a single learner.")
  }#IF

  # Create crossfitting tuples
  indxs <- get_sample_splits(cluster_variable,
                             sample_folds = sample_folds,
                             subsamples = subsamples)
  subsamples <- indxs$subsamples
  sample_folds <- length(subsamples)

  # Compute out-of-sample predictions for each learner.
  # Pass cv_folds = NULL to skip unnecessary CV subsample creation
  # since shortstacking computes its own weights from OOS residuals.
  res <- crosspred(y, X,
                   learners = learners,
                   ensemble_type = "average",
                   cv_folds = NULL,
                   subsamples = subsamples,
                   silent = silent,
                   auxiliary_X = auxiliary_X,
                   parallel = parallel)

  # Compute ensemble weights via subsample cross-fitted residual
  fakecv <- list()
  fakecv$cv_resid <- matrix(y, nobs, nlearners) -
    res$cf_fitted_bylearner
  weights <- ensemble_weights(y, X, learners = learners,
                              type = ensemble_type,
                              custom_weights = custom_ensemble_weights,
                              cv_results = fakecv)$weights

  # Compute predictions
  cf_fitted <- res$cf_fitted_bylearner %*% weights

  # Compute auxiliary predictions (optional)
  auxiliary_fitted <- rep(list(NULL), sample_folds)
  if (!is.null(auxiliary_X)) {
    for (k in seq_len(sample_folds)) {
      auxiliary_fitted[[k]] <- res$auxiliary_fitted_bylearner[[k]] %*% weights
    }#FOR
  }#IF

  # Per-learner OOS mspe and r-squared
  cf_resid_bylearner <- as.matrix(fakecv$cv_resid)
  oos_stats <- compute_mspe_r2(cf_resid_bylearner, y)
  mspe <- oos_stats$mspe
  r2 <- oos_stats$r2

  # return shortstacking output
  output <- list(cf_fitted = cf_fitted,
                 weights = weights, mspe = mspe, r2 = r2,
                 auxiliary_fitted = auxiliary_fitted,
                 cf_fitted_bylearner = res$cf_fitted_bylearner,
                 cf_resid_bylearner = cf_resid_bylearner,
                 auxiliary_fitted_bylearner = res$auxiliary_fitted_bylearner)
  return(output)
}#SHORTSTACKING
