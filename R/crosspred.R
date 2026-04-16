#' Cross-Fitted Predictions Using Stacking
#'
#' @family utilities
#'
#' @description Cross-fitted predictions using stacking.
#'
#' @details \code{crosspred} implements the cross-fitting step of the
#'     Double/Debiased Machine Learning procedure combined with
#'     stacking. It produces the cross-fitted nuisance estimates
#'     \eqn{\hat{\eta}(X_i)} used in the Neyman orthogonal scores of
#'     all \code{ddml_*} estimators.
#'
#' Let \eqn{\{I_1, \ldots, I_S\}} be an \eqn{S}-fold partition of
#'     \eqn{\{1, \ldots, n\}}, and denote the training set for fold
#'     \eqn{s} by
#'     \eqn{\mathcal{T}_s = \{1, \ldots, n\} \setminus I_s}.
#'     Given \eqn{J} base learners, the procedure operates on each
#'     cross-fitting fold \eqn{s} in three steps:
#'
#' \strong{Step 1 (Stacking weights).}
#'     Run \eqn{K}-fold cross-validation on \eqn{\mathcal{T}_s}
#'     (via \code{\link{crossval}}) to estimate the MSPE of each
#'     base learner, and solve for fold-specific stacking weights
#'     \eqn{\hat{w}_s = (\hat{w}_{1,s}, \ldots, \hat{w}_{J,s})'}.
#'
#' \strong{Step 2 (Fit).}
#'     Fit each base learner \eqn{j} on the full training set
#'     \eqn{\mathcal{T}_s}, yielding \eqn{\hat{f}_{j,s}(\cdot)}.
#'
#' \strong{Step 3 (Predict).}
#'     For each \eqn{i \in I_s}, compute the ensemble cross-fitted
#'     prediction
#'
#' \eqn{\hat{\eta}(X_i) = \sum_{j=1}^{J} \hat{w}_{j,s} \hat{f}_{j,s}(X_i).}
#'
#' Since every observation belongs to exactly one fold, the result is
#'     a complete \eqn{n}-vector of out-of-sample predictions.
#'     Crucially, both the stacking weights \eqn{\hat{w}_s} and the
#'     base learner fits \eqn{\hat{f}_{j,s}} depend only on
#'     \eqn{\mathcal{T}_s}, which does not contain observation
#'     \eqn{i}.
#'
#' When a single learner is used (\eqn{J = 1}), no stacking or inner
#'     cross-validation is performed: the learner is simply fitted on
#'     \eqn{\mathcal{T}_s} and predictions are made for \eqn{I_s}.
#'
#' @inheritParams crossval
#' @inheritParams ddml-intro
#' @param subsamples List of vectors with sample indices for cross-fitting.
#' @param cv_subsamples List of lists, each corresponding to a subsample
#'     containing vectors with subsample indices for cross-validation.
#' @param cv_subsamples_list Deprecated; use \code{cv_subsamples} instead.
#' @param auxiliary_X An optional list of matrices of length
#'     \code{sample_folds}, each containing additional observations to calculate
#'     predictions for.
#'
#' @return \code{crosspred} returns a list containing the following components:
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
#'         \item{\code{cv_resid_byfold}}{A list (length \code{sample_folds})
#'             of inner cross-validation residual matrices used for ensemble
#'             weight estimation. \code{NULL} when a single learner is used.}
#'         \item{\code{auxiliary_fitted}}{When \code{auxiliary_X} is not
#'             \code{NULL}, a list of matrices with additional predictions.}
#'         \item{\code{cf_fitted_bylearner}}{A matrix of out-of-sample
#'             predictions, each column corresponding to a base learner
#'             (in chronological order).}
#'         \item{\code{cf_resid_bylearner}}{A matrix of out-of-sample
#'             residuals (\code{y - cf_fitted_bylearner}), each column
#'             corresponding to a base learner.}
#'         \item{\code{auxiliary_fitted_bylearner}}{When \code{auxiliary_X}
#'             is not \code{NULL}, a list of matrices with additional
#'             predictions for each learner.}
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
#'                            learners = list(list(what = ols),
#'                                            list(what = mdl_glmnet)),
#'                            ensemble_type = c("average",
#'                                              "nnls1",
#'                                              "singlebest"),
#'                            sample_folds = 2,
#'                            cv_folds = 2,
#'                            silent = TRUE)
#' dim(crosspred_res$cf_fitted) # = length(y) by length(ensemble_type)
#' dim(crosspred_res$cf_fitted_bylearner) # = length(y) by length(learners)
crosspred <- function(y, X,
                      learners,
                      sample_folds = 10,
                      ensemble_type = "average",
                      cv_folds = 10,
                      custom_ensemble_weights = NULL,
                      cluster_variable = seq_along(y),
                      subsamples = NULL,
                      cv_subsamples = NULL,
                      cv_subsamples_list = NULL,
                      silent = FALSE,
                      auxiliary_X = NULL,
                      parallel = NULL) {
  # Backward compatibility for renamed parameter
  if (!is.null(cv_subsamples_list)) {
    if (!is.null(cv_subsamples))
      stop("Specify cv_subsamples or cv_subsamples_list, not both.",
           call. = FALSE)
    message("Note: cv_subsamples_list has been renamed to cv_subsamples.")
    cv_subsamples <- cv_subsamples_list
  }#IF

  # Data parameters
  nobs <- nrow(X)
  nlearners <- if (is_single_learner(learners)) 1L
    else length(learners)
  ncustom <- n_custom(custom_ensemble_weights)
  nensb <- length(ensemble_type) + ncustom

  # Create crossfitting and cv tuples
  indxs <- get_sample_splits(cluster_variable = cluster_variable,
                             sample_folds = sample_folds,
                             cv_folds = cv_folds,
                             subsamples = subsamples,
                             cv_subsamples = cv_subsamples)
  subsamples <- indxs$subsamples
  cv_subsamples <- indxs$cv_subsamples
  sample_folds <- length(subsamples)
  cv_folds <- if (!is.null(cv_subsamples))
    length(cv_subsamples[[1]]) else 0L
  if (cv_folds == 0L)
    cv_subsamples <- rep(list(NULL), sample_folds)

  # Dispatch fold computation
  fold_fun <- function(k) {
    crosspred_compute_fold(
      k = k, y = y, X = X,
      learners = learners,
      subsamples = subsamples,
      cv_subsamples_k = cv_subsamples[[k]],
      ensemble_type = ensemble_type,
      cv_folds = cv_folds,
      custom_ensemble_weights = custom_ensemble_weights,
      nensb = nensb, nlearners = nlearners,
      auxiliary_X = auxiliary_X)
  }#FOLD_FUN

  fold_results <- with_parallel(sample_folds, fold_fun,
                                parallel, silent)

  # Assemble results from fold_results list
  cf_fitted <- matrix(0, nobs, nensb)
  cf_fitted_bylearner <- matrix(0, nobs, nlearners)
  auxiliary_fitted <- rep(list(NULL), sample_folds)
  auxiliary_fitted_bylearner <- rep(list(NULL), sample_folds)
  weights <- array(0, dim = c(nlearners, nensb, sample_folds))

  cv_resid_byfold <- rep(list(NULL), sample_folds)
  for (res in fold_results) {
    k <- res$k
    cf_fitted[res$test_indices, ] <- res$cf_fitted_rows
    if (!is.null(res$weights_k)) weights[, , k] <- res$weights_k
    cv_resid_byfold[[k]] <- res$cv_resid_byfold_k
    auxiliary_fitted[[k]] <- res$auxiliary_fitted_k
    cf_fitted_bylearner[res$test_indices, ] <-
      res$cf_fitted_bylearner_rows
    auxiliary_fitted_bylearner[[k]] <- res$auxiliary_fitted_bylearner_k
  }#FOR

  # Assign dimnames to weights and cf_fitted
  wnames <- fold_results[[1]]$weight_colnames
  dimnames(weights) <- list(
    NULL, wnames,
    paste("sample fold ", seq_len(sample_folds)))
  colnames(cf_fitted) <- wnames

  # Compute per-learner OOS residuals
  cf_resid_bylearner <- drop(y) - cf_fitted_bylearner

  # Per-learner OOS mspe and r-squared (always available)
  oos_stats <- compute_mspe_r2(cf_resid_bylearner, y)
  mspe <- oos_stats$mspe
  r2 <- oos_stats$r2
  
  if (nlearners > 1) {
    # Ensemble OOS mspe and r-squared
    cf_resid_ens <- drop(y) - cf_fitted
    oos_stats_ens <- compute_mspe_r2(cf_resid_ens, y)
    
    mspe <- c(mspe, oos_stats_ens$mspe)
    r2 <- c(r2, oos_stats_ens$r2)
    names(mspe) <- names(r2) <- c(paste0("learner_", seq_len(nlearners)), 
                                  wnames)
  }#IF

  # Organize and return output
  output <- list(cf_fitted = cf_fitted,
                 weights = weights, mspe = mspe, r2 = r2,
                 cv_resid_byfold = cv_resid_byfold,
                 auxiliary_fitted = auxiliary_fitted,
                 cf_fitted_bylearner = cf_fitted_bylearner,
                 cf_resid_bylearner = cf_resid_bylearner,
                 auxiliary_fitted_bylearner = auxiliary_fitted_bylearner)
  return(output)
}#CROSSPRED

crosspred_compute_fold <- function(
    k, y, X, learners, subsamples, cv_subsamples_k,
    ensemble_type, cv_folds, custom_ensemble_weights,
    nensb, nlearners, auxiliary_X) {

  test_idx <- subsamples[[k]]
  train_idx <- -test_idx

  # Always route through ensemble() (handles J=1 trivially)
  mdl_fit <- ensemble(y[train_idx],
                      X[train_idx, , drop = FALSE],
                      ensemble_type, learners,
                      cv_folds, cv_subsamples_k,
                      custom_weights = custom_ensemble_weights,
                      silent = TRUE)
  cf_fitted_rows <-
    as.numeric(stats::predict(mdl_fit,
                              newdata = X[test_idx,
                                          , drop = FALSE]))

  # Ensemble metadata
  weights_k <- mdl_fit$weights
  cv_resid_byfold_k <- if (!is.null(mdl_fit$cv_results)) {
    mdl_fit$cv_results$cv_resid
  }
  weight_colnames <- colnames(mdl_fit$weights)

  # Auxiliary predictions (optional)
  auxiliary_fitted_k <- NULL
  if (!is.null(auxiliary_X)) {
    auxiliary_fitted_k <- stats::predict(mdl_fit,
                                         auxiliary_X[[k]])
  }#IF

  # By-learner predictions
  cf_fitted_bylearner_rows <- stats::predict(
    mdl_fit, newdata = X[test_idx, , drop = FALSE],
    type = "bylearner")
  auxiliary_fitted_bylearner_k <- NULL
  if (!is.null(auxiliary_X)) {
    auxiliary_fitted_bylearner_k <- stats::predict(
      mdl_fit, auxiliary_X[[k]], type = "bylearner")
  }#IF

  list(
    k = k,
    test_indices = test_idx,
    cf_fitted_rows = cf_fitted_rows,
    weights_k = weights_k,
    cv_resid_byfold_k = cv_resid_byfold_k,
    weight_colnames = weight_colnames,
    auxiliary_fitted_k = auxiliary_fitted_k,
    cf_fitted_bylearner_rows = cf_fitted_bylearner_rows,
    auxiliary_fitted_bylearner_k = auxiliary_fitted_bylearner_k
  )
}#CROSSPRED_COMPUTE_FOLD
