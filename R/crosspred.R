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
  # Unpack parallel options
  p <- parse_parallel(parallel)
  num_cores <- p$num_cores
  parallel_export <- p$export
  parallel_packages <- p$packages

  # Backward compatibility for renamed parameter
  if (!is.null(cv_subsamples_list)) {
    if (!is.null(cv_subsamples))
      stop("Specify cv_subsamples or cv_subsamples_list, not both.",
           call. = FALSE)
    message("Note: cv_subsamples_list has been renamed to cv_subsamples.")
    cv_subsamples <- cv_subsamples_list
  }#IF

  # Normalize learner specs before parallel dispatch
  learners <- normalize_learners(learners)

  # Data parameters
  nobs <- nrow(X)
  nlearners <- length(learners)
  calc_ensemble <- !is_single_learner(learners)
  ncustom <- ncol(custom_ensemble_weights)
  ncustom <- ifelse(is.null(ncustom), 0, ncustom)
  nensb <- length(ensemble_type) + ncustom
  # Check whether ddml uses conventional stacking w/ data driven weights
  w_cv <- any(ensemble_type %in% c("nnls", "nnls1", "singlebest", "ols")) &
    (!is.function(learners[[1]]))

  # Create crossfitting and cv tuples
  indxs <- get_sample_splits(cluster_variable = cluster_variable,
                             sample_folds = sample_folds,
                             cv_folds = cv_folds,
                             subsamples = subsamples,
                             cv_subsamples = cv_subsamples)
  subsamples <- indxs$subsamples
  cv_subsamples <- indxs$cv_subsamples
  sample_folds <- length(subsamples)
  cv_folds <- if (!is.null(cv_subsamples)) length(cv_subsamples[[1]]) else 0L

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
      calc_ensemble = calc_ensemble,
      nensb = nensb, nlearners = nlearners,
      auxiliary_X = auxiliary_X)
  }#FOLD_FUN

  cl <- NULL
  if (num_cores > 1) {
    cl <- tryCatch(
      setup_parallel_cluster(num_cores, parallel_export,
                             parallel_packages),
      error = function(e) {
        warning("Parallel setup failed: ",
                conditionMessage(e),
                ". Falling back to sequential.",
                call. = FALSE)
        NULL
      }
    )
    if (!is.null(cl))
      on.exit(parallel::stopCluster(cl), add = TRUE)
  }#IF

  if (silent) {
    op <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(op), add = TRUE)
  }#IF
  fold_results <- pbapply::pblapply(seq_len(sample_folds),
                                    fold_fun, cl = cl)

  # Assemble results from fold_results list
  cf_fitted <- matrix(0, nobs, nensb^(calc_ensemble))
  cf_fitted_bylearner <- matrix(0, nobs, nlearners)
  auxiliary_fitted <- rep(list(NULL), sample_folds)
  auxiliary_fitted_bylearner <- rep(list(NULL), sample_folds)
  mspe_inner <- matrix(0, nlearners^(calc_ensemble), sample_folds)
  colnames(mspe_inner) <- paste("sample fold ", seq_len(sample_folds))
  weights <- array(0, dim = c(nlearners, nensb, sample_folds))

  cv_resid_byfold <- rep(list(NULL), sample_folds)
  for (res in fold_results) {
    k <- res$k
    cf_fitted[res$test_indices, ] <- res$cf_fitted_rows
    if (!is.null(res$weights_k)) weights[, , k] <- res$weights_k
    if (!is.null(res$mspe_k)) mspe_inner[, k] <- res$mspe_k
    cv_resid_byfold[[k]] <- res$cv_resid_byfold_k
    auxiliary_fitted[[k]] <- res$auxiliary_fitted_k
    cf_fitted_bylearner[res$test_indices, ] <-
      res$cf_fitted_bylearner_rows
    auxiliary_fitted_bylearner[[k]] <- res$auxiliary_fitted_bylearner_k
  }#FOR

  # Assign dimnames to weights
  if (calc_ensemble) {
    wnames <- fold_results[[1]]$weight_colnames
    dimnames(weights) <- list(
      NULL, wnames,
      paste("sample fold ", seq_len(sample_folds)))
  }#IF

  # Compute per-learner OOS residuals
  cf_resid_bylearner <- drop(y) - cf_fitted_bylearner

  # Per-learner OOS mspe and r-squared (always available)
  mspe <- colMeans(cf_resid_bylearner^2)
  y_var <- as.numeric(stats::var(y))
  r2 <- if (y_var > 0) 1 - mspe / y_var else
    rep(NA_real_, length(mspe))

  # Organize and return output
  if (!calc_ensemble) weights <- cv_resid_byfold <- NULL
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
    calc_ensemble, nensb, nlearners,
    auxiliary_X) {

  test_idx <- subsamples[[k]]
  train_idx <- -test_idx

  if (!calc_ensemble) {
    learners$args$X <- X[train_idx, ]
    learners$args$y <- y[train_idx]
    mdl_fit <- do.call(do.call, learners)
    cf_fitted_rows <-
      as.numeric(stats::predict(mdl_fit,
                                X[test_idx, ]))
  } else {
    mdl_fit <- ensemble(y[train_idx],
                        X[train_idx, , drop = FALSE],
                        ensemble_type, learners,
                        cv_folds, cv_subsamples_k,
                        custom_weights = custom_ensemble_weights,
                        silent = TRUE)
    cf_fitted_rows <-
      as.numeric(stats::predict(mdl_fit,
                                newdata = X[test_idx, , drop = FALSE]))
  }#IFELSE

  # Ensemble metadata
  weights_k <- if (calc_ensemble) mdl_fit$weights else NULL
  mspe_k <- if (calc_ensemble &&
      !is.null(mdl_fit$cv_results)) {
    mdl_fit$cv_results$mspe
  } else {
    NULL
  }
  r2_k <- if (calc_ensemble &&
      !is.null(mdl_fit$cv_results)) {
    mdl_fit$cv_results$r2
  } else {
    NULL
  }
  cv_resid_byfold_k <- if (calc_ensemble &&
      !is.null(mdl_fit$cv_results)) {
    mdl_fit$cv_results$cv_resid
  } else {
    NULL
  }
  weight_colnames <- if (calc_ensemble) {
    colnames(mdl_fit$weights)
  } else {
    NULL
  }

  # Auxiliary predictions (optional)
  auxiliary_fitted_k <- NULL
  if (!is.null(auxiliary_X)) {
    auxiliary_fitted_k <- stats::predict(mdl_fit, auxiliary_X[[k]])
  }#IF

  # By-learner predictions
  auxiliary_fitted_bylearner_k <- NULL
  if (!calc_ensemble) {
    cf_fitted_bylearner_rows <- matrix(cf_fitted_rows, ncol = 1)
    if (!is.null(auxiliary_X)) {
      auxiliary_fitted_bylearner_k <- matrix(auxiliary_fitted_k,
                                             ncol = 1)
    }#IF
  } else {
    mdl_fit_bylearner <- mdl_fit
    mdl_fit_bylearner$weights <- diag(1, nlearners)
    cf_fitted_bylearner_rows <- stats::predict(
      mdl_fit_bylearner,
      newdata = X[test_idx, , drop = FALSE]
    )
    if (!is.null(auxiliary_X)) {
      auxiliary_fitted_bylearner_k <- stats::predict(
        mdl_fit_bylearner,
        auxiliary_X[[k]]
      )
    }#IF
  }#IFELSE

  list(
    k = k,
    test_indices = test_idx,
    cf_fitted_rows = cf_fitted_rows,
    weights_k = weights_k,
    mspe_k = mspe_k,
    r2_k = r2_k,
    cv_resid_byfold_k = cv_resid_byfold_k,
    weight_colnames = weight_colnames,
    auxiliary_fitted_k = auxiliary_fitted_k,
    cf_fitted_bylearner_rows = cf_fitted_bylearner_rows,
    auxiliary_fitted_bylearner_k = auxiliary_fitted_bylearner_k
  )
}#CROSSPRED_COMPUTE_FOLD
