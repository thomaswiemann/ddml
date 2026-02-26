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
#'             corresponding to predictive variables in \code{X} that are
#'             passed to the base learner.}
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
#' @param cv_subsamples List of lists, each corresponding to a subsample
#'     containing vectors with subsample indices for cross-validation.
#' @param cv_subsamples_list Deprecated; use \code{cv_subsamples} instead.
#' @param auxiliary_X An optional list of matrices of length
#'     \code{sample_folds}, each containing additional observations to calculate
#'     predictions for.
#' @param parallel An optional named list with parallel processing
#'     options. When \code{NULL} (the default), computation is
#'     sequential. Supported fields:
#'     \describe{
#'         \item{\code{cores}}{Number of cores to use.}
#'         \item{\code{export}}{Character vector of object names to
#'             export to parallel workers (for custom learners that
#'             reference global objects).}
#'         \item{\code{packages}}{Character vector of additional
#'             package names to load on workers (for custom learners
#'             that use packages not imported by \code{ddml}).}
#'     }
#'
#' @return \code{crosspred} returns a list containing the following components:
#'     \describe{
#'         \item{\code{oos_fitted}}{A matrix of out-of-sample predictions,
#'             each column corresponding to an ensemble type (in chronological
#'             order).}
#'         \item{\code{weights}}{An array, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedures.}
#'         \item{\code{is_fitted}}{When
#'             \code{compute_insample_predictions = TRUE}.
#'             a list of matrices with in-sample predictions by sample fold.}
#'         \item{\code{auxiliary_fitted}}{When \code{auxiliary_X} is not
#'             \code{NULL}, a list of matrices with additional predictions.}
#'         \item{\code{oos_fitted_bylearner}}{When
#'             \code{compute_predictions_bylearner = TRUE}, a matrix of
#'             out-of-sample predictions, each column corresponding to a base
#'             learner (in chronological order).}
#'         \item{\code{is_fitted_bylearner}}{When
#'             \code{compute_insample_predictions = TRUE} and
#'             \code{compute_predictions_bylearner = TRUE}, a list of matrices
#'             with in-sample predictions by sample fold.}
#'         \item{\code{auxiliary_fitted_bylearner}}{When \code{auxiliary_X}
#'             is not \code{NULL} and
#'             \code{compute_predictions_bylearner = TRUE}, a list of
#'             matrices with additional predictions for each learner.}
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
#'                            compute_predictions_bylearner = TRUE,
#'                            sample_folds = 2,
#'                            cv_folds = 2,
#'                            silent = TRUE)
#' dim(crosspred_res$oos_fitted) # = length(y) by length(ensemble_type)
#' dim(crosspred_res$oos_fitted_bylearner) # = length(y) by length(learners)
crosspred <- function(y, X, Z = NULL,
                      learners,
                      sample_folds = 10,
                      ensemble_type = "average",
                      cv_folds = 10,
                      custom_ensemble_weights = NULL,
                      compute_insample_predictions = FALSE,
                      compute_predictions_bylearner = FALSE,
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
      stop("Specify cv_subsamples or cv_subsamples_list, not both.")
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
    (class(learners[[1]]) != "function")

  # Create crossfitting and cv tuples
  indxs <- get_sample_splits(cluster_variable = cluster_variable,
                             sample_folds = sample_folds,
                             cv_folds = if (w_cv) cv_folds,
                             subsamples = subsamples,
                             cv_subsamples = cv_subsamples)
  subsamples <- indxs$subsamples
  cv_subsamples <- indxs$cv_subsamples
  sample_folds <- length(subsamples)
  cv_folds <- if (!is.null(cv_subsamples)) length(cv_subsamples[[1]]) else 0L

  # Dispatch fold computation
  fold_fun <- function(k) {
    crosspred_compute_fold(
      k = k, y = y, X = X, Z = Z,
      learners = learners,
      subsamples = subsamples,
      cv_subsamples_k = cv_subsamples[[k]],
      ensemble_type = ensemble_type,
      cv_folds = cv_folds,
      custom_ensemble_weights = custom_ensemble_weights,
      calc_ensemble = calc_ensemble,
      nensb = nensb, nlearners = nlearners,
      compute_insample_predictions = compute_insample_predictions,
      compute_predictions_bylearner = compute_predictions_bylearner,
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
  oos_fitted <- matrix(0, nobs, nensb^(calc_ensemble))
  oos_fitted_bylearner <- matrix(0, nobs, nlearners)
  is_fitted <- rep(list(NULL), sample_folds)
  is_fitted_bylearner <- rep(list(NULL), sample_folds)
  auxiliary_fitted <- rep(list(NULL), sample_folds)
  auxiliary_fitted_bylearner <- rep(list(NULL), sample_folds)
  mspe <- matrix(0, nlearners^(calc_ensemble), sample_folds)
  r2 <- matrix(NA_real_, nlearners^(calc_ensemble), sample_folds)
  colnames(mspe) <- colnames(r2) <-
    paste("sample fold ", seq_len(sample_folds))
  weights <- array(0, dim = c(nlearners, nensb, sample_folds))

  for (res in fold_results) {
    k <- res$k
    oos_fitted[res$oos_indices, ] <- res$oos_fitted_rows
    if (!is.null(res$weights_k)) weights[, , k] <- res$weights_k
    if (!is.null(res$mspe_k)) mspe[, k] <- res$mspe_k
    if (!is.null(res$r2_k)) r2[, k] <- res$r2_k
    is_fitted[[k]] <- res$is_fitted_k
    auxiliary_fitted[[k]] <- res$auxiliary_fitted_k
    if (!is.null(res$oos_fitted_bylearner_rows)) {
      oos_fitted_bylearner[res$oos_indices, ] <-
        res$oos_fitted_bylearner_rows
    }#IF
    is_fitted_bylearner[[k]] <- res$is_fitted_bylearner_k
    auxiliary_fitted_bylearner[[k]] <- res$auxiliary_fitted_bylearner_k
  }#FOR

  # Assign dimnames to weights
  if (calc_ensemble) {
    wnames <- fold_results[[1]]$weight_colnames
    dimnames(weights) <- list(
      NULL, wnames,
      paste("sample fold ", seq_len(sample_folds)))
  }#IF

  # Reorganize is_fitted for multiple ensembles
  if (compute_insample_predictions && calc_ensemble && nensb > 1) {
    new_is_fitted <- rep(list(rep(list(1), sample_folds)), nensb)
    for (i in seq_len(nensb)) {
      for (k in seq_len(sample_folds)) {
        new_is_fitted[[i]][[k]] <-
          is_fitted[[k]][, i, drop = FALSE]
      }#FOR
    }#FOR
    is_fitted <- new_is_fitted
  }#IF
  # Compute per-learner OOS residuals (not computed for FPLIV w/ LIE...)
  oos_resid_bylearner <- if (is.numeric(y) && !is.list(y)) {
    drop(y) - oos_fitted_bylearner
  } else {
    NULL
  }#IFELSE

  # Organize and return output
  if (!calc_ensemble) weights <- mspe <- r2 <- NULL
  output <- list(oos_fitted = oos_fitted,
                 weights = weights, mspe = mspe, r2 = r2,
                 is_fitted = is_fitted,
                 auxiliary_fitted = auxiliary_fitted,
                 oos_fitted_bylearner = oos_fitted_bylearner,
                 oos_resid_bylearner = oos_resid_bylearner,
                 is_fitted_bylearner = is_fitted_bylearner,
                 auxiliary_fitted_bylearner = auxiliary_fitted_bylearner)
  return(output)
}#CROSSPRED

crosspred_compute_fold <- function(
    k, y, X, Z, learners, subsamples, cv_subsamples_k,
    ensemble_type, cv_folds, custom_ensemble_weights,
    calc_ensemble, nensb, nlearners,
    compute_insample_predictions, compute_predictions_bylearner,
    auxiliary_X) {

  test_idx <- subsamples[[k]]
  train_idx <- -test_idx

  if (!calc_ensemble) {
    learners$args$X <- cbind(X[train_idx, ], Z[train_idx, ])
    if ("list" %in% class(y)) {
      learners$args$y <- y[[k]]
    } else {
      learners$args$y <- y[train_idx]
    }#IFELSE
    mdl_fit <- do.call(do.call, learners)
    oos_fitted_rows <-
      as.numeric(stats::predict(mdl_fit,
                                cbind(X[test_idx, ],
                                      Z[test_idx, ])))
  } else {
    if ("list" %in% class(y)) {
      y_ <- y[[k]]
    } else {
      y_ <- y[train_idx]
    }#IFELSE
    mdl_fit <- ensemble(y_,
                        X[train_idx, , drop = FALSE],
                        Z[train_idx, , drop = FALSE],
                        ensemble_type, learners,
                        cv_folds, cv_subsamples_k,
                        custom_weights = custom_ensemble_weights,
                        silent = TRUE)
    oos_fitted_rows <-
      as.numeric(predict.ensemble(mdl_fit,
                                  newdata = X[test_idx, , drop = FALSE],
                                  newZ = Z[test_idx, , drop = FALSE]))
  }#IFELSE

  # Ensemble metadata
  weights_k <- if (calc_ensemble) mdl_fit$weights else NULL
  mspe_k <- if (calc_ensemble && !is.null(mdl_fit$cv_res)) {
    mdl_fit$cv_res$mspe
  } else {
    NULL
  }
  r2_k <- if (calc_ensemble && !is.null(mdl_fit$cv_res)) {
    mdl_fit$cv_res$r2
  } else {
    NULL
  }
  weight_colnames <- if (calc_ensemble) {
    colnames(mdl_fit$weights)
  } else {
    NULL
  }

  # In-sample predictions (optional)
  is_fitted_k <- NULL
  if (compute_insample_predictions) {
    if (!calc_ensemble) {
      is_fitted_k <- stats::predict(mdl_fit,
                                    cbind(X[train_idx, ],
                                          Z[train_idx, ]))
    } else {
      is_fitted_k <- predict.ensemble(mdl_fit,
                                      newdata = X[train_idx, , drop = FALSE],
                                      newZ = Z[train_idx, , drop = FALSE])
    }#IFELSE
  }#IF

  # Auxiliary predictions (optional)
  auxiliary_fitted_k <- NULL
  if (!is.null(auxiliary_X)) {
    auxiliary_fitted_k <- stats::predict(mdl_fit, auxiliary_X[[k]])
  }#IF

  # By-learner predictions (optional)
  oos_fitted_bylearner_rows <- NULL
  is_fitted_bylearner_k <- NULL
  auxiliary_fitted_bylearner_k <- NULL
  if (compute_predictions_bylearner) {
    mdl_fit$weights <- diag(1, nlearners)
    oos_fitted_bylearner_rows <-
      as.numeric(predict.ensemble(mdl_fit,
                                  newdata = X[test_idx, , drop = FALSE],
                                  newZ = Z[test_idx, , drop = FALSE]))
    if (compute_insample_predictions) {
      is_fitted_bylearner_k <-
        predict.ensemble(mdl_fit,
                         newdata = X[train_idx, , drop = FALSE],
                         newZ = Z[train_idx, , drop = FALSE])
    }#IF
    if (!is.null(auxiliary_X)) {
      auxiliary_fitted_bylearner_k <-
        stats::predict(mdl_fit, auxiliary_X[[k]])
    }#IF
  }#IF

  list(
    k = k,
    oos_indices = test_idx,
    oos_fitted_rows = oos_fitted_rows,
    weights_k = weights_k,
    mspe_k = mspe_k,
    r2_k = r2_k,
    weight_colnames = weight_colnames,
    is_fitted_k = is_fitted_k,
    auxiliary_fitted_k = auxiliary_fitted_k,
    oos_fitted_bylearner_rows = oos_fitted_bylearner_rows,
    is_fitted_bylearner_k = is_fitted_bylearner_k,
    auxiliary_fitted_bylearner_k = auxiliary_fitted_bylearner_k
  )
}#CROSSPRED_COMPUTE_FOLD
