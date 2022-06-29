#' Compute cross-sample predictions.
#'
#' Compute cross-sample predictions.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param Z An optional instrument matrix.
#' @param models May take one of two forms, depending on whether a single model
#'     should be used for prediction, or whether an ensemble procedure should be
#'     employed.
#'     If a single model should be used, \code{models} is a list with two named
#'     elements:
#'     \itemize{
#'         \item{\code{what} The function used to trained the model. The
#'             function must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}.}
#'     }
#'     If an ensemble should be used, \code{models} is a list of lists, each
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
#' @param ens_type A string indicating the type of ensemble. Multiple types may
#'     be passed in form of a vector of strings.
#' @param cv_folds The number for cross-validation folds.
#' @param sample_folds The number of split sample folds used for calculation of
#'     the out-of-sample predictions.
#' @param subsamples An optional list of vectors, each containing indices of
#'     a test-sample. If not used-provided, the split sample folds are randomly
#'     drawn.
#' @param compute_is_predictions A boolean for whether in-sample predictions
#'     should be calculated.
#' @param setup_parallel An list containing details on the parallelization of
#'    \code{crossval}.
#' @param silent A boolean indicating whether current models and folds should be
#'     printed to the console.
#'
#' @return \code{crosspred} returns a list containig the following
#'     components:
#' \describe{
#' \item{\code{oos_fitted}}{A matrix of out-of-sample predictions, each column
#'     corresponding to an ensemble as passed via \code{ens_type}.}
#' \item{\code{is_fitted}}{A list of matrices of out-of-sample predictions, each
#'     list element corresponding to a training sample, each column
#'     corresponding to an ensemble as passed via \code{ens_type}.}
#' }
#'
#' @examples
#' # Add example here.
#'
#' @export crosspred
crosspred <- function(y, X, Z = NULL,
                      models,
                      ens_type = c("average"),
                      cv_folds = 5,
                      sample_folds = 2,
                      subsamples = NULL,
                      compute_is_predictions = FALSE,
                      setup_parallel = list(type = 'dynamic', cores = 1),
                      silent = F) {
  # Data parameters
  nobs <- nrow(X)
  nmodels <- length(models)
  calc_ensemble <- !("what" %in% names(models))
  # Draw samples if not user-supplied
  if (is.null(subsamples)) {
    subsamples <- split(c(1:nobs), sample(rep(c(1:sample_folds),
                                ceiling(nobs / sample_folds)))[1:nobs])
  }#IF
  sample_folds <- length(subsamples)
  # Initialize output matrices
  oos_fitted <- matrix(0, nobs, length(ens_type)^(calc_ensemble))
  is_fitted <- rep(list(NULL), sample_folds)
  mspe <- anyiv_cv <- anyiv <- matrix(0, nmodels^(calc_ensemble), sample_folds)
  weights <- array(0, dim = c(nmodels, length(ens_type), sample_folds))
  # Loop over training samples
  for (k in 1:sample_folds) {
    # Compute fit on training data. Check whether a single model or an ensemble
    #     should be computed. Check whether the user-supplied response is
    #     training-sample specific.
    if (!calc_ensemble) {
      # When a single model should be fitted, call the constructor function.
      #     Begin with assigning features and response to model arguments.
      #     Note: this is effectively copying the data -- fix needed.
      models$args$X <- cbind(X[-subsamples[[k]], ],
                              Z[-subsamples[[k]], ])
      if ("list" %in% class(y)) {
        models$args$y <- y[[k]]
      } else {
        models$args$y <- y[-subsamples[[k]]]
      }#IFELSE
      # Compute model
      mdl_fit <- do.call(do.call, models)
      # Compute out-of-sample predictions
      oos_fitted[subsamples[[k]], ] <-
        as.numeric(predict(mdl_fit, cbind(X[subsamples[[k]], ],
                                          Z[subsamples[[k]], ])))
      # Check whether instruments were selected (optional).
      if (!is.null(Z)) {
        index_iv <- (ncol(X) + 1):length(c(ncol(X), ncol(Z)))
        anyiv[1, k] <- any_iv(obj = mdl_fit,
                              index_iv = index_iv,
                              names_iv = colnames(Z))
      }#IF
    } else if (calc_ensemble) {
      # When multiple models are passed, fit an ensemble on the training data.
      if ("list" %in% class(y)) {
        y_ <- y[[k]]
      } else {
        y_ <- y[-subsamples[[k]]]
      }#IFELSE
      mdl_fit <- ensemble(y_, X[-subsamples[[k]], , drop = F],
                          Z[-subsamples[[k]], , drop = F],
                          ens_type, models, cv_folds,
                          setup_parallel = setup_parallel,
                          silent = silent)
      # Compute out-of-sample predictions
      oos_fitted[subsamples[[k]], ] <-
        as.numeric(predict(mdl_fit,
                           newX = X[subsamples[[k]], ,
                                    drop = F],
                           newZ = Z[subsamples[[k]], ,
                                    drop = F]))
      # Record ensemble weights
      weights[, , k] <- mdl_fit$weights
      # Record model MSPEs
      if (!is.null(mdl_fit$cv_res))mspe[,k]<-colMeans(mdl_fit$cv_res$oos_resid^2)
      # Record which models select IVs
      if (!is.null(mdl_fit$cv_res)) anyiv_cv[mdl_fit$mdl_w_iv, k] <- 1
      # Check whether instruments were selected (optional).
      if (!is.null(Z)) {
        for (m in 1:nmodels) {
          # Check if model is included in ensemble
          if (mdl_fit$weights[m] == 0) next
          # Check for X, Z assignment
          assign_X <- mdl_fit$models[[m]]$assign_X
          assign_Z <- mdl_fit$models[[m]]$assign_Z
          index_iv <- (length(assign_X) + 1):length(c(assign_X, assign_Z))
          anyiv[m, k] <- any_iv(obj = mdl_fit$mdl_fits[[m]],
                                index_iv = index_iv,
                                names_iv = colnames(Z[, assign_Z, drop = F]))
        }#FOR
      }#IF
    }#IFELSE
    # Compute in-sample predictions (optional)
    if (compute_is_predictions) {
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
  nensb <- length(ens_type)
  if (compute_is_predictions & calc_ensemble & nensb > 1) {
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
  if (!calc_ensemble) weights <- mspe <- anyiv_cv <- NULL
  output <- list(oos_fitted = oos_fitted, is_fitted = is_fitted, anyiv = anyiv,
                 weights = weights, mspe = mspe, anyiv_cv = anyiv_cv)
  return(output)
}#CROSSPRED
