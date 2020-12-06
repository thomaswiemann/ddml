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
  calc_ensemble <- !("what" %in% names(models))
  # Draw samples if not user-supplied
  if (is.null(subsamples)) {
    subsamples <- split(c(1:nobs), sample(rep(c(1:sample_folds),
                                ceiling(nobs / sample_folds)))[1:nobs])
  }#IF
  sample_folds <- length(subsamples)
  # Loop over training samples
  oos_fitted <- matrix(0, nobs, length(ens_type)^(calc_ensemble))
  is_fitted <- rep(list(NULL), sample_folds)
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
      oos_fitted[subsamples[[k]],] <- predict(mdl_fit,
                                            cbind(X[subsamples[[k]], ],
                                                  Z[subsamples[[k]], ]))
    } else if (calc_ensemble) {
      # When multiple models are passed, fit an ensemble on the training data.
      if ("list" %in% class(y)) {
        y_ <- y[[k]]
      } else {
        y_ <- y[-subsamples[[k]]]
      }#IFELSE
      mdl_fit <- ensemble(y_, X[-subsamples[[k]], ], Z[-subsamples[[k]], ],
                          ens_type, models, cv_folds,
                          setup_parallel = setup_parallel,
                          silent = silent)
      # Compute out-of-sample predictions
      oos_fitted[subsamples[[k]],] <- predict(mdl_fit,
                                            newX = X[subsamples[[k]], ],
                                            newZ = Z[subsamples[[k]], ])
    }#IFELSE
    # Compute in-sample predictions (optional)
    if (compute_is_predictions) {
      if (!calc_ensemble) {
        is_fitted[[k]] <- predict(mdl_fit, cbind(X[-subsamples[[k]], ],
                                                 Z[-subsamples[[k]], ]))
      } else if (calc_ensemble) {
        is_fitted[[k]] <- predict(mdl_fit, newX = X[-subsamples[[k]], ],
                                  newZ = Z[-subsamples[[k]], ])
      }#IFELSE
    }#IF
  }#FOR
  # Organize and return output
  output <- list(oos_fitted = oos_fitted, is_fitted = is_fitted)
  return(output)
}#CROSSPRED
