#' Compute in-sample predictions.
#'
#' Compute in-sample predictions.
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
#' @param setup_parallel An list containing details on the parallelization of
#'    \code{crossval}.
#' @param silent A boolean indicating whether current models and folds should be
#'     printed to the console.
#'
#' @return \code{crosspred} returns a list containig the following
#'     components:
#' \describe{
#' \item{\code{is_fitted}}{A list of matrices of out-of-sample predictions, each
#'     list element corresponding to a training sample, each column
#'     corresponding to an ensemble as passed via \code{ens_type}.}
#' }
#'
#' @export ispred
ispred <- function(y, X, Z = NULL,
                   models,
                   ens_type = c("average"),
                   cv_folds = 5,
                   setup_parallel = list(type = 'dynamic', cores = 1),
                   silent = F) {
  # Data parameters
  nobs <- nrow(X)
  nmodels <- length(models)
  calc_ensemble <- !("what" %in% names(models))
  nensb <- length(ens_type)
  # If instruments are included, check whether columns names are assigned. If
  #    not, assign artificial column names. This is required for checking
  #    whether any IVs were inlcuded in the estimation.
  if(!is.null(Z) & is.null(colnames(Z))) colnames(Z) <- c(1:ncol(Z))
  # Initialize output matrices
  is_fitted <- matrix(0, nobs, nensb^(calc_ensemble))
  mspe <- anyiv_cv <- anyiv <- matrix(0, nmodels^(calc_ensemble), 1)
  weights <- matrix(0, nmodels, nensb)
  # Compute fit on training data. Check whether a single model or an ensemble
  #     should be computed. Check whether the user-supplied response is
  #     training-sample specific.
  if (!calc_ensemble) {
    # When a single model should be fitted, call the constructor function.
    #     Begin with assigning features and response to model arguments.
    #     Note: this is effectively copying the data -- fix needed.
    models$args$X <- cbind(X, Z)
    models$args$y <- y
    # Compute model
    mdl_fit <- do.call(do.call, models)
    # Compute in-sample predictions
    is_fitted <- predict(mdl_fit, cbind(X, Z))
    # Check whether instruments were selected (optional).
    if (!is.null(Z)) {
      index_iv <- (ncol(X) + 1):length(c(ncol(X), ncol(Z)))
      anyiv[1] <- any_iv(obj = mdl_fit,
                         index_iv = index_iv,
                         names_iv = colnames(Z))
    }#IF
  } else if (calc_ensemble) {
    # When multiple models are passed, fit an ensemble on the training data.
    mdl_fit <- ensemble(y, X, Z,
                        ens_type, models, cv_folds,
                        setup_parallel = setup_parallel,
                        silent = silent)
    # Compute in-sample predictions
    is_fitted <- predict(mdl_fit, newX = X, newZ = Z)
    # Record ensemble weights
    weights <- mdl_fit$weights
    # Record model MSPEs
    if (!is.null(mdl_fit$cv_res)) mspe <- colSums(mdl_fit$cv_res$oos_resid^2)
    # Record which models select IVs
    if (!is.null(mdl_fit$cv_res)) anyiv_cv[mdl_fit$mdl_w_iv] <- 1
    # Check whether instruments were selected (optional).
    if (!is.null(Z)) {
      for (m in 1:nmodels) {
        # Check for X, Z assignment
        assign_X <- mdl_fit$models[[m]]$assign_X
        assign_Z <- mdl_fit$models[[m]]$assign_Z
        index_iv <- (length(assign_X) + 1):length(c(assign_X, assign_Z))
        anyiv[m] <- any_iv(obj = mdl_fit$mdl_fits[[m]],
                            index_iv = index_iv,
                            names_iv = colnames(Z[, assign_Z, drop = F]))


      }#FOR
    }#IF
  }#IFELSE
  # Organize and return output
  if (!calc_ensemble) weights <- mspe <- anyiv_cv <- NULL
  output <- list(is_fitted = is_fitted, weights = weights,
                 mspe = mspe, anyiv_cv = anyiv_cv, anyiv = anyiv)
  return(output)
}#ISPRED



