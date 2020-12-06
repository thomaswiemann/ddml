#' Compute LIE-conform DDML IV estimators.
#'
#' Compute LIE-conform DDML IV estimators.
#'
#' @param y A response vector.
#' @param D An endogeneous variable vector.
#' @param Z An instrumental variable vector or matrix.
#' @param X A control matrix.
#' @param models May take one of two forms, depending on whether a single model
#'     should be used for residualization, or whether an ensemble procedure
#'     should be employed.
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
#' @param setup_parallel An list containing details on the parallelization of
#'    \code{crossval}.
#' @param silent A boolean indicating whether current models and folds should be
#'     printed to the console.
#'
#' @return \code{ddml_iv} returns an object of S3 class
#'     \code{c("ddml_iv", "tsls")}.
#'
#' An object of class \code{c("ddml_iv", "tsls")}is a list containig the
#'     following components:
#' @return \code{crosspred} returns a list containig the following
#'     components:
#' \describe{
#' \item{\code{coef}}{A vector with the LIE-conform DDML IV coefficent in the
#'     first entry.}
#' \item{\code{y}}{The residualized response vector.}
#' \item{\code{X_}}{A matrix combining the residualized endogeneous variable D
#'     and a vector of ones.}
#' \item{\code{Z_}}{{A matrix combining a vector of ones and the constructed
#'     instrument estimated from \code{E[D|X,Z]-E[D|X]}.}}
#' \item{\code{FS}}{The first stage coefficient matrix.}
#' }
#'
#' @examples
#' # Add example here.
#'
#' @export ddml_iv
ddml_iv <- function(y, D, Z, X = matrix(1, nobs(y)),
                    models,
                    ens_type = c("average"),
                    cv_folds = 5,
                    sample_folds = 2,
                    subsamples = NULL,
                    setup_parallel = list(type = 'dynamic', cores = 1),
                    silent = F) {
  # Data parameters
  nobs <- length(y)
  # Draw samples if not user-supplied
  if (is.null(subsamples)) {
    subsamples <- split(c(1:nobs), sample(rep(c(1:sample_folds),
                                              ceiling(nobs / sample_folds)))[1:nobs])
  }#IF
  sample_folds <- length(subsamples)
  # Compute estimates of E[D|X,Z]
  D_XZ_res <- crosspred(D, X, Z,
                       models, ens_type, cv_folds,
                       sample_folds, subsamples, compute_is_predictions = T,
                       setup_parallel, silent)
  # Compute LIE-conform estimates of E[D|X]
  D_X_res <- crosspred(D_XZ_res$is_fitted, X, Z = NULL,
                      models, ens_type, cv_folds,
                      sample_folds, subsamples, compute_is_predictions = F,
                      setup_parallel, silent)
  # Compute estimates of E[y|X]
  y_X_res <- crosspred(y, X, Z = NULL,
                      models, ens_type, cv_folds,
                      sample_folds, subsamples, compute_is_predictions = F,
                      setup_parallel, silent)
  # Residualize
  y_r <- y - y_X_res$oos_fitted
  D_r <- D - D_X_res$oos_fitted
  V_r <- D_XZ_res$oos_fitted - D_X_res$oos_fitted
  # Compute IV estimate with constructed variables
  ddml_fit <- tsls(y_r, D_r, V_r)
  # Amend class and return
  class(ddml_fit) <- c("ddml_iv", class(ddml_fit))
  return(ddml_fit)
}#DDML_IV
