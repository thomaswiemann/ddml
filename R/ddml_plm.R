#' Compute DDML PLM estimators.
#'
#' Compute DDML PLM estimators.
#'
#' @param y A response vector.
#' @param D An endogeneous variable vector.
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
#' @return \code{ddml_plm} returns an object of S3 class
#'     \code{ddml_plm}.
#'
#' An object of class \code{ddml_plm}is a list containig the
#'     following components:
#' \describe{
#' \item{\code{coef}}{A vector with the DDML PLM coefficent in the
#'     first entry.}
#'     \item{\code{...}}{...}
#' }
#'
#' @examples
#' # Add example here.
#'
#' @export ddml_plm
ddml_plm <- function(y, D, X,
                    models,
                    ens_type = c("average"),
                    cv_folds = 5,
                    sample_folds = 2,
                    subsamples = NULL,
                    setup_parallel = list(type = 'dynamic', cores = 1),
                    silent = F) {
  # Data parameters
  nobs <- length(y)
  nmodels <- length(models); nensb <- length(ens_type)

  # Draw samples if not user-supplied
  if (is.null(subsamples)) {
    subsamples <- split(c(1:nobs),
                        sample(rep(c(1:sample_folds),
                                   ceiling(nobs / sample_folds)))[1:nobs])
  }#IF
  sample_folds <- length(subsamples)

  # Compute estimates of E[y|X]
  y_X_res <- crosspred(y, X, Z = NULL,
                       models, ens_type, cv_folds,
                       sample_folds, subsamples,
                       compute_is_predictions = F,
                       setup_parallel, silent)

  # Compute estimates of E[D|X].
  D_X_res <- crosspred(D, X, Z = NULL,
                       models, ens_type, cv_folds,
                       sample_folds, subsamples,
                       compute_is_predictions = F,
                       setup_parallel, silent)

  # Check whether multiple ensembles are computed simultaneously
  sim_ens <- length(ens_type) > 1

  # If a single ensemble is calculated, no loops are required.
  if (!sim_ens) {

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - D_X_res$oos_fitted

    # Compute IV estimate with constructed variables
    ols_fit <- ols(y_r, D_r)

    # Organize complementary ensemble output
    coef <- ols_fit$coef[1]
    weights <- list(y_X = y_X_res$weights,
                    D_X = D_X_res$weights)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (sim_ens) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    nensb <- length(ens_type)
    coef <- matrix(0, 1, nensb)
    mspe <- ols_fit <- rep(list(1), nensb)
    weights <- replicate(2, array(0, dim = c(nmodels, nensb, sample_folds)),
                         simplify = F)
    weights[[2]] <- D_X_res$weights; weights[[1]] <- y_X_res$weights
    names(weights) <- c("y_X", "D_X")
    for (j in 1:nensb) {
      # Residualize
      D_r <- D - D_X_res$oos_fitted[, j]

      # Residualize y
      y_r <- y - y_X_res$oos_fitted[, j]

      # Compute IV estimate with constructed variables
      ols_fit_j <- ols(y_r, D_r)

      # Organize complementary ensemble output
      coef[j] <- ols_fit_j$coef[1]
      ols_fit[[j]] <- ols_fit_j
      weights[[2]][, j, ] <- D_X_res$weights[, j, ]
    }#FOR
    # Name output appropriately by ensemble type
    names(ols_fit) <- ens_type
  }#IF

  # Store complementary ensemble output
  mspe <- list(D_X = D_X_res$mspe,
               y_X = y_X_res$mspe)

  # Organize output
  ddml_fit <- list(coef = coef,
                   weights = weights,
                   mspe = mspe,
                   models = models,
                   ols_fit = ols_fit,
                   subsamples = subsamples,
                   ens_type = ens_type,
                   nobs = nobs, y = y, D = D)

  # Amend class and return
  class(ddml_fit) <- c("ddml_plm")
  return(ddml_fit)
}#DDML_PLM
