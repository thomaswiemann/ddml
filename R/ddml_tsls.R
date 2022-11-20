#' Compute DDML IV estimators.
#'
#' Compute DDML IV estimators.
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
#' @param models_FS Same as \code{models}. May be used to consider an seperate
#'     set of models to be used for the estimation of E[D|X,Z] (first stage).
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
#'     \code{ddml_iv}.
#'
#' An object of class \code{ddml_iv}is a list containig the
#'     following components:
#' \describe{
#' \item{\code{coef}}{A vector with the LIE-conform DDML IV coefficent in the
#'     first entry.}
#'     \item{\code{...}}{...}
#' }
#'
#'
#' @export ddml_tsls
ddml_tsls <- function(y, D, Z, X = matrix(1, nobs(y)),
                    models,
                    models_DX = models,
                    models_ZX = models,
                    ens_type = c("average"),
                    cv_folds = 5,
                    sample_folds = 2,
                    subsamples = NULL,
                    setup_parallel = list(type = 'dynamic', cores = 1),
                    silent = F) {
  # Data parameters
  nobs <- length(y)
  nmodels <- length(models)

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
                       sample_folds, subsamples, compute_is_predictions = F,
                       setup_parallel, silent)

  # Compute estimates of E[D|X].
  D_X_res <- crosspred(D, X, Z = NULL,
                       models_DX, ens_type, cv_folds,
                       sample_folds, subsamples,
                       compute_is_predictions = F,
                       setup_parallel, silent)

  # Compute estimates of E[Z|X].
  Z_X_res <- crosspred(Z, X, Z = NULL,
                       models_ZX, ens_type, cv_folds,
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
    V_r <- Z - Z_X_res$oos_fitted

    # Compute IV estimate with constructed variables
    iv_fit <- tsls(y_r, D_r, V_r)

    # Organize complementary ensemble output
    coef <- iv_fit$coef[1]
    weights <- list(Z_X = Z_X_res$weights,
                    D_X = D_X_res$weights,
                    y_X = y_X_res$weights)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (sim_ens) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    nensb <- length(ens_type)
    coef <- matrix(0, 1, nensb)
    mspe <- iv_fit <- rep(list(1), nensb)
    weights <- replicate(3, array(0, dim = c(nmodels, nensb, sample_folds)),
                         simplify = F)
    weights[[1]] <- Z_X_res$weights; weights[[3]] <- y_X_res$weights
    weights[[2]] <- D_X_res$weights
    names(weights) <- c("Z_X", "D_X", "y_X")
    for (j in 1:nensb) {
      # Residualize
      y_r <- y - y_X_res$oos_fitted[, j]
      D_r <- D - D_X_res$oos_fitted[, j]
      V_r <- Z - Z_X_res$oos_fitted[, j]

      # Compute IV estimate with constructed variables
      iv_fit_j <- tsls(y_r, D_r, V_r)

      # Organize complementary ensemble output
      coef[j] <- iv_fit_j$coef[1]
      iv_fit[[j]] <- iv_fit_j
    }#FOR
    # Name output appropriately by ensemble type
    names(iv_fit) <- ens_type
  }#IF

  # Store complementary ensemble output
  mspe <- list(Z_X = Z_X_res$mspe,
               D_X = D_X_res$mspe,
               y_X = y_X_res$mspe)

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights,
                   mspe = mspe,
                   models = models,
                   models_ZX = models_ZX, models_DX = models_DX,
                   iv_fit = iv_fit,
                   subsamples = subsamples,
                   ens_type = ens_type,
                   nobs = nobs, y = y, D = D)

  # Amend class and return
  class(ddml_fit) <- c("ddml_tsls")
  return(ddml_fit)
}#DDML_TSLS
