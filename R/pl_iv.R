#' Compute LIE-conform partially linear IV estimators w/o sample splitting.
#'
#' Compute LIE-conform partially linear IV estimators w/o sample splitting.
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
#' @param setup_parallel An list containing details on the parallelization of
#'    \code{crossval}.
#' @param silent A boolean indicating whether current models and folds should be
#'     printed to the console.
#'
#' @return \code{pl_iv} returns an object of S3 class
#'     \code{pl_iv}.
#'
#' An object of class \code{pl_iv}is a list containig the
#'     following components:
#' \describe{
#' \item{\code{coef}}{A vector with ...}
#'     \item{\code{...}}{...}
#' }
#'
#' @export pl_iv
pl_iv <- function(y, D, Z, X = matrix(1, nobs(y)),
                  models,
                  models_DXZ = models,
                  models_DX = models_DXZ,
                  ens_type = c("average"),
                  cv_folds = 5,
                  enforce_LIE = TRUE,
                  setup_parallel = list(type = 'dynamic', cores = 1),
                  silent = F) {
  # Data parameters
  nobs <- length(y)
  nmodels <- length(models)
  nmodels_DX <- length(models_DX)
  nmodels_DXZ <- length(models_DXZ)
  calc_ensemble <- !("what" %in% names(models))
  calc_ensemble_DX <- !("what" %in% names(models_DX))
  calc_ensemble_DXZ <- !("what" %in% names(models_DXZ))
  nensb <- length(ens_type)

  # Compute estimates of E[y|X]
  y_X_res <- ispred(y, X, Z = NULL,
                    models, ens_type,
                    cv_folds, setup_parallel, silent)

  # Compute estimates of E[D|X,Z]. Also calculate in-sample predictions when
  #     the LIE is enforced.
  D_XZ_res <- ispred(D, X, Z,
                     models_DXZ, ens_type,
                     cv_folds, setup_parallel, silent)

  # When the LIE is not enforced, estimating E[D|X] is straightforward.
  if (!enforce_LIE) {
    D_X_res <- ispred(D, X, Z = NULL,
                      models_DX, ens_type,
                      cv_folds, setup_parallel, silent)
  }#IF

  # Check whether multiple ensembles are computed simultaneously
  sim_ens <- calc_ensemble && nensb > 1
  # If a single ensemble is calculated, no loops are required.
  if (!sim_ens) {
    # Check whether the law of iterated expectations (LIE) should be enforced.
    #     When the LIE is enforced (recommended), the estimates of E[D|X,Z] are
    #     used for the calculation of the estimates of E[D|X].
    if (enforce_LIE) {
      # Compute LIE-conform estimates of E[D|X]
      D_X_res <- ispred(D_XZ_res$is_fitted, X, Z = NULL,
                        models_DX, ens_type,
                        cv_folds, setup_parallel, silent)
    }#IFELSE

    # Residualize
    y_r <- y - y_X_res$is_fitted
    D_r <- D - D_X_res$is_fitted
    V_r <- D_XZ_res$is_fitted - D_X_res$is_fitted

    # Compute IV estimate with constructed variables
    iv_fit <- tsls(y_r, D_r, V_r)

    # Organize complementary ensemble output
    coef <- iv_fit$coef[1]
    weights <- list(D_XZ = D_XZ_res$weights,
                    D_X = D_X_res$weights,
                    y_X = y_X_res$weights)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (sim_ens) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    coef <- matrix(0, 1, nensb)
    mspe <- anyiv_cv <- anyiv <- iv_fit <- rep(list(1), nensb)
    weights <- replicate(3, matrix(0, nmodels, nensb),
                         simplify = F)
    weights[[1]] <- D_XZ_res$weights; weights[[3]] <- y_X_res$weights
    names(weights) <- c("D_XZ", "D_X", "y_X")
    for (j in 1:nensb) {
      # When the LIE is enforced, compute LIE-conform estimates of E[D|X].
      #     Otherwise use the previously calculated estimates of E[D|X].
      if (enforce_LIE) {
        D_X_res <-  ispred(D_XZ_res$is_fitted[, j], X, Z = NULL,
                           models_DX, ens_type[j],
                           cv_folds, setup_parallel, silent)
        # Residualize
        D_r <- D - D_X_res$is_fitted
        V_r <- D_XZ_res$is_fitted[, j] - D_X_res$is_fitted
      } else {
        # Residualize
        D_r <- D - D_X_res$is_fitted[, j]
        V_r <- D_XZ_res$is_fitted[, j] - D_X_res$is_fitted[, j]
      }#IFELSE

      # Residualize y
      y_r <- y - y_X_res$is_fitted[, j]

      # Compute IV estimate with constructed variables
      iv_fit_j <- tsls(y_r, D_r, V_r)

      # Organize complementary ensemble output
      coef[j] <- iv_fit_j$coef[1]
      iv_fit[[j]] <- iv_fit_j
      if (enforce_LIE) {
        weights[[2]][, j] <- D_X_res$weights
      } else {
        weights[[2]][, j] <- D_X_res$weights[, j]
      }#IFELSE

    }#FOR
    # Name output appropriately by ensemble type
   names(iv_fit) <- ens_type
  }#IF

  # Store complementary ensemble output
  mspe <- list(D_XZ = D_XZ_res$mspe,
               D_X = D_X_res$mspe,
               y_X = y_X_res$mspe)
  anyiv_cv <- D_XZ_res$anyiv_cv
  anyiv <- D_XZ_res$anyiv

  # Organize output
  pl_fit <- list(coef = coef, weights = weights,
                 mspe = mspe, anyiv_cv = anyiv_cv,
                 anyiv = anyiv,
                 models = models,
                 models_DX = models_DX,
                 models_DXZ = models_DXZ,
                 iv_fit = iv_fit,
                 ens_type = ens_type,
                 enforce_LIE = enforce_LIE,
                 nobs = nobs, y = y, D = D)

  # Amend class and return
  class(pl_fit) <- c("pl_iv")
  return(pl_fit)
}#PL_IV
