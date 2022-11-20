#' Compute PLM estimators.
#'
#' Compute PLM estimators.
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
#' @param silent A boolean indicating whether current models and folds should be
#'     printed to the console.
#'
#' @return \code{plm} returns an object of S3 class
#'     \code{plm}.
#'
#' An object of class \code{plm}is a list containig the
#'     following components:
#' \describe{
#' \item{\code{coef}}{A vector with the DDML PLM coefficent in the
#'     first entry.}
#'     \item{\code{...}}{...}
#' }
#'
#' @export plm
plm <- function(y, D, X,
                     models,
                     ens_type = c("average"),
                     cv_folds = 5,
                     setup_parallel = list(type = 'dynamic', cores = 1),
                     silent = F) {
  # Data parameters
  nobs <- length(y)
  nmodels <- length(models); nensb <- length(ens_type)

  # Compute estimates of E[y|X]
  y_X_res <- ispred(y, X, Z = NULL,
                    models, ens_type,
                    cv_folds, setup_parallel, silent)

  # Compute estimates of E[D|X].
  D_X_res <- ispred(D, X, Z = NULL,
                    models, ens_type,
                    cv_folds, setup_parallel, silent)

  # Check whether multiple ensembles are computed simultaneously
  sim_ens <- length(ens_type) > 1

  # If a single ensemble is calculated, no loops are required.
  if (!sim_ens) {

    # Residualize
    y_r <- y - y_X_res$is_fitted
    D_r <- D - D_X_res$is_fitted

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
    weights <- replicate(2, array(0, dim = c(nmodels, nensb)),
                         simplify = F)
    weights[[2]] <- D_X_res$weights; weights[[1]] <- y_X_res$weights
    names(weights) <- c("y_X", "D_X")
    for (j in 1:nensb) {
      # Residualize
      D_r <- D - D_X_res$is_fitted[, j]

      # Residualize y
      y_r <- y - y_X_res$is_fitted[, j]

      # Compute IV estimate with constructed variables
      ols_fit_j <- ols(y_r, D_r)

      # Organize complementary ensemble output
      coef[j] <- ols_fit_j$coef[1]
      ols_fit[[j]] <- ols_fit_j
      weights[[2]][, j] <- D_X_res$weights[, j]
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
                   ens_type = ens_type,
                   nobs = nobs, y = y, D = D)

  # Amend class and return
  class(ddml_fit) <- c("plm")
  return(ddml_fit)
}#PLM
