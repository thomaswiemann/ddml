#' DDML estimator for the Extended Partially Linear IV model.
#'
#' DDML estimator for the Extended Partially Linear IV model.
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
#' @param ensemble_type A string indicating the type of ensemble. Multiple types may
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
#' @return \code{ddml_epliv} returns an object of S3 class
#'     \code{ddml_epliv}.
#'
#' An object of class \code{ddml_epliv}is a list containig the
#'     following components:
#' \describe{
#' \item{\code{coef}}{A vector with the LIE-conform DDML IV coefficent in the
#'     first entry.}
#'     \item{\code{...}}{...}
#' }
#'
#' @examples
#' # Add example here.
#'
#' @export ddml_epliv
ddml_epliv <- function(y, D, Z, X,
                    learners,
                    learners_DXZ = learners,
                    learners_DX = learners,
                    sample_folds = 2,
                    ensemble_type = c("average"),
                    cv_folds = 5,
                    enforce_LIE = TRUE,
                    subsamples = NULL,
                    silent = F) {
  # Data parameters
  nobs <- length(y)

  # Create sample fold tuple
  if (is.null(subsamples)) {
    sampleframe <- rep(1:sample_folds, ceiling(nobs/sample_folds))
    sample_groups <- sample(sampleframe, size=nobs, replace=FALSE)
    subsamples <- sapply(c(1:sample_folds),
                         function(x) {which(sample_groups == x)})
  }#IF
  # In case subsamples are user-provided
  sample_folds <- length(subsamples)

  # Draw samples if not user-supplied
  if (is.null(subsamples)) {
    subsamples <- split(c(1:nobs),
                        sample(rep(c(1:sample_folds),
                                   ceiling(nobs / sample_folds)))[1:nobs])
  }#IF
  sample_folds <- length(subsamples)

  # Compute estimates of E[y|X]
  y_X_res <- crosspred(y, X, Z = NULL,
                       learners, ensemble_type, cv_folds,
                       sample_folds, subsamples,
                       compute_insample_predictions = F,
                       silent)

  # Compute estimates of E[D|X,Z]. Also calculate in-sample predictions when
  #     the LIE is enforced.
  D_XZ_res <- crosspred(D, X, Z,
                        learners_DXZ, ensemble_type, cv_folds,
                        sample_folds, subsamples,
                        compute_insample_predictions = enforce_LIE,
                        silent)

  # When the LIE is not enforced, estimating E[D|X] is straightforward.
  if (!enforce_LIE) {
    D_X_res <- crosspred(D, X, Z = NULL,
                         learners_DX, ensemble_type, cv_folds,
                         sample_folds, subsamples,
                         compute_insample_predictions = F,
                         silent)
  }#IF

  # Check whether multiple ensembles are computed simultaneously.
  multiple_ensembles <- length(ensemble_type) > 1

  # If a single ensemble is calculated, no loops are required.
  if (!multiple_ensembles) {
    # Check whether the law of iterated expectations (LIE) should be enforced.
    #     When the LIE is enforced (recommended), the estimates of E[D|X,Z] are
    #     used for the calculation of the estimates of E[D|X].
    if (enforce_LIE) {
      # Compute LIE-conform estimates of E[D|X]
      D_X_res <- crosspred(D_XZ_res$is_fitted, X, Z = NULL,
                           learners_DX, ensemble_type, cv_folds,
                           sample_folds, subsamples,
                           compute_insample_predictions = F,
                           silent)
    }#IFELSE

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - D_X_res_direct$oos_fitted
    V_r <- D_XZ_res$oos_fitted - D_X_res$oos_fitted

    # Compute IV estimate with constructed variables
    iv_fit <- tsls(y_r, D_r, V_r)

    # Organize complementary ensemble output
    coef <- iv_fit$coef[1]
    weights <- list(y_X = y_X_res$weights,
                    D_X = D_X_res$weights,
                    D_XZ = D_XZ_res$weights)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (multiple_ensembles) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    nensb <- length(ensemble_type)
    coef <- matrix(0, 1, nensb)
    mspe <- iv_fit <- rep(list(1), nensb)
    weights <- list(array(0, dim = c(length(learners), nenseb, sample_folds),
                          simplify = F),
                    array(0, dim = c(length(learners_DX), nenseb, sample_folds),
                          simplify = F),
                    array(0, dim = c(length(learners_DXZ), nenseb, sample_folds),
                          simplify = F))
    weights[[1]] <- y_X_res$weights; weights[[3]] <- D_XZ_res$weights;
    names(weights) <- c("y_X", "D_X", "D_XZ")
    for (j in 1:nensb) {
      # When the LIE is enforced, compute LIE-conform estimates of E[D|X].
      #     Otherwise use the previously calculated estimates of E[D|X].
      if (enforce_LIE) {
        D_X_res <- crosspred(D_XZ_res$is_fitted[[j]], X, Z = NULL,
                             learners_DX, ensemble_type[j], cv_folds,
                             sample_folds, subsamples,
                             compute_insample_predictions = F,
                             silent)
      }#IF

      # Residualize
      if (enforce_LIE) {
        D_r <- D - D_X_res$oos_fitted
        V_r <- D_XZ_res$oos_fitted[, j] - D_X_res$oos_fitted
      } else {
        D_r <- D - D_X_res$oos_fitted[, j]
        V_r <- D_XZ_res$oos_fitted[, j] - D_X_res$oos_fitted[, j]
      }#IFELSE

      # Residualize y
      y_r <- y - y_X_res$oos_fitted[, j]

      # Compute IV estimate with constructed variables
      iv_fit_j <- tsls(y_r, D_r, V_r)

      # Organize complementary ensemble output
      coef[j] <- iv_fit_j$coef[1]
      iv_fit[[j]] <- iv_fit_j
      if (enforce_LIE) {
        weights[[2]][, j, ] <- D_X_res$weights
      } else {
        weights[[2]][, j, ] <- D_X_res$weights[, j, ]
      }#IFELSE

    }#FOR
    # Name output appropriately by ensemble type
    names(iv_fit) <- ensemble_type
  }#IF

  # Store complementary ensemble output
  mspe <- list(y_X = y_X_res$mspe,
               D_X = D_X_res$mspe,
               D_XZ = D_XZ_res$mspe)

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights,
                   mspe = mspe,
                   learners = learners,
                   learners_DXZ = learners_DXZ,
                   learners_DX = learners_DX,
                   iv_fit = iv_fit,
                   subsamples = subsamples,
                   ensemble_type = ensemble_type,
                   enforce_LIE = enforce_LIE,
                   nobs = nobs, y = y, D = D)

  # Amend class and return
  class(ddml_fit) <- c("ddml_epliv")
  return(ddml_fit)
}#ddml_epliv
