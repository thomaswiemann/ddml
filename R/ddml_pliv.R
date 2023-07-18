#' Estimator for the Partially Linear IV Model.
#'
#' @family ddml
#'
#' @seealso [AER::ivreg()]
#'
#' @description Estimator for the partially linear IV model.
#'
#' @details \code{ddml_pliv} provides a double/debiased machine learning
#'     estimator for the parameter of interest \eqn{\theta_0} in the partially
#'     linear IV model given by
#'
#' \eqn{Y = \theta_0D + g_0(X) + U,}
#'
#' where \eqn{(Y, D, X, Z, U)} is a random vector such that
#'     \eqn{E[Cov(U, Z\vert X)] = 0} and \eqn{E[Cov(D, Z\vert X)] \neq 0}, and
#'     \eqn{g_0} is an unknown nuisance function.
#'
#' @inheritParams ddml_plm
#' @param Z The instrumental variable.
#' @param learners_ZX Optional argument to allow for different estimators of
#'     \eqn{E[Z|X]}. Setup is identical to \code{learners}.
#'
#' @return \code{ddml_pliv} returns an object of S3 class
#'     \code{ddml_pliv}. An object of class \code{ddml_pliv} is a list
#'     containing the following components:
#'     \describe{
#'         \item{\code{coef}}{A vector with the \eqn{\theta_0} estimates.}
#'         \item{\code{weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{iv_fit}}{Object of class \code{ivreg} from the IV
#'             regression of \eqn{Y - \hat{E}[Y|X]} on
#'             \eqn{D - \hat{E}[D|X]} using \eqn{Z - \hat{E}[Z|X]} as the
#'             instrument.}
#'         \item{\code{learners},\code{learners_DX},\code{learners_ZX},
#'             \code{subsamples},\code{cv_subsamples_list},\code{ensemble_type}
#'             }{Pass-through of selected user-provided arguments. See above.}
#'     }
#' @export
#'
#' @examples
#' # Construct data from the included SIPP_1991 data
#' y = as.matrix(SIPP_1991$net_tfa)
#' D = as.matrix(SIPP_1991$p401)
#' Z = as.matrix(SIPP_1991$e401)
#' X = as.matrix(SIPP_1991[, c("age", "inc", "educ", "fsize",
#'                             "marr", "twoearn", "db", "pira", "hown")])
#' # Estimate the partially linear IV model using a single base learner: Ridge.
#' pliv_fit <- ddml_pliv(y, D, Z, X,
#'                       learners = list(what = mdl_glmnet,
#'                                       args = list(alpha = 0)),
#'                       sample_folds = 2,
#'                       silent = TRUE)
#' pliv_fit$coef
ddml_pliv <- function(y, D, Z, X,
                      learners,
                      learners_DX = learners,
                      learners_ZX = learners,
                      sample_folds = 2,
                      ensemble_type = "average",
                      shortstack = FALSE,
                      cv_folds = 5,
                      subsamples = NULL,
                      cv_subsamples_list = NULL,
                      silent = F) {
  # Data parameters
  nobs <- length(y)
  nlearners <- length(learners)

  # Create sample fold tuple
  if (is.null(subsamples)) {
    subsamples <- generate_subsamples(nobs, sample_folds)
  }#IF
  sample_folds <- length(subsamples)

  # Create cv-subsamples tuple
  if (is.null(cv_subsamples_list) & !shortstack) {
    cv_subsamples_list <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      nobs_k <- nobs - length(subsamples[[k]])
      cv_subsamples_list[[k]] <- generate_subsamples(nobs_k, cv_folds)
    }# FOR
  }#IF

  # Compute estimates of E[y|X]
  y_X_res <- get_CEF(y, X,
                     learners = learners, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     subsamples = subsamples,
                     cv_subsamples_list = cv_subsamples_list,
                     silent = silent, progress = "E[Y|X]: ")

  # Compute estimates of E[Z|X].
  Z_X_res <- get_CEF(y, X,
                     learners = learners_ZX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     subsamples = subsamples,
                     cv_subsamples_list = cv_subsamples_list,
                     silent = silent, progress = "E[Z|X]: ")

  # Compute estimates of E[D|X].
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     subsamples = subsamples,
                     cv_subsamples_list = cv_subsamples_list,
                     silent = silent, progress = "E[D|X]: ")

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- length(ensemble_type) > 1

  # If a single ensemble is calculated, no loops are required.
  if (!multiple_ensembles) {

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - D_X_res$oos_fitted
    V_r <- Z - Z_X_res$oos_fitted

    # Compute IV estimate with constructed variables
    iv_fit <- AER::ivreg(y_r ~ D_r | V_r)

    # Organize complementary ensemble output
    coef <- stats::coef(iv_fit)[2]
    weights <- list(y_X = y_X_res$weights,
                    D_X = D_X_res$weights,
                    Z_X = Z_X_res$weights)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (multiple_ensembles) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    nensb <- length(ensemble_type)
    coef <- matrix(0, 1, nensb)
    iv_fit <- rep(list(1), nensb)
    nlearners <- length(learners)
    nlearners_DX <- length(learners_DX); nlearners_ZX <- length(learners_ZX)
    # Initialize weights
    weights <- list()
    weights[[1]] <- y_X_res$weights; weights[[2]] <- D_X_res$weights
    weights[[3]] <- Z_X_res$weights
    # Assign names for more legible output
    colnames(coef) <- names(iv_fit) <- ensemble_type
    names(weights) <- c("y_X", "D_X", "Z_X")
    # Compute coefficients for each ensemble
    for (j in 1:nensb) {
      # Residualize
      y_r <- y - y_X_res$oos_fitted[, j]
      D_r <- D - D_X_res$oos_fitted[, j]
      V_r <- Z - Z_X_res$oos_fitted[, j]

      # Compute IV estimate with constructed variables
      iv_fit_j <- AER::ivreg(y_r ~ D_r | V_r)

      # Organize complementary ensemble output
      coef[j] <- stats::coef(iv_fit_j)[2]
      iv_fit[[j]] <- iv_fit_j
    }#FOR
    # Name output appropriately by ensemble type
    names(iv_fit) <- ensemble_type
  }#IF

  # Store complementary ensemble output
  mspe <- list(y_X = y_X_res$mspe,
               D_X = D_X_res$mspe,
               Z_X = Z_X_res$mspe)

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights, mspe = mspe,
                   learners = learners,
                   learners_ZX = learners_ZX,
                   learners_DX = learners_DX,
                   iv_fit = iv_fit,
                   subsamples = subsamples,
                   cv_subsamples_list = cv_subsamples_list,
                   ensemble_type = ensemble_type)

  # Amend class and return
  class(ddml_fit) <- c("ddml_pliv")
  return(ddml_fit)
}#DDML_PLIV
