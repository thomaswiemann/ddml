#' Estimator for the Flexible Partially Linear IV Model.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml_fpliv()], [AER::ivreg()]
#'
#' @description Estimator for the flexible partially linear IV model.
#'
#' @details \code{ddml_fpliv} provides a double/debiased machine learning
#'     estimator for the parameter of interest \eqn{\theta_0} in the partially
#'     linear IV model given by
#'
#' \eqn{Y = \theta_0D + g_0(X) + U,}
#'
#' where \eqn{(Y, D, X, Z, U)} is a random vector such that
#'     \eqn{E[U\vert X, Z] = 0} and \eqn{E[Var(E[D\vert X, Z]\vert X)] \neq 0},
#'     and \eqn{g_0} is an unknown nuisance function.
#'
#' @inheritParams ddml_pliv
#' @param Z A (sparse) matrix of instruments.
#' @param learners_DXZ,learners_DX Optional arguments to allow for different
#'     estimators of \eqn{E[D \vert X, Z]}, \eqn{E[D \vert X]}. Setup is
#'     identical to \code{learners}.
#' @param custom_ensemble_weights_DXZ,custom_ensemble_weights_DX Optional
#'     arguments to allow for different
#'     custom ensemble weights for \code{learners_DXZ},\code{learners_DX}. Setup
#'     is identical to \code{custom_ensemble_weights}. Note:
#'     \code{custom_ensemble_weights} and
#'     \code{custom_ensemble_weights_DXZ},\code{custom_ensemble_weights_DX} must
#'     have the same number of columns.
#' @param enforce_LIE Indicator equal to 1 if the law of iterated expectations
#'     is enforced in the first stage.
#'
#' @return \code{ddml_fpliv} returns an object of S3 class
#'     \code{ddml_fpliv}. An object of class \code{ddml_fpliv} is a list
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
#'             regression of \eqn{Y - \hat{E}[Y\vert X]} on
#'             \eqn{D - \hat{E}[D\vert X]} using
#'             \eqn{\hat{E}[D\vert X,Z] - \hat{E}[D\vert X]} as the instrument.}
#'         \item{\code{learners},\code{learners_DX},\code{learners_DXZ},
#'             \code{cluster_variable},\code{subsamples},
#'             \code{cv_subsamples_list},\code{ensemble_type}}{Pass-through of
#'             selected user-provided arguments. See above.}
#'     }
#' @export
#'
#' @references
#' Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). "Model Averaging and 
#'     Double Machine Learning." Journal of Applied Econometrics, 40(3): 249-269.
#'
#' Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
#'     Robins J (2018). "Double/debiased machine learning for treatment and
#'     structural parameters." The Econometrics Journal, 21(1), C1-C68.
#'
#' Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2), 241-259.
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' Z = AE98[, "samesex", drop = FALSE]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the partially linear IV model using a single base learner: Ridge.
#' fpliv_fit <- ddml_fpliv(y, D, Z, X,
#'                         learners = list(what = mdl_glmnet,
#'                                         args = list(alpha = 0)),
#'                         sample_folds = 2,
#'                         silent = TRUE)
#' summary(fpliv_fit)
ddml_fpliv <- function(y, D, Z, X,
                       learners,
                       learners_DXZ = learners,
                       learners_DX = learners,
                       sample_folds = 10,
                       ensemble_type = "nnls",
                       shortstack = FALSE,
                       cv_folds = 10,
                       enforce_LIE = TRUE,
                       custom_ensemble_weights = NULL,
                       custom_ensemble_weights_DXZ = custom_ensemble_weights,
                       custom_ensemble_weights_DX = custom_ensemble_weights,
                       cluster_variable = seq_along(y),
                       subsamples = NULL,
                       cv_subsamples_list = NULL,
                       silent = FALSE) {
  # Data parameters
  nobs <- length(y)
  nensb_raw <- length(ensemble_type) # number of ensembles w/o custom weights

  # Check for multivariate endogenous variables
  D <- as.matrix(D)
  nD <- ncol(D)

  # Create sample and cv-fold tuples
  cf_indxs <- get_crossfit_indices(cluster_variable = cluster_variable,
                                   sample_folds = sample_folds,
                                   cv_folds = cv_folds,
                                   subsamples = subsamples,
                                   cv_subsamples_list = cv_subsamples_list)

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|X]
  y_X_res <- get_CEF(y, X, Z = NULL,
                     learners = learners, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights,
                     subsamples = cf_indxs$subsamples,
                     cv_subsamples_list = cf_indxs$cv_subsamples_list,
                     compute_insample_predictions = F,
                     silent = silent, progress = "E[Y|X]: ")

  # Compute estimates of E[D|X,Z]. Also calculate in-sample predictions when
  #     the LIE is enforced.
  D_XZ_res_list <- list()
  for (k in 1:nD) {
    D_XZ_res_list[[k]] <- get_CEF(D[, k, drop = F], X, Z,
                                  learners = learners_DXZ,
                                  ensemble_type = ensemble_type,
                                  shortstack = shortstack,
                                  custom_ensemble_weights =
                                    custom_ensemble_weights_DXZ,
                                  subsamples = cf_indxs$subsamples,
                                  cv_subsamples_list =
                                    cf_indxs$cv_subsamples_list,
                                  compute_insample_predictions = enforce_LIE,
                                  silent = silent,
                                  progress = paste0("E[D", k, "|X,Z]: "))
  }#FOR

  # When the LIE is not enforced, estimating E[D|X] is straightforward.
  if (!enforce_LIE) {
    D_X_res_list <- list()
    for (k in 1:nD) {
      D_X_res_list[[k]] <- get_CEF(D[, k, drop = F], X, Z = NULL,
                                   learners = learners_DX,
                                   ensemble_type = ensemble_type,
                                   shortstack = shortstack,
                                   custom_ensemble_weights =
                                     custom_ensemble_weights_DX,
                                   subsamples = cf_indxs$subsamples,
                                   cv_subsamples_list =
                                     cf_indxs$cv_subsamples_list,
                                   compute_insample_predictions = F,
                                   silent = silent,
                                   progress = paste0("E[D", k, "|X]: "))
    }#FOR
  }#IF

  # Update ensemble type to account for (optional) custom weights
  ensemble_type <- dimnames(y_X_res$weights)[[2]]
  nensb <- length(ensemble_type)

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # If a single ensemble is calculated, no loops are required.
  if (!multiple_ensembles) {
    # Check whether the law of iterated expectations (LIE) should be enforced.
    #     When the LIE is enforced (recommended), the estimates of E[D|X,Z] are
    #     used for the calculation of the estimates of E[D|X].
    if (enforce_LIE) {
      D_X_res_list <- list()
      for (k in 1:nD) {
        D_X_res_list[[k]] <-
          get_CEF(D_XZ_res_list[[k]]$is_fitted, X, Z = NULL,
                  learners = learners_DX,
                  ensemble_type = ensemble_type,
                  shortstack = shortstack,
                  subsamples = cf_indxs$subsamples,
                  cv_subsamples_list = cf_indxs$cv_subsamples_list,
                  compute_insample_predictions = F,
                  silent = silent,
                  progress = paste0("E[D", k, "|X]: "),
                  shortstack_y = D_XZ_res_list[[k]]$oos_fitted)
      }#FOR
    }#IFELSE

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - get_oosfitted(D_X_res_list)
    V_r <- get_oosfitted(D_XZ_res_list) - get_oosfitted(D_X_res_list)

    # Compute IV estimate with constructed variables
    iv_fit <- AER::ivreg(y_r ~ D_r | V_r)

    # Organize complementary ensemble output
    coef <- stats::coef(iv_fit)[-1]
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (multiple_ensembles) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    coef <- matrix(0, nD, nensb)
    iv_fit <- rep(list(1), nensb)
    nlearners <- length(learners)
    nlearners_DX <- length(learners_DX); nlearners_DXZ <- length(learners_DXZ)
    # Assign names for more legible output
    colnames(coef) <- names(iv_fit) <- ensemble_type

    # Create intermediate weight matrices for enforce_LIE = TRUE
    if (enforce_LIE) {
      # weights
      weights_DX <- array(0, dim = c(nlearners_DX, nensb, sample_folds))
      dimnames(weights_DX) <- dimnames(y_X_res$weights)
      weights_DX <- replicate(nD, weights_DX, simplify = FALSE)
    }#IF

    # Compute coefficients for each ensemble
    for (j in 1:nensb) {
      # When the LIE is enforced, compute LIE-conform estimates of E[D|X].
      #     Otherwise use the previously calculated estimates of E[D|X].
      if (enforce_LIE) {
        D_X_res_list <- list()
        for (k in 1:nD) {
          progress_jk <- paste0("E[D", k, "|X] (", ensemble_type[j], "): ")

          # Check whether j is a custom ensemble specification. Necessary to
          #     assign correct corresponding custom_weights vector.
          if (j <= nensb_raw) { # j is not a custom specification
            D_X_res_list[[k]] <-
              get_CEF(D_XZ_res_list[[k]]$is_fitted[[j]], X, Z = NULL,
                      learners = learners_DX,
                      ensemble_type = ensemble_type[j],
                      shortstack = shortstack,
                      subsamples = cf_indxs$subsamples,
                      cv_subsamples_list = cf_indxs$cv_subsamples_list,
                      compute_insample_predictions = F,
                      silent = silent,
                      progress = progress_jk,
                      shortstack_y = D_XZ_res_list[[k]]$oos_fitted[, j])
          } else { # j is a custom specification
            D_X_res_list[[k]] <-
              get_CEF(D_XZ_res_list[[k]]$is_fitted[[j]], X, Z = NULL,
                      learners = learners_DX,
                      ensemble_type = "average",
                      shortstack = shortstack,
                      custom_ensemble_weights =
                        custom_ensemble_weights_DX[, j - nensb_raw, drop = F],
                      subsamples = cf_indxs$subsamples,
                      cv_subsamples_list = cf_indxs$cv_subsamples_list,
                      compute_insample_predictions = F,
                      silent = silent,
                      progress = progress_jk,
                      shortstack_y = D_XZ_res_list[[k]]$oos_fitted[, j])
            # Remove "average" oos_fitted and weights
            D_X_res_list[[k]]$oos_fitted <- D_X_res_list[[k]]$oos_fitted[, -1]
            D_X_res_list[[k]]$weights <- D_X_res_list[[k]]$weights[, -1, ,
                                                                   drop = F]
          }#IFELSE


        }#FOR
      }#IF

      # Residualize
      if (enforce_LIE) {
        D_r <- D - get_oosfitted(D_X_res_list)
        V_r <- get_oosfitted(D_XZ_res_list, j) - get_oosfitted(D_X_res_list)
      } else {
        D_r <- D - get_oosfitted(D_X_res_list, j)
        V_r <- get_oosfitted(D_XZ_res_list, j) - get_oosfitted(D_X_res_list, j)
      }#IFELSE

      # Residualize y
      y_r <- y - y_X_res$oos_fitted[, j]

      # Compute IV estimate with constructed variables
      iv_fit_j <- AER::ivreg(y_r ~ D_r | V_r)

      # Organize complementary ensemble output
      coef[, j] <- stats::coef(iv_fit_j)[-1]
      iv_fit[[j]] <- iv_fit_j
      if (enforce_LIE) {
        for (k in 1:nD) weights_DX[[k]][, j, ] <- D_X_res_list[[k]]$weights
      }#IF
    }#FOR
  }#IF

  # Store complementary ensemble output
  weights <- list(y_X = y_X_res$weights)
  mspe <- list(y_X = y_X_res$mspe)
  for (k in 1:nD){
    if (enforce_LIE & multiple_ensembles) {
      weights[[paste0("D", k, "_X")]] <- weights_DX[[k]]
    } else {
      weights[[paste0("D", k, "_X")]] <- D_X_res_list[[k]]$weights
    }#IFELSE
    #mspe[[paste0("D", k, "_X")]] <- D_X_res_list[[k]]$mspe
  }#FOR
  for (k in 1:nD){
    weights[[paste0("D", k, "_XZ")]] <- D_XZ_res_list[[k]]$weights
    mspe[[paste0("D", k, "_XZ")]] <- D_XZ_res_list[[k]]$mspe
  }#FOR

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights, mspe = mspe,
                   learners = learners,
                   learners_DXZ = learners_DXZ,
                   learners_DX = learners_DX,
                   iv_fit = iv_fit,
                   cluster_variable = cluster_variable,
                   subsamples = subsamples,
                   cv_subsamples_list = cv_subsamples_list,
                   ensemble_type = ensemble_type,
                   enforce_LIE = enforce_LIE)

  # Print estimation progress
  if (!silent) cat("DDML estimation completed. \n")

  # Amend class and return
  class(ddml_fit) <- "ddml_fpliv"
  return(ddml_fit)
}#DDML_FPLIV

#' @rdname summary.ddml_plm
#'
#' @export
summary.ddml_fpliv <- function(object, ...) {
  # Check whether stacking was used, replace ensemble type if TRUE
  single_learner <- ("what" %in% names(object$learners))
  if (single_learner) object$ensemble_type <- "single base learner"
  # Compute and print inference results
  coefficients <- organize_inf_results(fit_obj_list = object$iv_fit,
                                       ensemble_type = object$ensemble_type,
                                       cluster_variable =
                                         object$cluster_variable,
                                       ...)
  class(coefficients) <- c("summary.ddml_fpliv", class(coefficients))
  coefficients
}#SUMMARY.DDML_FPLIV

#' @rdname print.summary.ddml_plm
#'
#' @export
print.summary.ddml_fpliv <- function(x, digits = 3, ...) {
  cat("FPLIV estimation results: \n \n")
  class(x) <- class(x)[-1]
  print(x, digits = digits)
}#PRINT.SUMMARY.DDML_FPLIV
