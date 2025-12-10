#' Estimator for the Partially Linear IV Model.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml_pliv()], [AER::ivreg()]
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
#' @param Z A matrix of instruments.
#' @param learners May take one of two forms, depending on whether a single
#'     learner or stacking with multiple learners is used for estimation of the
#'     conditional expectation functions.
#'     If a single learner is used, \code{learners} is a list with two named
#'     elements:
#'     \itemize{
#'         \item{\code{what} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}.}
#'     }
#'     If stacking with multiple learners is used, \code{learners} is a list of
#'     lists, each containing four named elements:
#'     \itemize{
#'         \item{\code{fun} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}.}
#'         \item{\code{assign_X} An optional vector of column indices
#'             corresponding to control variables in \code{X} that are passed to
#'             the base learner.}
#'         \item{\code{assign_Z} An optional vector of column indices
#'             corresponding to instruments in \code{Z} that are passed to the
#'             base learner.}
#'     }
#'     Omission of the \code{args} element results in default arguments being
#'     used in \code{fun}. Omission of \code{assign_X} (and/or \code{assign_Z})
#'     results in inclusion of all variables in \code{X} (and/or \code{Z}).
#' @param learners_DX,learners_ZX Optional arguments to allow for different
#'     base learners for estimation of \eqn{E[D|X]}, \eqn{E[Z|X]}. Setup is
#'     identical to \code{learners}.
#' @param custom_ensemble_weights_DX,custom_ensemble_weights_ZX Optional
#'     arguments to allow for different
#'     custom ensemble weights for \code{learners_DX},\code{learners_ZX}. Setup
#'     is identical to \code{custom_ensemble_weights}. Note:
#'     \code{custom_ensemble_weights} and
#'     \code{custom_ensemble_weights_DX},\code{custom_ensemble_weights_ZX} must
#'     have the same number of columns.
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
#'             regression of \eqn{Y - \hat{E}[Y\vert X]} on
#'             \eqn{D - \hat{E}[D\vert X]} using \eqn{Z - \hat{E}[Z\vert X]} as
#'             the instrument. See also [AER::ivreg()] for details.}
#'         \item{\code{learners},\code{learners_DX},\code{learners_ZX},
#'             \code{cluster_variable}, \code{subsamples},
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
#' Kleiber C, Zeileis A (2008). Applied Econometrics with R. Springer-Verlag,
#'     New York.
#'
#' Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2), 241-259.
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' Z = AE98[, "samesex"]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the partially linear IV model using a single base learner, ridge.
#' pliv_fit <- ddml_pliv(y, D, Z, X,
#'                       learners = list(what = mdl_glmnet,
#'                                       args = list(alpha = 0)),
#'                       sample_folds = 2,
#'                       silent = TRUE)
#' summary(pliv_fit)
ddml_pliv <- function(y, D, Z, X,
                      learners,
                      learners_DX = learners,
                      learners_ZX = learners,
                      sample_folds = 10,
                      ensemble_type = "nnls",
                      shortstack = FALSE,
                      cv_folds = 10,
                      custom_ensemble_weights = NULL,
                      custom_ensemble_weights_DX = custom_ensemble_weights,
                      custom_ensemble_weights_ZX = custom_ensemble_weights,
                      cluster_variable = seq_along(y),
                      subsamples = NULL,
                      cv_subsamples_list = NULL,
                      silent = FALSE) {
  # Data parameters
  nobs <- length(y)
  nlearners <- length(learners)
  ncustom <- ncol(custom_ensemble_weights)
  ncustom <- ifelse(is.null(ncustom), 0, ncustom)

  # Check for multivariate endogenous variables
  D <- as.matrix(D)
  nD <- ncol(D)

  # Check for multivariate instruments
  Z <- as.matrix(Z)
  nZ <- ncol(Z)

  # Create sample and cv-fold tuples
  cf_indxs <- get_crossfit_indices(cluster_variable = cluster_variable,
                                   sample_folds = sample_folds,
                                   cv_folds = cv_folds,
                                   subsamples = subsamples,
                                   cv_subsamples_list = cv_subsamples_list)

  # Compute estimates of E[y|X]
  y_X_res <- get_CEF(y, X,
                     learners = learners, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights,
                     subsamples = cf_indxs$subsamples,
                     cv_subsamples_list = cf_indxs$cv_subsamples_list,
                     silent = silent, progress = "E[Y|X]: ")

  # Compute estimates of E[Z|X], loop through instruments
  Z_X_res_list <- list()
  for (k in 1:nZ) {
    Z_X_res_list[[k]] <- get_CEF(Z[, k, drop = F], X,
                                 learners = learners_ZX,
                                 ensemble_type = ensemble_type,
                                 shortstack = shortstack,
                                 custom_ensemble_weights =
                                   custom_ensemble_weights_ZX,
                                 subsamples = cf_indxs$subsamples,
                                 cv_subsamples_list =
                                   cf_indxs$cv_subsamples_list,
                                 silent = silent,
                                 progress = paste0("E[Z", k, "|X]: "))
  }#FOR

  # Compute estimates of E[D|X], loop through endogenous variables
  D_X_res_list <- list()
  for (k in 1:nD) {
    D_X_res_list[[k]] <- get_CEF(D[, k, drop = F], X,
                                 learners = learners_DX,
                                 ensemble_type = ensemble_type,
                                 shortstack = shortstack,
                                 custom_ensemble_weights =
                                   custom_ensemble_weights_DX,
                                 subsamples = cf_indxs$subsamples,
                                 cv_subsamples_list =
                                   cf_indxs$cv_subsamples_list,
                                 silent = silent,
                                 progress = paste0("E[D", k, "|X]: "))
  }#FOR

  # Update ensemble type to account for (optional) custom weights
  ensemble_type <- dimnames(y_X_res$weights)[[2]]
  nensb <- length(ensemble_type)

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # If a single ensemble is calculated, no loops are required.
  if (!multiple_ensembles) {

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - get_oosfitted(D_X_res_list)
    V_r <- Z - get_oosfitted(Z_X_res_list)

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
    nlearners_DX <- length(learners_DX); nlearners_ZX <- length(learners_ZX)

    # Compute coefficients for each ensemble
    for (j in 1:nensb) {
      # Residualize
      y_r <- y - y_X_res$oos_fitted[, j]
      D_r <- D - get_oosfitted(D_X_res_list, j)
      V_r <- Z - get_oosfitted(Z_X_res_list, j)

      # Compute IV estimate with constructed variables
      iv_fit_j <- AER::ivreg(y_r ~ D_r | V_r)

      # Organize complementary ensemble output
      coef[, j] <- stats::coef(iv_fit_j)[-1]
      iv_fit[[j]] <- iv_fit_j
    }#FOR
    # Assign names for more legible output
    colnames(coef) <- names(iv_fit) <- ensemble_type
    rownames(coef) <- names(iv_fit_j$coefficients)[-1]
  }#IF

  # Store complementary ensemble output
  weights <- list(y_X = y_X_res$weights)
  mspe <- list(y_X = y_X_res$mspe)
  for (k in 1:nD){
    weights[[paste0("D", k, "_X")]] <- D_X_res_list[[k]]$weights
    mspe[[paste0("D", k, "_X")]] <- D_X_res_list[[k]]$mspe
  }#FOR
  for (k in 1:nZ){
    weights[[paste0("Z", k, "_X")]] <- Z_X_res_list[[k]]$weights
    mspe[[paste0("Z", k, "_X")]] <- Z_X_res_list[[k]]$mspe
  }#FOR

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights, mspe = mspe,
                   learners = learners,
                   learners_ZX = learners_ZX,
                   learners_DX = learners_DX,
                   iv_fit = iv_fit,
                   cluster_variable = cluster_variable,
                   subsamples = subsamples,
                   cv_subsamples_list = cv_subsamples_list,
                   ensemble_type = ensemble_type)

  # Amend class and return
  class(ddml_fit) <- "ddml_pliv"
  return(ddml_fit)
}#DDML_PLIV

#' @rdname summary.ddml_plm
#'
#' @export
summary.ddml_pliv <- function(object, ...) {
  # Check whether stacking was used, replace ensemble type if TRUE
  single_learner <- ("what" %in% names(object$learners))
  if (single_learner) object$ensemble_type <- "single base learner"
  # Compute and print inference results
  coefficients <- organize_inf_results(fit_obj_list = object$iv_fit,
                                       ensemble_type = object$ensemble_type,
                                       cluster_variable =
                                         object$cluster_variable,
                                       ...)
  class(coefficients) <- c("summary.ddml_pliv", class(coefficients))
  coefficients
}#SUMMARY.DDML_PLIV

#' @rdname print.summary.ddml_plm
#'
#' @export
print.summary.ddml_pliv <- function(x, digits = 3, ...) {
  cat("PLIV estimation results: \n \n")
  class(x) <- class(x)[-1]
  print(x, digits = digits)
}#PRINT.SUMMARY.DDML_PLIV
