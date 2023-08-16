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
#' @param learners_ZX Optional argument to allow for different estimators of
#'     \eqn{E[Z\vert X]}. Setup is identical to \code{learners}.
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
#'             \code{subsamples},\code{cv_subsamples_list},\code{ensemble_type}
#'             }{Pass-through of selected user-provided arguments. See above.}
#'     }
#' @export
#'
#' @references
#' Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2023). "ddml: Double/debiased
#'     machine learning in Stata." \url{https://arxiv.org/abs/2301.09397}
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
#'
#' # Estimate the partially linear IV model using short-stacking with base
#' #     ols, lasso, and ridge.
#' pliv_fit <- ddml_pliv(y, D, Z, X,
#'                       learners = list(list(fun = ols),
#'                                       list(fun = mdl_glmnet),
#'                                       list(fun = mdl_glmnet,
#'                                            args = list(alpha = 0))),
#'                       ensemble_type = 'nnls',
#'                       shortstack = TRUE,
#'                       sample_folds = 2,
#'                       silent = TRUE)
#' summary(pliv_fit)
ddml_pliv <- function(y, D, Z, X,
                      learners,
                      learners_DX = learners,
                      learners_ZX = learners,
                      sample_folds = 2,
                      ensemble_type = "nnls",
                      shortstack = FALSE,
                      cv_folds = 5,
                      subsamples = NULL,
                      cv_subsamples_list = NULL,
                      silent = F) {
  # Data parameters
  nobs <- length(y)
  nlearners <- length(learners)

  # Check for multivariate endogenous variables
  D <- as.matrix(D)
  nD <- ncol(D)

  # Check for multivariate instruments
  Z <- as.matrix(Z)
  nZ <- ncol(Z)

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

  # Compute estimates of E[Z|X], loop through instruments
  Z_X_res_list <- list()
  for (k in 1:nZ) {
    Z_X_res_list[[k]] <- get_CEF(Z[, k, drop = F], X,
                                 learners = learners_ZX,
                                 ensemble_type = ensemble_type,
                                 shortstack = shortstack,
                                 subsamples = subsamples,
                                 cv_subsamples_list = cv_subsamples_list,
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
                                 subsamples = subsamples,
                                 cv_subsamples_list = cv_subsamples_list,
                                 silent = silent,
                                 progress = paste0("E[D", k, "|X]: "))
  }#FOR

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- length(ensemble_type) > 1

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
    nensb <- length(ensemble_type)
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
  cat("PLIV estimation results: \n \n")
  organize_inf_results(fit_obj_list = object$iv_fit,
                       ensemble_type = object$ensemble_type,
                       ...)
}#SUMMARY.DDML_PLIV
