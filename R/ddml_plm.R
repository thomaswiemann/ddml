#' Estimator for the Partially Linear Model.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml_plm()]
#'
#' @description Estimator for the partially linear model.
#'
#' @details \code{ddml_plm} provides a double/debiased machine learning
#'     estimator for the parameter of interest \eqn{\theta_0} in the partially
#'     linear model given by
#'
#' \eqn{Y = \theta_0D + g_0(X) + U,}
#'
#' where \eqn{(Y, D, X, U)} is a random vector such that
#'     \eqn{E[Cov(U, D\vert X)] = 0} and \eqn{E[Var(D\vert X)] \neq 0}, and
#'     \eqn{g_0} is an unknown nuisance function.
#'
#' @param y The outcome variable.
#' @param D A matrix of endogenous variables.
#' @param X A (sparse) matrix of control variables.
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
#'     }
#'     Omission of the \code{args} element results in default arguments being
#'     used in \code{fun}. Omission of \code{assign_X} results in inclusion of
#'     all variables in \code{X}.
#' @param learners_DX Optional argument to allow for different estimators of
#'     \eqn{E[D|X]}. Setup is identical to \code{learners}.
#' @param sample_folds Number of cross-fitting folds.
#' @param ensemble_type Ensemble method to combine base learners into final
#'     estimate of the conditional expectation functions. Possible values are:
#'     \itemize{
#'         \item{\code{"nnls"} Non-negative least squares.}
#'         \item{\code{"nnls1"} Non-negative least squares with the constraint
#'             that all weights sum to one.}
#'         \item{\code{"singlebest"} Select base learner with minimum MSPE.}
#'         \item{\code{"ols"} Ordinary least squares.}
#'         \item{\code{"average"} Simple average over base learners.}
#'     }
#'     Multiple ensemble types may be passed as a vector of strings.
#' @param shortstack Boolean to use short-stacking.
#' @param cv_folds Number of folds used for cross-validation in ensemble
#'     construction.
#' @param custom_ensemble_weights A numerical matrix with user-specified
#'     ensemble weights. Each column corresponds to a custom ensemble
#'     specification, each row corresponds to a base learner in \code{learners}
#'     (in chronological order). Optional column names are used to name the
#'     estimation results corresponding the custom ensemble specification.
#' @param custom_ensemble_weights_DX Optional argument to allow for different
#'     custom ensemble weights for \code{learners_DX}. Setup is identical to
#'     \code{custom_ensemble_weights}. Note: \code{custom_ensemble_weights} and
#'     \code{custom_ensemble_weights_DX} must have the same number of columns.
#' @param cluster_variable A vector of cluster indices.
#' @param subsamples List of vectors with sample indices for cross-fitting.
#' @param cv_subsamples_list List of lists, each corresponding to a subsample
#'     containing vectors with subsample indices for cross-validation.
#' @param silent Boolean to silence estimation updates.
#'
#' @return \code{ddml_plm} returns an object of S3 class
#'     \code{ddml_plm}. An object of class \code{ddml_plm} is a list containing
#'     the following components:
#'     \describe{
#'         \item{\code{coef}}{A vector with the \eqn{\theta_0} estimates.}
#'         \item{\code{weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{ols_fit}}{Object of class \code{lm} from the second
#'             stage regression of \eqn{Y - \hat{E}[Y|X]} on
#'             \eqn{D - \hat{E}[D|X]}.}
#'         \item{\code{learners},\code{learners_DX},\code{cluster_variable},
#'             \code{subsamples}, \code{cv_subsamples_list},
#'             \code{ensemble_type}}{Pass-through of selected user-provided
#'             arguments. See above.}
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
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the partially linear model using a single base learner, ridge.
#' plm_fit <- ddml_plm(y, D, X,
#'                     learners = list(what = mdl_glmnet,
#'                                     args = list(alpha = 0)),
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(plm_fit)
#'
#' # Estimate the partially linear model using short-stacking with base learners
#' #     ols, lasso, and ridge. We can also use custom_ensemble_weights
#' #     to estimate the ATE using every individual base learner.
#' weights_everylearner <- diag(1, 3)
#' colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
#' plm_fit <- ddml_plm(y, D, X,
#'                     learners = list(list(fun = ols),
#'                                     list(fun = mdl_glmnet),
#'                                     list(fun = mdl_glmnet,
#'                                          args = list(alpha = 0))),
#'                     ensemble_type = 'nnls',
#'                     custom_ensemble_weights = weights_everylearner,
#'                     shortstack = TRUE,
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(plm_fit)
ddml_plm <- function(y, D, X,
                     learners,
                     learners_DX = learners,
                     sample_folds = 10,
                     ensemble_type = "nnls",
                     shortstack = FALSE,
                     cv_folds = 10,
                     custom_ensemble_weights = NULL,
                     custom_ensemble_weights_DX = custom_ensemble_weights,
                     cluster_variable = seq_along(y),
                     subsamples = NULL,
                     cv_subsamples_list = NULL,
                     silent = FALSE) {
  # Data parameters
  nobs <- length(y)
  ncustom <- ncol(custom_ensemble_weights)
  ncustom <- ifelse(is.null(ncustom), 0, ncustom)
  nensb <- length(ensemble_type) + ncustom

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
  y_X_res <- get_CEF(y, X,
                     learners = learners,
                     ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights,
                     subsamples = cf_indxs$subsamples,
                     cv_subsamples_list = cf_indxs$cv_subsamples_list,
                     silent = silent, progress = "E[Y|X]: ")

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

    # Compute OLS estimate with constructed variables
    ols_fit <- stats::lm(y_r ~ D_r)

    # Organize complementary ensemble output
    coef <- stats::coef(ols_fit)[-1]
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (multiple_ensembles) {
    # Iterate over ensemble type. Compute DDML estimate for each.
    coef <- matrix(0, nD, nensb)
    mspe <- ols_fit <- rep(list(1), nensb)
    nlearners <- length(learners); nlearners_DX <- length(learners_DX)

    # Compute coefficients for each ensemble
    for (j in 1:nensb) {
      # Residualize
      D_r <- D - get_oosfitted(D_X_res_list, j)

      # Residualize y
      y_r <- y - y_X_res$oos_fitted[, j]

      # Compute OLS estimate with constructed variables
      ols_fit_j <- stats::lm(y_r ~ D_r)

      # Organize complementary ensemble output
      coef[, j] <- stats::coef(ols_fit_j)[-1]
      ols_fit[[j]] <- ols_fit_j
    }#FOR

    # Assign names for more legible output
    colnames(coef) <- names(ols_fit) <- dimnames(y_X_res$weights)[[2]]
    rownames(coef) <- names(ols_fit_j$coefficients)[-1]
  }#IF

  # Store complementary ensemble output
  weights <- list(y_X = y_X_res$weights)
  mspe <- list(y_X = y_X_res$mspe)
  for (k in 1:nD){
    weights[[paste0("D", k, "_X")]] <- D_X_res_list[[k]]$weights
    mspe[[paste0("D", k, "_X")]] <- D_X_res_list[[k]]$mspe
  }#FOR

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights, mspe = mspe,
                   learners = learners,
                   learners_DX = learners_DX,
                   ols_fit = ols_fit,
                   cluster_variable = cluster_variable,
                   subsamples = cf_indxs$subsamples,
                   cv_subsamples_list = cf_indxs$cv_subsamples_list,
                   ensemble_type = ensemble_type)

  # Print estimation progress
  if (!silent) cat("DDML estimation completed. \n")

  # Amend class and return
  class(ddml_fit) <- "ddml_plm"
  return(ddml_fit)
}#DDML_PLM

#' Inference Methods for Partially Linear Estimators.
#'
#' @seealso [sandwich::vcovHC()], [sandwich::vcovCL()]
#'
#' @description Inference methods for partially linear estimators. Simple
#'     wrapper for [sandwich::vcovHC()] and [sandwich::vcovCL()]. Default
#'     standard errors are heteroskedasiticty-robust. If the \code{ddml}
#'     estimator was computed using a \code{cluster_variable}, the standard
#'     errors are also cluster-robust by default.
#'
#' @param object An object of class \code{ddml_plm}, \code{ddml_pliv}, or
#'     \code{ddml_fpliv} as fitted by [ddml::ddml_plm()], [ddml::ddml_pliv()],
#'     and [ddml::ddml_fpliv()], respectively.
#' @param ... Additional arguments passed to \code{vcovHC} and \code{vcovCL}.
#'     See [sandwich::vcovHC()] and [sandwich::vcovCL()] for a complete list of
#'     arguments.
#'
#' @return An array with inference results for each \code{ensemble_type}.
#'
#' @references
#' Zeileis A (2004). "Econometric Computing with HC and HAC Covariance Matrix
#'     Estimators.” Journal of Statistical Software, 11(10), 1-17.
#'
#' Zeileis A (2006). “Object-Oriented Computation of Sandwich Estimators.”
#'     Journal of Statistical Software, 16(9), 1-16.
#'
#' Zeileis A, Köll S, Graham N (2020). “Various Versatile Variances: An
#'     Object-Oriented Implementation of Clustered Covariances in R.” Journal of
#'     Statistical Software, 95(1), 1-36.
#'
#' @export
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the partially linear model using a single base learner, ridge.
#' plm_fit <- ddml_plm(y, D, X,
#'                     learners = list(what = mdl_glmnet,
#'                                     args = list(alpha = 0)),
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(plm_fit)
summary.ddml_plm <- function(object, ...) {
  # Check whether stacking was used, replace ensemble type if TRUE
  single_learner <- ("what" %in% names(object$learners))
  if (single_learner) object$ensemble_type <- "single base learner"
  # Compute and print inference results
  coefficients <- organize_inf_results(fit_obj_list = object$ols_fit,
                                       ensemble_type = object$ensemble_type,
                                       cluster_variable =
                                         object$cluster_variable,
                                       ...)
  class(coefficients) <- c("summary.ddml_plm", class(coefficients))
  coefficients
}#SUMMARY.DDML_PLM

#' Print Methods for Treatment Effect Estimators.
#'
#' @description Print methods for treatment effect estimators.
#'
#' @param x An object of class \code{summary.ddml_plm},
#'     \code{summary.ddml_pliv}, and \code{summary.ddml_fpliv}, as
#'     returned by [ddml::summary.ddml_plm()], [ddml::summary.ddml_pliv()],
#'     and [ddml::summary.ddml_fpliv()], respectively.
#' @param digits Number of significant digits used for priniting.
#' @param ... Currently unused.
#'
#' @return NULL.
#'
#' @export
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the partially linear model using a single base learner, ridge.
#' plm_fit <- ddml_plm(y, D, X,
#'                     learners = list(what = mdl_glmnet,
#'                                     args = list(alpha = 0)),
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(plm_fit)
print.summary.ddml_plm <- function(x, digits = 3, ...) {
  cat("PLM estimation results: \n \n")
  class(x) <- class(x)[-1]
  print(x, digits = digits)
}#PRINT.SUMMARY.DDML_PLM
