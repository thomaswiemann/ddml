#' Estimator of the Local Average Treatment Effect.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml_late()]
#'
#' @description Estimator of the local average treatment effect.
#'
#' @details \code{ddml_late} provides a double/debiased machine learning
#'     estimator for the local average treatment effect in the interactive model
#'     given by
#'
#' \eqn{Y = g_0(D, X) + U,}
#'
#' where \eqn{(Y, D, X, Z, U)} is a random vector such that
#'     \eqn{\operatorname{supp} D = \operatorname{supp} Z = \{0,1\}},
#'     \eqn{E[U\vert X, Z] = 0}, \eqn{E[Var(E[D\vert X, Z]\vert X)] \neq 0},
#'     \eqn{\Pr(Z=1\vert X) \in (0, 1)} with probability 1,
#'     \eqn{p_0(1, X) \geq p_0(0, X)} with probability 1 where
#'     \eqn{p_0(Z, X) \equiv \Pr(D=1\vert Z, X)}, and
#'     \eqn{g_0} is an unknown nuisance function.
#'
#' In this model, the local average treatment effect is defined as
#'
#' \eqn{\theta_0^{\textrm{LATE}} \equiv
#'     E[g_0(1, X) - g_0(0, X)\vert p_0(1, X) > p(0, X)]}.
#'
#' @inheritParams ddml_ate
#' @param Z Binary instrumental variable.
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
#' @param learners_DXZ,learners_ZX Optional arguments to allow for different
#'     estimators of \eqn{E[D \vert X, Z]}, \eqn{E[Z \vert X]}. Setup is
#'     identical to \code{learners}.
#' @param custom_ensemble_weights_DXZ,custom_ensemble_weights_ZX Optional
#'     arguments to allow for different
#'     custom ensemble weights for \code{learners_DXZ},\code{learners_ZX}. Setup
#'     is identical to \code{custom_ensemble_weights}. Note:
#'     \code{custom_ensemble_weights} and
#'     \code{custom_ensemble_weights_DXZ},\code{custom_ensemble_weights_ZX} must
#'     have the same number of columns.
#' @param subsamples_byZ List of two lists corresponding to the two instrument
#'     levels. Each list contains vectors with sample indices for
#'     cross-fitting.
#' @param cv_subsamples_byZ List of two lists, each corresponding to one of the
#'     two instrument levels. Each of the two lists contains lists, each
#'     corresponding to a subsample and contains vectors with subsample indices
#'     for cross-validation.
#'
#' @return \code{ddml_late} returns an object of S3 class
#'     \code{ddml_late}. An object of class \code{ddml_late} is a list
#'     containing the following components:
#'     \describe{
#'         \item{\code{late}}{A vector with the average treatment effect
#'             estimates.}
#'         \item{\code{weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{psi_a}, \code{psi_b}}{Matrices needed for the computation
#'             of scores. Used in [ddml::summary.ddml_late()].}
#'         \item{\code{oos_pred}}{List of matrices, providing the reduced form
#'             predicted values.}
#'         \item{\code{learners},\code{learners_DXZ},\code{learners_ZX},
#'             \code{cluster_variable},\code{subsamples_Z0},
#'             \code{subsamples_Z1},\code{cv_subsamples_list_Z0},
#'             \code{cv_subsamples_list_Z1},\code{ensemble_type}}{Pass-through
#'             of selected user-provided arguments. See above.}
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
#' Imbens G, Angrist J (1004). "Identification and Estimation of Local Average
#'     Treatment Effects." Econometrica, 62(2), 467-475.
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
#' # Estimate the local average treatment effect using a single base learner,
#' #     ridge.
#' late_fit <- ddml_late(y, D, Z, X,
#'                       learners = list(what = mdl_glmnet,
#'                                       args = list(alpha = 0)),
#'                       sample_folds = 2,
#'                       silent = TRUE)
#' summary(late_fit)
#'
#' # Estimate the local average treatment effect using short-stacking with base
#' #     learners ols, lasso, and ridge. We can also use custom_ensemble_weights
#' #     to estimate the ATE using every individual base learner.
#' weights_everylearner <- diag(1, 3)
#' colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
#' late_fit <- ddml_late(y, D, Z, X,
#'                       learners = list(list(fun = ols),
#'                                       list(fun = mdl_glmnet),
#'                                       list(fun = mdl_glmnet,
#'                                            args = list(alpha = 0))),
#'                       ensemble_type = 'nnls',
#'                       custom_ensemble_weights = weights_everylearner,
#'                       shortstack = TRUE,
#'                       sample_folds = 2,
#'                       silent = TRUE)
#' summary(late_fit)
ddml_late <- function(y, D, Z, X,
                      learners,
                      learners_DXZ = learners,
                      learners_ZX = learners,
                      sample_folds = 10,
                      ensemble_type = "nnls",
                      shortstack = FALSE,
                      cv_folds = 10,
                      custom_ensemble_weights = NULL,
                      custom_ensemble_weights_DXZ = custom_ensemble_weights,
                      custom_ensemble_weights_ZX = custom_ensemble_weights,
                      cluster_variable = seq_along(y),
                      subsamples_byZ = NULL,
                      cv_subsamples_byZ = NULL,
                      trim = 0.01,
                      silent = FALSE) {
  # Data parameters
  nobs <- length(y)
  is_Z0 <- which(Z == 0)

  # Create sample and cv-fold tuples
  cf_indxs <- get_crossfit_indices(cluster_variable = cluster_variable, D = Z,
                                   sample_folds = sample_folds,
                                   cv_folds = cv_folds,
                                   subsamples_byD = subsamples_byZ,
                                   cv_subsamples_byD = cv_subsamples_byZ)

  # Create tuple for extrapolated fitted values
  aux_indxs <- get_auxiliary_indx(cf_indxs$subsamples_byD, Z)

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|Z=0,X]
  y_X_Z0_res <- get_CEF(y[is_Z0], X[is_Z0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        subsamples = cf_indxs$subsamples_byD[[1]],
                        cv_subsamples_list = cf_indxs$cv_subsamples_byD[[1]],
                        silent = silent, progress = "E[Y|Z=0,X]: ",
                        auxiliary_X = get_auxiliary_X(aux_indxs[[1]], X))

  # Compute estimates of E[y|Z=1,X]
  y_X_Z1_res <- get_CEF(y[-is_Z0], X[-is_Z0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        subsamples = cf_indxs$subsamples_byD[[2]],
                        cv_subsamples_list = cf_indxs$cv_subsamples_byD[[2]],
                        silent = silent, progress = "E[Y|Z=1,X]: ",
                        auxiliary_X = get_auxiliary_X(aux_indxs[[2]], X))

  # Check for perfect non-compliance
  if (all(D[Z==0] == 0)) {
    # Artificially construct values for subsample with Z=0
    D_X_Z0_res <- list(NULL)
    D_X_Z0_res$oos_fitted <- rep(0, length(is_Z0))
    D_X_Z0_res$auxiliary_fitted <-
      lapply(y_X_Z0_res$auxiliary_fitted, function (x) {x * 0})
    if (!silent) cat("E[D|Z=0,X]: perfect non-compliance -- Done! \n")
  } else {
    # Compute estimates of E[D|Z=0,X]
    D_X_Z0_res <- get_CEF(D[is_Z0], X[is_Z0, , drop = F],
                          learners = learners_DXZ,
                          ensemble_type = ensemble_type,
                          shortstack = shortstack,
                          custom_ensemble_weights = custom_ensemble_weights_DXZ,
                          subsamples = cf_indxs$subsamples_byD[[1]],
                          cv_subsamples_list = cf_indxs$cv_subsamples_byD[[1]],
                          silent = silent, progress = "E[Y|Z=0,X]: ",
                          auxiliary_X = get_auxiliary_X(aux_indxs[[1]], X))
  }#IFELSE

  # Check for perfect compliance
  if (all(D[Z==1] == 1)) {
    # Artificially construct values for subsample with Z=0
    D_X_Z1_res <- list(NULL)
    D_X_Z1_res$oos_fitted <- rep(0, nobs - length(is_Z0))
    D_X_Z1_res$auxiliary_fitted <-
      lapply(y_X_Z1_res$auxiliary_fitted, function (x) {x * 0})
    if (!silent) cat("E[D|Z=1,X]: perfect compliance -- Done! \n")
  } else {
    # Compute estimates of E[D|Z=1,X]
    D_X_Z1_res <- get_CEF(D[-is_Z0], X[-is_Z0, , drop = F],
                          learners = learners_DXZ,
                          ensemble_type = ensemble_type,
                          shortstack = shortstack,
                          custom_ensemble_weights = custom_ensemble_weights_DXZ,
                          subsamples = cf_indxs$subsamples_byD[[2]],
                          cv_subsamples_list = cf_indxs$cv_subsamples_byD[[2]],
                          silent = silent, progress = "E[Y|Z=0,X]: ",
                          auxiliary_X = get_auxiliary_X(aux_indxs[[2]], X))
  }#IFELSE

  # Compute estimates of E[Z|X]
  Z_X_res <- get_CEF(Z, X,
                     learners = learners_ZX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights_ZX,
                     subsamples = cf_indxs$subsamples,
                     cv_subsamples_list = cf_indxs$cv_subsamples_list,
                     compute_insample_predictions = F,
                     silent = silent, progress = "E[Z|X]: ")

  # Update ensemble type to account for (optional) custom weights
  ensemble_type <- dimnames(y_X_Z0_res$weights)[[2]]
  nensb <- ifelse(is.null(ensemble_type), 1, length(ensemble_type))

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # Construct reduced form variables
  l_X_byZ <- extrapolate_CEF(D = Z,
                             CEF_res_byD = list(list(y_X_Z0_res, d=0),
                                                list(y_X_Z1_res, d=1)),
                             aux_indxs = aux_indxs)
  p_X_byZ <- extrapolate_CEF(D = Z,
                             CEF_res_byD = list(list(D_X_Z0_res, d=0),
                                                list(D_X_Z1_res, d=1)),
                             aux_indxs = aux_indxs)
  r_X <- Z_X_res$oos_fitted

  # Trim propensity scores, return warnings
  r_X_tr <- trim_propensity_scores(r_X, trim, ensemble_type)

  # Compute the ATE using the constructed variables
  y_copy <- matrix(rep(y, nensb), nobs, nensb)
  D_copy <- matrix(rep(D, nensb), nobs, nensb)
  Z_copy <- matrix(rep(Z, nensb), nobs, nensb)
  psi_b <- Z_copy * (y_copy - l_X_byZ[, , 2]) / r_X_tr -
    (1 - Z_copy) * (y_copy - l_X_byZ[, , 1]) / (1 - r_X_tr) +
    l_X_byZ[, , 2] - l_X_byZ[, , 1]
  psi_a <- -(Z_copy * (D_copy - p_X_byZ[, , 2]) / r_X_tr -
    (1 - Z_copy) * (D_copy - p_X_byZ[, , 1]) / (1 - r_X_tr) +
      p_X_byZ[, , 2] - p_X_byZ[, , 1])
  numerator <- colMeans(psi_b)
  denominator <- colMeans(psi_a)
  late <- -numerator / denominator
  names(late) <- ensemble_type

  # Organize complementary ensemble output
  weights <- list(y_X_Z0 = y_X_Z0_res$weights,
                  y_X_Z1 = y_X_Z1_res$weights,
                  D_X_Z0 = D_X_Z0_res$weights,
                  D_X_Z1 = D_X_Z1_res$weights,
                  Z_X = Z_X_res$weights)

  # Store complementary ensemble output
  mspe <- list(y_X_Z0 = y_X_Z0_res$mspe,
               y_X_Z1 = y_X_Z1_res$mspe,
               D_X_Z0 = D_X_Z0_res$mspe,
               D_X_Z1 = D_X_Z1_res$mspe,
               Z_X = Z_X_res$mspe)

  # Organize reduced form predicted values
  oos_pred <- list(EY_Z0_X = l_X_byZ[, , 1], EY_Z1_X = l_X_byZ[, , 2],
                   ED_Z0_X = p_X_byZ[, , 1], ED_Z1_X = p_X_byZ[, , 2],
                   EZ_X = r_X)

  # Organize output
  ddml_fit <- list(late = late, weights = weights, mspe = mspe,
                   psi_a = psi_a, psi_b = psi_b,
                   oos_pred = oos_pred,
                   learners = learners,
                   learners_DXZ = learners_DXZ,
                   learners_ZX = learners_ZX,
                   cluster_variable = cluster_variable,
                   subsamples_byZ = subsamples_byZ,
                   cv_subsamples_byZ = cv_subsamples_byZ,
                   ensemble_type = ensemble_type)

  # Print estimation progress
  if (!silent) cat("DDML estimation completed. \n")

  # Amend class and return
  class(ddml_fit) <- "ddml_late"
  return(ddml_fit)
}#DDML_LATE

#' @rdname summary.ddml_ate
#'
#' @export
summary.ddml_late <- function(object, ...) {
  # Check whether stacking was used, replace ensemble type if TRUE
  single_learner <- ("what" %in% names(object$learners))
  if (single_learner) object$ensemble_type <- " "
  # Compute and print inference results
  coefficients <- organize_interactive_inf_results(coef = object$late,
                                                   psi_a = object$psi_a,
                                                   psi_b = object$psi_b,
                                                   ensemble_type =
                                                     object$ensemble_type,
                                                   cluster_variable =
                                                     object$cluster_variable)
  class(coefficients) <- c("summary.ddml_late", class(coefficients))
  coefficients
}#SUMMARY.DDML_LATE

#' @rdname print.summary.ddml_ate
#'
#' @export
print.summary.ddml_late <- function(x, digits = 3, ...) {
  cat("LATE estimation results: \n \n")
  class(x) <- class(x)[-1]
  print(x, digits = digits)
}#PRINT.SUMMARY.DDML_LATE

