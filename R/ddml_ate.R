#' Estimators of Average Treatment Effects.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml_ate()], [ddml::summary.ddml_att()]
#'
#' @description Estimators of the average treatment effect and the average
#'     treatment effect on the treated.
#'
#' @details \code{ddml_ate} and \code{ddml_att} provide double/debiased machine
#'     learning  estimators for the average treatment effect and the average
#'     treatment effect on the treated, respectively, in the interactive model
#'     given by
#'
#' \eqn{Y = g_0(D, X) + U,}
#'
#' where \eqn{(Y, D, X, U)} is a random vector such that
#'     \eqn{\operatorname{supp} D = \{0,1\}}, \eqn{E[U\vert D, X] = 0}, and
#'     \eqn{\Pr(D=1\vert X) \in (0, 1)} with probability 1,
#'     and \eqn{g_0} is an unknown nuisance function.
#'
#' In this model, the average treatment effect is defined as
#'
#' \eqn{\theta_0^{\textrm{ATE}} \equiv E[g_0(1, X) - g_0(0, X)]}.
#'
#' and the average treatment effect on the treated is defined as
#'
#' \eqn{\theta_0^{\textrm{ATT}} \equiv E[g_0(1, X) - g_0(0, X)\vert D = 1]}.
#'
#' @inheritParams ddml_plm
#' @param D The binary endogenous variable of interest.
#' @param subsamples_byD List of two lists corresponding to the two treatment
#'     levels. Each list contains vectors with sample indices for
#'     cross-fitting.
#' @param cv_subsamples_byD List of two lists, each corresponding to one of the
#'     two treatment levels. Each of the two lists contains lists, each
#'     corresponding to a subsample and contains vectors with subsample indices
#'     for cross-validation.
#' @param trim Number in (0, 1) for trimming the estimated propensity scores at
#'     \code{trim} and \code{1-trim}.
#'
#' @return \code{ddml_ate} and \code{ddml_att} return an object of S3 class
#'     \code{ddml_ate} and \code{ddml_att}, respectively. An object of class
#'     \code{ddml_ate} or \code{ddml_att} is a list containing
#'     the following components:
#'     \describe{
#'         \item{\code{ate} / \code{att}}{A vector with the average treatment
#'             effect / average treatment effect on the treated estimates.}
#'         \item{\code{weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{psi_a}, \code{psi_b}}{Matrices needed for the computation
#'             of scores. Used in [ddml::summary.ddml_ate()] or
#'             [ddml::summary.ddml_att()].}
#'         \item{\code{oos_pred}}{List of matrices, providing the reduced form
#'             predicted values.}
#'         \item{\code{learners},\code{learners_DX},\code{cluster_variable},
#'             \code{subsamples_D0},\code{subsamples_D1},
#'             \code{cv_subsamples_list_D0},\code{cv_subsamples_list_D1},
#'             \code{ensemble_type}}{Pass-through of
#'             selected user-provided arguments. See above.}
#'     }
#' @export
#'
#' @references
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
#' # Estimate the average treatment effect using a single base learner, ridge.
#' ate_fit <- ddml_ate(y, D, X,
#'                     learners = list(what = mdl_glmnet,
#'                                     args = list(alpha = 0)),
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(ate_fit)
#'
#' # Estimate the average treatment effect using short-stacking with base
#' #     learners ols, lasso, and ridge. We can also use custom_ensemble_weights
#' #     to estimate the ATE using every individual base learner.
#' weights_everylearner <- diag(1, 3)
#' colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
#' ate_fit <- ddml_ate(y, D, X,
#'                     learners = list(list(fun = ols),
#'                                     list(fun = mdl_glmnet),
#'                                     list(fun = mdl_glmnet,
#'                                          args = list(alpha = 0))),
#'                     ensemble_type = 'nnls',
#'                     custom_ensemble_weights = weights_everylearner,
#'                     shortstack = TRUE,
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(ate_fit)
ddml_ate <- function(y, D, X,
                     learners,
                     learners_DX = learners,
                     sample_folds = 10,
                     ensemble_type = "nnls",
                     shortstack = FALSE,
                     cv_folds = 10,
                     custom_ensemble_weights = NULL,
                     custom_ensemble_weights_DX = custom_ensemble_weights,
                     cluster_variable = seq_along(y),
                     subsamples_byD = NULL,
                     cv_subsamples_byD = NULL,
                     trim = 0.01,
                     silent = FALSE) {
  # Data parameters
  nobs <- length(y)
  is_D0 <- which(D == 0)

  # Create sample and cv-fold tuples
  cf_indxs <- get_crossfit_indices(cluster_variable = cluster_variable, D = D,
                                   sample_folds = sample_folds,
                                   cv_folds = cv_folds,
                                   subsamples_byD = subsamples_byD,
                                   cv_subsamples_byD = cv_subsamples_byD)

  # Create tuple for extrapolated fitted values
  aux_indxs <- get_auxiliary_indx(cf_indxs$subsamples_byD, D)

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|D=0,X]
  y_X_D0_res <- get_CEF(y[is_D0], X[is_D0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        subsamples = cf_indxs$subsamples_byD[[1]],
                        cv_subsamples_list = cf_indxs$cv_subsamples_byD[[1]],
                        silent = silent, progress = "E[Y|D=0,X]: ",
                        auxiliary_X = get_auxiliary_X(aux_indxs[[1]], X))

  # Compute estimates of E[y|D=1,X]
  y_X_D1_res <- get_CEF(y[-is_D0], X[-is_D0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        subsamples = cf_indxs$subsamples_byD[[2]],
                        cv_subsamples_list = cf_indxs$cv_subsamples_byD[[2]],
                        silent = silent, progress = "E[Y|D=1,X]: ",
                        auxiliary_X = get_auxiliary_X(aux_indxs[[2]], X))

  # Compute estimates of E[D|X]
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights_DX,
                     subsamples = cf_indxs$subsamples,
                     cv_subsamples_list = cf_indxs$cv_subsamples_list,
                     silent = silent, progress = "E[D|X]: ")

  # Update ensemble type to account for (optional) custom weights
  ensemble_type <- dimnames(y_X_D0_res$weights)[[2]]
  nensb <- ifelse(is.null(ensemble_type), 1, length(ensemble_type))

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # Construct reduced form variables
  g_X_byD <- extrapolate_CEF(D = D,
                             CEF_res_byD = list(list(y_X_D0_res, d=0),
                                                list(y_X_D1_res, d=1)),
                             aux_indxs = aux_indxs)
  m_X <- D_X_res$oos_fitted

  # Trim propensity scores, return warnings
  m_X_tr <- trim_propensity_scores(m_X, trim, ensemble_type)

  # Compute the ATE using the constructed variables
  y_copy <- matrix(rep(y, nensb), nobs, nensb)
  D_copy <- matrix(rep(D, nensb), nobs, nensb)
  psi_b <- D_copy * (y_copy - g_X_byD[, , 2]) / m_X_tr -
    (1 - D_copy) * (y_copy -  g_X_byD[, , 1]) / (1 - m_X_tr) +
    g_X_byD[, , 2] - g_X_byD[, , 1]
  ate <- colMeans(psi_b)
  names(ate) <- ensemble_type

  # Also set psi_a scores for easier computation of summary.ddml_ate
  psi_a <- matrix(-1, nobs, nensb)

  # Organize complementary ensemble output
  weights <- list(y_X_D0 = y_X_D0_res$weights,
                  y_X_D1 = y_X_D1_res$weights,
                  D_X = D_X_res$weights)

  # Store complementary ensemble output
  mspe <- list(y_X_D0 = y_X_D0_res$mspe,
               y_X_D1 = y_X_D1_res$mspe,
               D_X = D_X_res$mspe)

  # Organize reduced form predicted values
  oos_pred <- list(EY_D0_X = g_X_byD[, , 1],
                   EY_D1_X = g_X_byD[, , 2],
                   ED_X = m_X)

  # Organize output
  ddml_fit <- list(ate = ate, weights = weights, mspe = mspe,
                   psi_a = psi_a, psi_b = psi_b,
                   oos_pred = oos_pred,
                   learners = learners,
                   learners_DX = learners_DX,
                   cluster_variable = cluster_variable,
                   subsamples_byD = subsamples_byD,
                   cv_subsamples_byD = cv_subsamples_byD,
                   ensemble_type = ensemble_type)

  # Print estimation progress
  if (!silent) cat("DDML estimation completed. \n")

  # Amend class and return
  class(ddml_fit) <- "ddml_ate"
  return(ddml_fit)
}#DDML_ATE

#' Inference Methods for Treatment Effect Estimators.
#'
#' @description Inference methods for treatment effect estimators. By default,
#'     standard errors are heteroskedasiticty-robust. If the \code{ddml}
#'     estimator was computed using a \code{cluster_variable}, the standard
#'     errors are also cluster-robust by default.
#'
#' @param object An object of class \code{ddml_ate}, \code{ddml_att}, and
#'     \code{ddml_late}, as fitted by [ddml::ddml_ate()], [ddml::ddml_att()],
#'     and [ddml::ddml_late()], respectively.
#' @param ... Currently unused.
#'
#' @return A matrix with inference results.
#'
#' @export
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the average treatment effect using a single base learner, ridge.
#' ate_fit <- ddml_ate(y, D, X,
#'                     learners = list(what = mdl_glmnet,
#'                                     args = list(alpha = 0)),
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(ate_fit)
summary.ddml_ate <- function(object, ...) {
  # Check whether stacking was used, replace ensemble type if TRUE
  single_learner <- ("what" %in% names(object$learners))
  if (single_learner) object$ensemble_type <- " "
  # Compute and return inference results
  coefficients <- organize_interactive_inf_results(coef = object$ate,
                                                   psi_a = object$psi_a,
                                                   psi_b = object$psi_b,
                                                   ensemble_type =
                                                     object$ensemble_type,
                                                   cluster_variable =
                                                     object$cluster_variable)
  class(coefficients) <- c("summary.ddml_ate", class(coefficients))
  coefficients
}#SUMMARY.DDML_ATE

#' Print Methods for Treatment Effect Estimators.
#'
#' @description Print methods for treatment effect estimators.
#'
#' @param x An object of class \code{summary.ddml_ate},
#'     \code{summary.ddml_att}, and \code{ddml_late}, as returned by
#'     [ddml::summary.ddml_ate()], [ddml::summary.ddml_att()], and
#'     [ddml::summary.ddml_late()], respectively.
#' @param digits The number of significant digits used for printing.
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
#' # Estimate the average treatment effect using a single base learner, ridge.
#' ate_fit <- ddml_ate(y, D, X,
#'                     learners = list(what = mdl_glmnet,
#'                                     args = list(alpha = 0)),
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(ate_fit)
print.summary.ddml_ate <- function(x, digits = 3, ...) {
  cat("ATE estimation results: \n \n")
  class(x) <- class(x)[-1]
  print(x, digits = digits)
}#PRINT.SUMMARY.DDML_ATE
