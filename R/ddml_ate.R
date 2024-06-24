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
#' @param subsamples_D0,subsamples_D1 List of vectors with sample indices for
#'     cross-fitting, corresponding to untreated and treated observations,
#'     respectively.
#' @param cv_subsamples_list_D0,cv_subsamples_list_D1 List of lists, each
#'     corresponding to a subsample containing vectors with subsample indices
#'     for cross-validation. Arguments are separated for untreated and treated
#'     observations, respectively.
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
#'         \item{\code{learners},\code{learners_DX},
#'             \code{subsamples_D0},\code{subsamples_D1},
#'             \code{cv_subsamples_list_D0},\code{cv_subsamples_list_D1},
#'             \code{ensemble_type}}{Pass-through of
#'             selected user-provided arguments. See above.}
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
                     sample_folds = 2,
                     ensemble_type = "nnls",
                     shortstack = FALSE,
                     cv_folds = 5,
                     custom_ensemble_weights = NULL,
                     custom_ensemble_weights_DX = custom_ensemble_weights,
                     subsamples_D0 = NULL,
                     subsamples_D1 = NULL,
                     cv_subsamples_list_D0 = NULL,
                     cv_subsamples_list_D1 = NULL,
                     silent = FALSE) {
  # Data parameters
  nobs <- length(y)
  is_D0 <- which(D == 0)
  nobs_D0 <- length(is_D0)
  nobs_D1 <- nobs - nobs_D0

  # Create sample fold tuple by treatment
  if (is.null(subsamples_D0) | is.null(subsamples_D1)) {
    subsamples_D0 <- generate_subsamples(nobs_D0, sample_folds)
    subsamples_D1 <- generate_subsamples(nobs_D1, sample_folds)
  }#IF
  sample_folds <- length(subsamples_D0)

  # Create cv-subsamples tuple by treatment
  if (is.null(cv_subsamples_list_D0) | is.null(cv_subsamples_list_D1)) {
    cv_subsamples_list_D0 <- rep(list(NULL), sample_folds)
    cv_subsamples_list_D1 <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      nobs_D0_k <- nobs_D0 - length(subsamples_D0[[k]])
      nobs_D1_k <- nobs_D1 - length(subsamples_D1[[k]])
      cv_subsamples_list_D0[[k]] <- generate_subsamples(nobs_D0_k, cv_folds)
      cv_subsamples_list_D1[[k]] <- generate_subsamples(nobs_D1_k, cv_folds)
    }# FOR
  }#IF

  # Merge subsamples across treatment and create auxilliary control matrix
  subsamples <- subsamples_D0
  cv_subsamples_list <- cv_subsamples_list_D0
  auxilliary_X_D0 <- rep(list(NULL), sample_folds)
  auxilliary_X_D1 <- rep(list(NULL), sample_folds)
  for (k in 1:sample_folds) {
    # Sample folds
    subsamples[[k]] <- sort(c((1:nobs)[is_D0][subsamples_D0[[k]]],
                              (1:nobs)[-is_D0][subsamples_D1[[k]]]))
    # CV folds
    nobs_k <- nobs - length(subsamples[[k]])
    is_D0_k <- which(D[-subsamples[[k]]] == 0)
    is_D1_k <- which(D[-subsamples[[k]]] == 1)
    for (j in 1:cv_folds) {
      indx_D0 <- is_D0_k[cv_subsamples_list_D0[[k]][[j]]]
      indx_D1 <- is_D1_k[cv_subsamples_list_D1[[k]][[j]]]
      cv_subsamples_list[[k]][[j]] <- sort(c(indx_D0, indx_D1))
    }#FOR

    # Auxilliary X
    auxilliary_X_D1[[k]] <- X[-is_D0, , drop=F][subsamples_D1[[k]], , drop=F]
    auxilliary_X_D0[[k]] <- X[is_D0, , drop=F][subsamples_D0[[k]], , drop=F]
  }#FOR

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|D=0,X]
  y_X_D0_res <- get_CEF(y[is_D0], X[is_D0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        cv_subsamples_list = cv_subsamples_list_D0,
                        subsamples = subsamples_D0,
                        silent = silent, progress = "E[Y|D=0,X]: ",
                        auxilliary_X = auxilliary_X_D1)

  # Compute estimates of E[y|D=1,X]
  y_X_D1_res <- get_CEF(y[-is_D0], X[-is_D0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        cv_subsamples_list = cv_subsamples_list_D1,
                        subsamples = subsamples_D1,
                        silent = silent, progress = "E[Y|D=1,X]: ",
                        auxilliary_X = auxilliary_X_D0)

  # Compute estimates of E[D|X]
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights_DX,
                     cv_subsamples_list = cv_subsamples_list,
                     subsamples = subsamples,
                     silent = silent, progress = "E[D|X]: ")

  # Update ensemble type to account for (optional) custom weights
  ensemble_type <- dimnames(y_X_D0_res$weights)[[2]]
  nensb <- ifelse(is.null(ensemble_type), 1, length(ensemble_type))

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # Construct reduced form variables
  g_D0 <- g_D1 <- matrix(0, nobs, nensb)
  g_D0[is_D0, ] <- y_X_D0_res$oos_fitted
  g_D1[-is_D0, ] <- y_X_D1_res$oos_fitted
  if (!multiple_ensembles) {
    for (k in 1:sample_folds) {
      g_D1[is_D0][subsamples_D0[[k]]] <- y_X_D1_res$auxilliary_fitted[[k]]
      g_D0[-is_D0][subsamples_D1[[k]]] <- y_X_D0_res$auxilliary_fitted[[k]]
    }#FOR
  } else {
    for (k in 1:sample_folds) {
      g_D1[is_D0, ][subsamples_D0[[k]], ] <- y_X_D1_res$auxilliary_fitted[[k]]
      g_D0[-is_D0, ][subsamples_D1[[k]], ] <- y_X_D0_res$auxilliary_fitted[[k]]
    }#FOR
  }#IF
  m_X <- D_X_res$oos_fitted

  # Compute the ATE using the constructed variables
  y_copy <- matrix(rep(y, nensb), nobs, nensb)
  D_copy <- matrix(rep(D, nensb), nobs, nensb)
  psi_b <- D_copy * (y_copy - g_D1) / m_X +
    (1 - D_copy) * (y_copy - g_D0) / (1 - m_X) + g_D1 - g_D0
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

  # Organize output
  ddml_fit <- list(ate = ate, weights = weights, mspe = mspe,
                   psi_a = psi_a, psi_b = psi_b,
                   learners = learners,
                   learners_DX = learners_DX,
                   subsamples_D0 = subsamples_D0,
                   subsamples_D1 = subsamples_D1,
                   cv_subsamples_list_D0 = cv_subsamples_list_D0,
                   cv_subsamples_list_D1 = cv_subsamples_list_D1,
                   ensemble_type = ensemble_type)

  # Print estimation progress
  if (!silent) cat("DDML estimation completed. \n")

  # Amend class and return
  class(ddml_fit) <- "ddml_ate"
  return(ddml_fit)
}#DDML_ATE

#' Inference Methods for Treatment Effect Estimators.
#'
#' @description Inference methods for treatment effect estimators.
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
                                                     object$ensemble_type)
  summary_res <- list(coefficients = coefficients,
                      parameter = "ATE")
  class(summary_res) <- "summary.ddml_ate"
  summary_res
}#SUMMARY.DDML_ATE

#' Print Methods for Treatment Effect Estimators.
#'
#' @description Inference methods for treatment effect estimators.
#'
#' @param x An object of class \code{summary.ddml_ate},
#'     \code{summary.ddml_att}, and \code{ddml_late}, as returned by
#'     [ddml::summary.ddml_ate()], [ddml::summary.ddml_att()], and
#'     [ddml::summary.ddml_late()], respectively.
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
print.summary.ddml_ate <- function(x, ...) {
  cat("ATE estimation results: \n \n")
  print(x$coefficients)
}#PRINT.SUMMARY.DDML_ATE
