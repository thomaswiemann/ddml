#' Estimator of the Average Treatment Effect on the Treated.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml_att()]
#'
#' @description Estimator of the average treatment effect on the treated.
#'
#' @details \code{ddml_att} provides a double/debiased machine learning
#'     estimator for the average treatment effect on the treated in the
#'     interactive model given by
#'
#' \eqn{Y = g_0(D, X) + U,}
#'
#' where \eqn{(Y, D, X, U)} is a random vector such that
#'     \eqn{\operatorname{supp} D = \{0,1\}}, \eqn{E[U\vert D, X] = 0}, and
#'     \eqn{\Pr(D=1\vert X) \in (0, 1)} with probability 1,
#'     and \eqn{g_0} is an unknown nuisance function.
#'
#' In this model, the average treatment effect on the treated is defined as
#'
#' \eqn{\theta_0^{\textrm{ATT}} \equiv E[g_0(1, X) - g_0(0, X) \vert D = 1]}.
#'
#' @inheritParams ddml_ate
#'
#' @return \code{ddml_att} returns an object of S3 class
#'     \code{ddml_att}. An object of class \code{ddml_att} is a list containing
#'     the following components:
#'     \describe{
#'         \item{\code{att}}{A vector with the average treatment effect
#'             on the treated estimates.}
#'         \item{\code{weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{psi_a}, \code{psi_b}}{Matrices needed for the computation
#'             of scores. Used in [ddml::summary.ddml_att()].}
#'         \item{\code{learners},\code{learners_DX},
#'             \code{subsamples_D0},\code{cv_subsamples_list_D0},
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
#' # Estimate the average treatment effect on the treated using a single base
#' #     learner, ridge.
#' att_fit <- ddml_att(y, D, X,
#'                     learners = list(what = mdl_glmnet,
#'                                     args = list(alpha = 0)),
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(att_fit)
#'
#' # Estimate the average treatment effect on the treated using short-stacking
#' #     with base learners ols, lasso, and ridge.
#' att_fit <- ddml_att(y, D, X,
#'                     learners = list(list(fun = ols),
#'                                     list(fun = mdl_glmnet),
#'                                     list(fun = mdl_glmnet,
#'                                          args = list(alpha = 0))),
#'                     ensemble_type = 'nnls',
#'                     shortstack = TRUE,
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(att_fit)
ddml_att <- function(y, D, X,
                     learners,
                     learners_DX = learners,
                     sample_folds = 2,
                     ensemble_type = "nnls",
                     shortstack = FALSE,
                     cv_folds = 5,
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
  nensb <- length(ensemble_type)

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

    # Auxilliary X (only need treated observations)
    auxilliary_X_D1[[k]] <- X[-is_D0, , drop=F][subsamples_D1[[k]], , drop=F]
  }#FOR

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|D=0,X]
  y_X_D0_res <- get_CEF(y[is_D0], X[is_D0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        cv_subsamples_list = cv_subsamples_list_D0,
                        subsamples = subsamples_D0,
                        silent = silent, progress = "E[Y|D=0,X]: ",
                        auxilliary_X = auxilliary_X_D1)

  # Compute estimates of E[D|X]
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     cv_subsamples_list = cv_subsamples_list,
                     subsamples = subsamples,
                     silent = silent, progress = "E[D|X]: ")

  # Compute estimates of E[D] -- simple computation of averages here
  D_res <- get_CEF(D, matrix(1, nobs, 1),
                   learners = list(what = ols),
                   ensemble_type = "average",
                   shortstack = FALSE,
                   cv_subsamples_list = NULL,
                   subsamples = subsamples,
                   silent = TRUE)

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # Construct reduced form variables
  g_D0 <- matrix(0, nobs, nensb)
  g_D0[is_D0, ] <- y_X_D0_res$oos_fitted
  if (!multiple_ensembles) {
    for (k in 1:sample_folds) {
      g_D0[-is_D0][subsamples_D1[[k]]] <- y_X_D0_res$auxilliary_fitted[[k]]
    }#FOR
  } else {
    for (k in 1:sample_folds) {
      g_D0[-is_D0, ][subsamples_D1[[k]], ] <- y_X_D0_res$auxilliary_fitted[[k]]
    }#FOR
  }#IF
  m_X <- D_X_res$oos_fitted

  # Compute the ATT using the constructed variables
  y_copy <- matrix(rep(y, nensb), nobs, nensb)
  D_copy <- matrix(rep(D, nensb), nobs, nensb)
  p_copy <- matrix(rep(D_res$oos_fitted, nensb), nobs, nensb)
  psi_b <- D_copy * (y_copy - g_D0) / p_copy -
    m_X * (1 - D_copy) * (y_copy - g_D0) / (p_copy * (1 - m_X))
  psi_a <- -D_copy / p_copy
  att <- -colMeans(psi_b) / colMeans(psi_a)
  names(att) <- ensemble_type

  # Organize complementary ensemble output
  weights <- list(y_X_D0 = y_X_D0_res$weights,
                  D_X = D_X_res$weights)

  # Store complementary ensemble output
  mspe <- list(y_X_D0 = y_X_D0_res$mspe,
               D_X = D_X_res$mspe)

  # Organize output
  ddml_fit <- list(att = att, weights = weights, mspe = mspe,
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
  class(ddml_fit) <- "ddml_att"
  return(ddml_fit)
}#DDML_ATT

#' @rdname summary.ddml_ate
#'
#' @export
summary.ddml_att <- function(object, ...) {
  # Check whether stacking was used, replace ensemble type if TRUE
  single_learner <- ("what" %in% names(object$learners))
  if (single_learner) object$ensemble_type <- " "
  # Compute and print inference results
  cat("ATT estimation results: \n \n")
  organize_interactive_inf_results(coef = object$att,
                                   psi_a = object$psi_a,
                                   psi_b = object$psi_b,
                                   ensemble_type = object$ensemble_type)
}#SUMMARY.DDML_ATT
