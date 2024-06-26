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
#' @param subsamples_Z0,subsamples_Z1 List of vectors with sample indices for
#'     cross-fitting, corresponding to observations with \eqn{Z=0} and
#'     \eqn{Z=1}, respectively.
#' @param cv_subsamples_list_Z0,cv_subsamples_list_Z1 List of lists, each
#'     corresponding to a subsample containing vectors with subsample indices
#'     for cross-validation. Arguments are separated for observations with
#'     \eqn{Z=0} and \eqn{Z=1}, respectively.
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
#'             \code{subsamples_Z0},\code{subsamples_Z1},
#'             \code{cv_subsamples_list_Z0},\code{cv_subsamples_list_Z1},
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
                      sample_folds = 2,
                      ensemble_type = "nnls",
                      shortstack = FALSE,
                      cv_folds = 5,
                      custom_ensemble_weights = NULL,
                      custom_ensemble_weights_DXZ = custom_ensemble_weights,
                      custom_ensemble_weights_ZX = custom_ensemble_weights,
                      subsamples_Z0 = NULL,
                      subsamples_Z1 = NULL,
                      cv_subsamples_list_Z0 = NULL,
                      cv_subsamples_list_Z1 = NULL,
                      trim = 0.01,
                      silent = FALSE) {
  # Data parameters
  nobs <- length(y)
  is_Z0 <- which(Z == 0)
  nobs_Z0 <- length(is_Z0)
  nobs_Z1 <- nobs - nobs_Z0

  # Create sample fold tuple by treatment
  if (is.null(subsamples_Z0) | is.null(subsamples_Z1)) {
    subsamples_Z0 <- generate_subsamples(nobs_Z0, sample_folds)
    subsamples_Z1 <- generate_subsamples(nobs_Z1, sample_folds)
  }#IF
  sample_folds <- length(subsamples_Z0)

  # Create cv-subsamples tuple by treatment
  if (is.null(cv_subsamples_list_Z0) | is.null(cv_subsamples_list_Z1)) {
    cv_subsamples_list_Z0 <- rep(list(NULL), sample_folds)
    cv_subsamples_list_Z1 <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      nobs_Z0_k <- nobs_Z0 - length(subsamples_Z0[[k]])
      nobs_Z1_k <- nobs_Z1 - length(subsamples_Z1[[k]])
      cv_subsamples_list_Z0[[k]] <- generate_subsamples(nobs_Z0_k, cv_folds)
      cv_subsamples_list_Z1[[k]] <- generate_subsamples(nobs_Z1_k, cv_folds)
    }# FOR
  }#IF

  # Merge subsamples across treatment and create auxiliary control matrix
  subsamples <- subsamples_Z0
  cv_subsamples_list <- cv_subsamples_list_Z0
  auxilliary_X_Z0 <- rep(list(NULL), sample_folds)
  auxilliary_X_Z1 <- rep(list(NULL), sample_folds)
  for (k in 1:sample_folds) {
    # Sample folds
    subsamples[[k]] <- sort(c((1:nobs)[is_Z0][subsamples_Z0[[k]]],
                              (1:nobs)[-is_Z0][subsamples_Z1[[k]]]))
    # CV folds
    nobs_k <- nobs - length(subsamples[[k]])
    is_Z0_k <- which(Z[-subsamples[[k]]] == 0)
    is_Z1_k <- which(Z[-subsamples[[k]]] == 1)
    for (j in 1:cv_folds) {
      indx_Z0 <- is_Z0_k[cv_subsamples_list_Z0[[k]][[j]]]
      indx_Z1 <- is_Z1_k[cv_subsamples_list_Z1[[k]][[j]]]
      cv_subsamples_list[[k]][[j]] <- sort(c(indx_Z0, indx_Z1))
    }#FOR

    # Auxilliary X
    auxilliary_X_Z1[[k]] <- X[-is_Z0, , drop=F][subsamples_Z1[[k]], , drop=F]
    auxilliary_X_Z0[[k]] <- X[is_Z0, , drop=F][subsamples_Z0[[k]], , drop=F]
  }#FOR

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|Z=0,X]
  y_X_Z0_res <- get_CEF(y[is_Z0], X[is_Z0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        cv_subsamples_list = cv_subsamples_list_Z0,
                        subsamples = subsamples_Z0,
                        silent = silent, progress = "E[Y|Z=0,X]: ",
                        auxilliary_X = auxilliary_X_Z1)

  # Compute estimates of E[y|Z=1,X]
  y_X_Z1_res <- get_CEF(y[-is_Z0], X[-is_Z0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        cv_subsamples_list = cv_subsamples_list_Z1,
                        subsamples = subsamples_Z1,
                        silent = silent, progress = "E[Y|Z=1,X]: ",
                        auxilliary_X = auxilliary_X_Z0)

  # Check for perfect non-compliance
  if (all(D[Z==0] == 0)) {
    # Artificially construct values for subsample with Z=0
    D_X_Z0_res <- list(NULL)
    D_X_Z0_res$oos_fitted <- rep(0, nobs_Z0)
    D_X_Z0_res$auxilliary_fitted <-
      lapply(y_X_Z0_res$auxilliary_fitted, function (x) {x * 0})
    if (!silent) cat("E[D|Z=0,X]: perfect non-compliance -- Done! \n")
  } else {
    # Compute estimates of E[D|Z=0,X]
    D_X_Z0_res <- get_CEF(D[is_Z0], X[is_Z0, , drop = F],
                          learners = learners_DXZ,
                          ensemble_type = ensemble_type,
                          shortstack = shortstack,
                          custom_ensemble_weights = custom_ensemble_weights_DXZ,
                          cv_subsamples_list = cv_subsamples_list_Z0,
                          subsamples = subsamples_Z0,
                          silent = silent, progress = "E[D|Z=0,X]: ",
                          auxilliary_X = auxilliary_X_Z1)
  }#IFELSE

  # Check for perfect compliance
  if (all(D[Z==1] == 1)) {
    # Artificially construct values for subsample with Z=0
    D_X_Z1_res <- list(NULL)
    D_X_Z1_res$oos_fitted <- rep(0, nobs_Z1)
    D_X_Z1_res$auxilliary_fitted <-
      lapply(y_X_Z1_res$auxilliary_fitted, function (x) {x * 0})
    if (!silent) cat("E[D|Z=1,X]: perfect compliance -- Done! \n")
  } else {
    # Compute estimates of E[D|Z=1,X]
    D_X_Z1_res <- get_CEF(D[-is_Z0], X[-is_Z0, , drop = F],
                          learners = learners_DXZ,
                          ensemble_type = ensemble_type,
                          shortstack = shortstack,
                          custom_ensemble_weights = custom_ensemble_weights_DXZ,
                          cv_subsamples_list = cv_subsamples_list_Z1,
                          subsamples = subsamples_Z1,
                          silent = silent, progress = "E[D|Z=1,X]: ",
                          auxilliary_X = auxilliary_X_Z0)
  }#IFELSE

  # Compute estimates of E[Z|X]
  Z_X_res <- get_CEF(Z, X,
                     learners = learners_ZX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights_ZX,
                     cv_subsamples_list = cv_subsamples_list,
                     subsamples = subsamples,
                     compute_insample_predictions = F,
                     silent = silent, progress = "E[Z|X]: ")

  # Update ensemble type to account for (optional) custom weights
  ensemble_type <- dimnames(y_X_Z0_res$weights)[[2]]
  nensb <- ifelse(is.null(ensemble_type), 1, length(ensemble_type))

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # Construct reduced form variables
  l_Z0 <- l_Z1 <- p_Z0 <- p_Z1 <- matrix(0, nobs, nensb)
  l_Z0[is_Z0, ] <- y_X_Z0_res$oos_fitted
  l_Z1[-is_Z0, ] <- y_X_Z1_res$oos_fitted
  p_Z0[is_Z0, ] <- D_X_Z0_res$oos_fitted
  p_Z1[-is_Z0, ] <- D_X_Z1_res$oos_fitted
  if (!multiple_ensembles) {
    for (k in 1:sample_folds) {
      l_Z1[is_Z0][subsamples_Z0[[k]]] <- y_X_Z1_res$auxilliary_fitted[[k]]
      l_Z0[-is_Z0][subsamples_Z1[[k]]] <- y_X_Z0_res$auxilliary_fitted[[k]]
      p_Z1[is_Z0][subsamples_Z0[[k]]] <- D_X_Z1_res$auxilliary_fitted[[k]]
      p_Z0[-is_Z0][subsamples_Z1[[k]]] <- D_X_Z0_res$auxilliary_fitted[[k]]
    }#FOR
  } else {
    for (k in 1:sample_folds) {
      l_Z1[is_Z0, ][subsamples_Z0[[k]], ] <- y_X_Z1_res$auxilliary_fitted[[k]]
      l_Z0[-is_Z0, ][subsamples_Z1[[k]], ] <- y_X_Z0_res$auxilliary_fitted[[k]]
      p_Z1[is_Z0, ][subsamples_Z0[[k]], ] <- D_X_Z1_res$auxilliary_fitted[[k]]
      p_Z0[-is_Z0, ][subsamples_Z1[[k]], ] <- D_X_Z0_res$auxilliary_fitted[[k]]
    }#FOR
  }#IF
  r_X <- Z_X_res$oos_fitted

  # Trim propensity scores, return warnings
  r_X_tr <- trim_propensity_scores(r_X, trim, ensemble_type)

  # Compute the ATE using the constructed variables
  y_copy <- matrix(rep(y, nensb), nobs, nensb)
  D_copy <- matrix(rep(D, nensb), nobs, nensb)
  Z_copy <- matrix(rep(Z, nensb), nobs, nensb)
  psi_b <- Z_copy * (y_copy - l_Z1) / r_X_tr -
    (1 - Z_copy) * (y_copy - l_Z0) / (1 - r_X_tr) +
    l_Z1 - l_Z0
  psi_a <- -(Z_copy * (D_copy - p_Z1) / r_X_tr -
    (1 - Z_copy) * (D_copy - p_Z0) / (1 - r_X_tr) +
    p_Z1 - p_Z0)
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
  oos_pred <- list(EY_Z0_X = l_Z0, EY_Z1_X = l_Z1,
                   ED_Z0_X = p_Z0, ED_Z1_X = p_Z1,
                   EZ_X = r_X)

  # Organize output
  ddml_fit <- list(late = late, weights = weights, mspe = mspe,
                   psi_a = psi_a, psi_b = psi_b,
                   oos_pred = oos_pred,
                   learners = learners,
                   learners_DXZ = learners_DXZ,
                   learners_ZX = learners_ZX,
                   subsamples_Z0 = subsamples_Z0,
                   subsamples_Z1 = subsamples_Z1,
                   cv_subsamples_list_Z0 = cv_subsamples_list_Z0,
                   cv_subsamples_list_Z1 = cv_subsamples_list_Z1,
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
                                                     object$ensemble_type)
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

