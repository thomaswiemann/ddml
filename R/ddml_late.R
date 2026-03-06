#' Estimator of the Local Average Treatment Effect.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml()], [ddml::coef.ddml()],
#'     [ddml::confint.ddml()], [ddml::tidy.ddml()],
#'     [ddml::glance.ddml()], [ddml::diagnostics()]
#'
#' @description Estimator of the local average treatment effect.
#'
#' @details \code{ddml_late} provides a Double/Debiased Machine Learning
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
#' In this model, the local average treatment effect (LATE) is defined as
#'
#' \eqn{\theta_0^{\textrm{LATE}} \equiv
#'     E[g_0(1, X) - g_0(0, X)\vert p_0(1, X) > p_0(0, X)]}.
#'
#' The estimating equation is
#'     \eqn{E[m(W; \theta_0, \eta_0)] = 0}, where
#'     \eqn{m(W; \theta, \eta) = \psi_b(W; \eta) + \psi_a(W; \eta)\theta},
#'     \eqn{W = (Y, D, Z, X)}, and the Neyman orthogonal scores are
#'
#' \eqn{\psi_b(W; \eta) = \frac{Z(Y - \ell_1(X))}{m(X)} - \frac{(1-Z)(Y-\ell_0(X))}{1-m(X)} + \ell_1(X) - \ell_0(X)}
#'
#' \eqn{\psi_a(W; \eta) = -\left(\frac{Z(D - r_1(X))}{m(X)} - \frac{(1-Z)(D-r_0(X))}{1-m(X)} + r_1(X) - r_0(X)\right)}
#'
#'     with nuisance parameters \eqn{\eta = (\ell_0, \ell_1, r_0, r_1, m)} taking
#'     true values \eqn{\ell_{z,0}(X) = E[Y|Z=z, X]},
#'     \eqn{r_{z,0}(X) = E[D|Z=z, X]}, and \eqn{m_0(X) = E[Z|X]}.
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
#' @param splits An optional list of sample split objects. For
#'     \code{ddml_late}, recommended keys are \code{subsamples},
#'     \code{subsamples_byZ}, \code{cv_subsamples}, and
#'     \code{cv_subsamples_byZ}.
#' @param ... Deprecated arguments (\code{subsamples},
#'     \code{subsamples_byZ}, \code{cv_subsamples},
#'     \code{cv_subsamples_byZ}) are still accepted for backward
#'     compatibility but should be replaced with \code{splits}.
#'
#' @return \code{ddml_late} returns an object of S3 class
#'     \code{ddml_late}. An object of class \code{ddml_late} is a list
#'     containing the following components:
#'     \describe{
#'         \item{\code{coefficients}}{A matrix of estimated coefficients.}
#'         \item{\code{ensemble_weights}}{A list of matrices, providing the
#'             weight assigned to each base learner by the ensemble
#'             procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of
#'             each base learner computed by the cross-validation step
#'             in the ensemble construction.}
#'         \item{\code{r2}}{The out-of-sample R-squared.}
#'         \item{\code{psi_a}, \code{psi_b}}{Matrices needed for the
#'             computation of scores. Used in [ddml::summary.ddml()].}
#'         \item{\code{scores}}{A list of evaluated Neyman orthogonal
#'             scores.}
#'         \item{\code{J}}{A list of evaluated Jacobians.}
#'         \item{\code{fitted}}{A list of fitted nuisance estimators.
#'             See \code{\link{ddml_plm}} for more information.}
#'         \item{\code{splits}}{The data splitting structure.}
#'         \item{\code{learners},\code{learners_DXZ},
#'             \code{learners_ZX},\code{cluster_variable},
#'             \code{ensemble_type}}{Pass-through of selected
#'             user-provided arguments. See above.}
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
#' Imbens G, Angrist J (1994). "Identification and Estimation of Local Average
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
#' \donttest{
#' # Estimate the local average treatment effect using short-stacking with base
#' #     learners ols, lasso, and ridge. We can also use custom_ensemble_weights
#' #     to estimate the ATE using every individual base learner.
#' weights_everylearner <- diag(1, 3)
#' colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
#' late_fit <- ddml_late(y, D, Z, X,
#'                       learners = list(list(what = ols),
#'                                       list(what = mdl_glmnet),
#'                                       list(what = mdl_glmnet,
#'                                            args = list(alpha = 0))),
#'                       ensemble_type = 'nnls',
#'                       custom_ensemble_weights = weights_everylearner,
#'                       shortstack = TRUE,
#'                       sample_folds = 2,
#'                       silent = TRUE)
#' summary(late_fit)
#' }
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
                      stratify = TRUE,
                      trim = 0.01,
                      silent = FALSE,
                      parallel = NULL,
                      fitted = NULL,
                      splits = NULL,
                      save_crossval = TRUE,
                      ...) {
  cl <- match.call()

  # == Preliminaries ================================================

  dots <- list(...)
  messages <- resolve_messages(dots, "ddml_late", list(
    y_Z0 = "E[Y|Z=0,X]",
    y_Z1 = "E[Y|Z=1,X]",
    D_Z0 = "E[D|Z=0,X]",
    D_Z1 = "E[D|Z=1,X]",
    Z_X = "E[Z|X]"))

  validate_inputs(y = y, D = D, X = X, Z = Z,
                  learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type, trim = trim,
                  require_binary_D = FALSE)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DXZ,
                          learners_DXZ)
  validate_custom_weights(custom_ensemble_weights_ZX,
                          learners_ZX)

  nobs <- length(y)

  splits <- normalize_splits(
    splits = splits, by_label = "Z", ...)
  validate_fitted_splits_pair(fitted, splits, !shortstack)

  indxs <- get_sample_splits(
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    D = Z, stratify = stratify,
    subsamples = splits$subsamples,
    subsamples_byD = splits$subsamples_byZ,
    cv_subsamples = splits$cv_subsamples,
    cv_subsamples_byD = splits$cv_subsamples_byZ)
  check_subsamples(indxs$subsamples, indxs$subsamples_byD,
                   stratify, Z)

  t0 <- proc.time()[3]
  mode_str <- if (!is.null(parallel)) {
    p <- parse_parallel(parallel)
    paste0("parallel, ", p$num_cores, " cores")
  } else {
    "sequential"
  }#IFELSE
  if (!is.null(messages$start) && messages$start != "") {
    info_msg(sprintf(messages$start, mode_str),
             silent = silent)
  }#IF

  # == Reduced-form estimation ======================================

  # E[Z|X]
  Z_X_res <- get_CEF(Z, X,
                     learners = learners_ZX,
                     ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights =
                       custom_ensemble_weights_ZX,
                     subsamples = indxs$subsamples,
                     cv_subsamples = indxs$cv_subsamples,
                     compute_insample_predictions = FALSE,
                     silent = silent, label = messages$Z_X,
                     parallel = parallel,
                     fitted = fitted$Z_X)

  splits_Z <- list(
    subsamples = indxs$subsamples,
    subsamples_byD = indxs$subsamples_byD,
    cv_subsamples = indxs$cv_subsamples,
    cv_subsamples_byD = indxs$cv_subsamples_byD)

  fitted_Z_X <- build_fitted_entry(Z_X_res, save_crossval)

  # Second stage: Y ~ Z via ddml_ate
  ate_rf <- ddml_ate(
    y = y, D = Z, X = X,
    learners = learners, learners_DX = learners_ZX,
    sample_folds = sample_folds, cv_folds = cv_folds,
    custom_ensemble_weights = custom_ensemble_weights,
    custom_ensemble_weights_DX = custom_ensemble_weights_ZX,
    cluster_variable = cluster_variable,
    ensemble_type = ensemble_type, shortstack = shortstack,
    trim = trim, parallel = parallel, silent = silent,
    splits = splits_Z,
    fitted = list(y_X_D0 = fitted$y_X_Z0,
                  y_X_D1 = fitted$y_X_Z1,
                  D_X = fitted_Z_X),
    save_crossval = save_crossval,
    messages = list(start = "", finish = "",
                    y_D1 = messages$y_Z1,
                    y_D0 = messages$y_Z0, D_X = ""))

  # First stage: D ~ Z via ddml_ate
  ate_fs <- ddml_ate(
    y = D, D = Z, X = X,
    learners = learners_DXZ, learners_DX = learners_ZX,
    sample_folds = sample_folds, cv_folds = cv_folds,
    custom_ensemble_weights = custom_ensemble_weights_DXZ,
    custom_ensemble_weights_DX = custom_ensemble_weights_ZX,
    cluster_variable = cluster_variable,
    ensemble_type = ensemble_type, shortstack = shortstack,
    trim = trim, parallel = parallel, silent = silent,
    splits = splits_Z,
    fitted = list(y_X_D0 = fitted$D_X_Z0,
                  y_X_D1 = fitted$D_X_Z1,
                  D_X = fitted_Z_X),
    save_crossval = save_crossval,
    messages = list(start = "", finish = "",
                    y_D1 = messages$D_Z1,
                    y_D0 = messages$D_Z0, D_X = ""))

  ensemble_type <- ate_rf$ensemble_type
  nensb <- ncol(ate_rf$coefficients)

  # == Score construction ===========================================

  psi_b <- ate_rf$psi_b
  psi_a <- -ate_fs$psi_b

  # == Target parameter =============================================

  late <- as.vector(ate_rf$coefficients) / as.vector(ate_fs$coefficients)

  scores <- lapply(seq_len(nensb), function(j) {
    as.matrix(psi_a[, j] * late[j] + psi_b[, j])
  })
  J_list <- lapply(seq_len(nensb), function(j) {
    as.matrix(mean(psi_a[, j]))
  })

  coef_names <- "LATE"
  coef <- matrix(late, nrow = 1, ncol = nensb)
  rownames(coef) <- coef_names
  colnames(coef) <- ensemble_type

  # == Output =======================================================

  splits_export <- ate_rf$splits
  names(splits_export)[names(splits_export) ==
    "subsamples_byD"] <- "subsamples_byZ"
  if ("cv_subsamples_byD" %in% names(splits_export)) {
    names(splits_export)[names(splits_export) ==
      "cv_subsamples_byD"] <- "cv_subsamples_byZ"
  }#IF

  ddml_fit <- list(
    coefficients = coef,
    ensemble_weights = list(
      y_X_Z0 = ate_rf$ensemble_weights$y_X_D0,
      y_X_Z1 = ate_rf$ensemble_weights$y_X_D1,
      D_X_Z0 = ate_fs$ensemble_weights$y_X_D0,
      D_X_Z1 = ate_fs$ensemble_weights$y_X_D1,
      Z_X = ate_rf$ensemble_weights$D_X),
    mspe = list(
      y_X_Z0 = ate_rf$mspe$y_X_D0,
      y_X_Z1 = ate_rf$mspe$y_X_D1,
      D_X_Z0 = ate_fs$mspe$y_X_D0,
      D_X_Z1 = ate_fs$mspe$y_X_D1,
      Z_X = ate_rf$mspe$D_X),
    r2 = list(
      y_X_Z0 = ate_rf$r2$y_X_D0,
      y_X_Z1 = ate_rf$r2$y_X_D1,
      D_X_Z0 = ate_fs$r2$y_X_D0,
      D_X_Z1 = ate_fs$r2$y_X_D1,
      Z_X = ate_rf$r2$D_X),
    psi_a = psi_a, psi_b = psi_b,
    scores = scores, J = J_list,
    coef_names = coef_names,
    estimator_name = "Local Average Treatment Effect",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    learners = learners,
    learners_DXZ = learners_DXZ,
    learners_ZX = learners_ZX,
    cluster_variable = cluster_variable,
    fitted = list(
      y_X_Z0 = ate_rf$fitted$y_X_D0,
      y_X_Z1 = ate_rf$fitted$y_X_D1,
      D_X_Z0 = ate_fs$fitted$y_X_D0,
      D_X_Z1 = ate_fs$fitted$y_X_D1,
      Z_X = ate_rf$fitted$D_X),
    splits = splits_export,
    call = cl)

  elapsed <- round(proc.time()[3] - t0, 1)
  if (!is.null(messages$finish) && messages$finish != "") {
    info_msg(sprintf(messages$finish, elapsed),
             silent = silent)
  }#IF

  class(ddml_fit) <- c("ddml_late", "ddml")
  return(ddml_fit)
}#DDML_LATE
