#' Estimators of Average Treatment Effects.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml()], [ddml::coef.ddml()],
#'     [ddml::confint.ddml()], [ddml::tidy.ddml()],
#'     [ddml::glance.ddml()], [ddml::diagnostics()]
#'
#' @description Estimators of the average treatment effect and the average
#'     treatment effect on the treated.
#'
#' @details \code{ddml_ate} and \code{ddml_att} provide Double/Debiased Machine
#'     Learning estimators for the average treatment effect and the average
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
#' In this model, the average treatment effect (ATE) is defined as
#'
#' \eqn{\theta_0^{\textrm{ATE}} \equiv E[g_0(1, X) - g_0(0, X)]},
#'
#' and the average treatment effect on the treated (ATT) is defined as
#'
#' \eqn{\theta_0^{\textrm{ATT}} \equiv E[g_0(1, X) - g_0(0, X)\vert D = 1]}.
#'
#' The estimating equations are
#'     \eqn{E[m(W; \theta_0, \eta_0)] = 0}, where
#'     \eqn{m(W; \theta, \eta) = \psi_b(W; \eta) + \psi_a(W; \eta)\theta},
#'     \eqn{W = (Y, D, X)}, and the Neyman orthogonal scores are
#'
#' \eqn{\psi_b^{\textrm{ATE}}(W; \eta) = \frac{D(Y - \ell_1(X))}{m(X)} - \frac{(1-D)(Y-\ell_0(X))}{1-m(X)} + \ell_1(X) - \ell_0(X)}
#'
#' \eqn{\psi_a^{\textrm{ATE}}(W; \eta) = -1}
#'
#' \eqn{\psi_b^{\textrm{ATT}}(W; \eta) = \frac{D(Y - \ell_0(X))}{p} - \frac{m(X)(1-D)(Y-\ell_0(X))}{p(1-m(X))}}
#'
#' \eqn{\psi_a^{\textrm{ATT}}(W; \eta) = -\frac{D}{p}}
#'
#'     with nuisance parameters \eqn{\eta = (\ell_0, \ell_1, m, p)} taking
#'     true values \eqn{\ell_{d,0}(X) = E[Y|D=d, X]}, \eqn{m_0(X) = E[D|X]},
#'     and \eqn{p_0 = E[D]}.
#'
#' @inheritParams ddml_plm
#' @param D The binary endogenous variable of interest.
#' @param splits An optional list of sample split objects. For
#'     \code{ddml_ate}/\code{ddml_att}, recommended keys are
#'     \code{subsamples}, \code{subsamples_byD}, \code{cv_subsamples},
#'     and \code{cv_subsamples_byD}.
#' @param ... Deprecated arguments (\code{subsamples},
#'     \code{subsamples_byD}, \code{cv_subsamples},
#'     \code{cv_subsamples_byD}) are still accepted for backward
#'     compatibility but should be replaced with \code{splits}.
#' @param stratify Boolean for stratified cross-fitting: if \code{TRUE},
#'     subsamples are constructed to be balanced across treatment levels.
#' @param trim Number in (0, 1) for trimming the estimated propensity scores at
#'     \code{trim} and \code{1-trim}.
#' @param parallel An optional named list with parallel processing
#'     options. When \code{NULL} (the default), computation is
#'     sequential. Supported fields:
#'     \describe{
#'         \item{\code{cores}}{Number of cores to use.}
#'         \item{\code{export}}{Character vector of object names to
#'             export to parallel workers (for custom learners that
#'             reference global objects).}
#'         \item{\code{packages}}{Character vector of additional
#'             package names to load on workers (for custom learners
#'             that use packages not imported by \code{ddml}).}
#'     }
#' @param fitted An optional named list of per-equation cross-fitted
#'     predictions, typically obtained via \code{fit$fitted}. See
#'     \code{\link{ddml_plm}} for details and an example.
#' @param save_crossval Logical; store inner cross-validation
#'     residuals for exact weight recomputation on pass-through.
#'     See \code{\link{ddml_plm}} for details.
#'
#' @return \code{ddml_ate} and \code{ddml_att} return an object of S3 class
#'     \code{ddml_ate} and \code{ddml_att}, respectively. An object of class
#'     \code{ddml_ate} or \code{ddml_att} is a list containing
#'     the following components:
#'     \describe{
#'         \item{\code{coefficients}}{A matrix of estimated coefficients.}
#'         \item{\code{ensemble_weights}}{A list of matrices, providing the
#'             weight assigned to each base learner by the ensemble
#'             procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner computed by the cross-validation step in the
#'             ensemble construction.}
#'         \item{\code{r2}}{The out-of-sample R-squared.}
#'         \item{\code{psi_a}, \code{psi_b}}{Matrices needed for the
#'             computation of scores. Used in [ddml::summary.ddml()].}
#'         \item{\code{scores}}{A list of evaluated Neyman orthogonal
#'             scores.}
#'         \item{\code{J}}{A list of evaluated Jacobians.}
#'         \item{\code{fitted}}{A list of fitted nuisance estimators.
#'             See \code{\link{ddml_plm}} for more information.}
#'         \item{\code{splits}}{The data splitting structure.}
#'         \item{\code{learners},\code{learners_DX},
#'             \code{cluster_variable},
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
#'                     learners = list(list(what = ols),
#'                                     list(what = mdl_glmnet),
#'                                     list(what = mdl_glmnet,
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
  messages <- resolve_messages(dots, "ddml_ate", list(
    y_D0 = "E[Y|D=0,X]",
    y_D1 = "E[Y|D=1,X]",
    D_X = "E[D|X]"))

  validate_inputs(y = y, D = D, X = X, learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type, trim = trim,
                  require_binary_D = TRUE)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX, learners_DX)

  nobs <- length(y)

  splits <- normalize_splits(splits = splits, by_label = "D", ...)
  validate_fitted_splits_pair(fitted, splits, !shortstack)

  indxs <- get_sample_splits(
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    D = D, stratify = stratify,
    subsamples = splits$subsamples,
    subsamples_byD = splits$subsamples_byD,
    cv_subsamples = splits$cv_subsamples,
    cv_subsamples_byD = splits$cv_subsamples_byD)
  check_subsamples(indxs$subsamples, indxs$subsamples_byD,
                   stratify, D)

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

  # E[D|X]
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX,
                     ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights =
                       custom_ensemble_weights_DX,
                     subsamples = indxs$subsamples,
                     cv_subsamples = indxs$cv_subsamples,
                     silent = silent, label = messages$D_X,
                     parallel = parallel,
                     fitted = fitted$D_X)

  # Splits for apo(d=1): byD indices match D_ind = D
  apo_splits_1 <- list(
    subsamples = indxs$subsamples,
    subsamples_byd = indxs$subsamples_byD,
    cv_subsamples = indxs$cv_subsamples,
    cv_subsamples_byd = indxs$cv_subsamples_byD)
  # Splits for apo(d=0): swap byD indices since D_ind = 1-D
  apo_splits_0 <- list(
    subsamples = indxs$subsamples,
    subsamples_byd = list(indxs$subsamples_byD[[2]],
                          indxs$subsamples_byD[[1]]),
    cv_subsamples = indxs$cv_subsamples,
    cv_subsamples_byd = if (!is.null(indxs$cv_subsamples_byD))
      list(indxs$cv_subsamples_byD[[2]],
           indxs$cv_subsamples_byD[[1]]))

  # Pre-ensembled propensity for ddml_apo delegation.
  # For d=0, flip: P(D=0|X) = 1 - P(D=1|X).
  fitted_D_X_1 <- list(ensemble_fitted = D_X_res$oos_fitted)
  fitted_D_X_0 <- list(ensemble_fitted = 1 - D_X_res$oos_fitted)

  # E[g(1,X)] via ddml_apo
  apo_1 <- ddml_apo(
    y = y, D = D, X = X, d = 1, weights = NULL,
    learners = learners, learners_DX = learners_DX,
    sample_folds = sample_folds, cv_folds = cv_folds,
    custom_ensemble_weights = custom_ensemble_weights,
    custom_ensemble_weights_DX = custom_ensemble_weights_DX,
    cluster_variable = cluster_variable,
    ensemble_type = ensemble_type, shortstack = shortstack,
    trim = trim, parallel = parallel, silent = silent,
    splits = apo_splits_1,
    fitted = list(y_X = fitted$y_X_D1,
                  D_X = fitted_D_X_1),
    messages = list(start = "", finish = "",
                    y_X = messages$y_D1, D_X = ""))

  # E[g(0,X)] via ddml_apo
  apo_0 <- ddml_apo(
    y = y, D = D, X = X, d = 0, weights = NULL,
    learners = learners, learners_DX = learners_DX,
    sample_folds = sample_folds, cv_folds = cv_folds,
    custom_ensemble_weights = custom_ensemble_weights,
    custom_ensemble_weights_DX = custom_ensemble_weights_DX,
    cluster_variable = cluster_variable,
    ensemble_type = ensemble_type, shortstack = shortstack,
    trim = trim, parallel = parallel, silent = silent,
    splits = apo_splits_0,
    fitted = list(y_X = fitted$y_X_D0,
                  D_X = fitted_D_X_0),
    messages = list(start = "", finish = "",
                    y_X = messages$y_D0, D_X = ""))

  ensemble_type <- apo_1$ensemble_type
  nensb <- ncol(apo_1$coefficients)

  # == Score construction ===========================================
  psi_b <- apo_1$psi_b - apo_0$psi_b
  psi_a <- matrix(-1, nobs, nensb)

  # == Target parameter =============================================

  ate <- as.vector(apo_1$coefficients) - as.vector(apo_0$coefficients)

  scores <- lapply(seq_len(nensb), function(j) {
    as.matrix(psi_a[, j] * ate[j] + psi_b[, j])
  })
  J_list <- lapply(seq_len(nensb), function(j) {
    as.matrix(mean(psi_a[, j]))
  })

  coef_names <- "ATE"
  coef <- matrix(ate, nrow = 1, ncol = nensb)
  rownames(coef) <- coef_names
  colnames(coef) <- ensemble_type

  # == Output =======================================================

  ddml_fit <- list(
    coefficients = coef,
    ensemble_weights = list(
      y_X_D0 = apo_0$ensemble_weights$y_X,
      y_X_D1 = apo_1$ensemble_weights$y_X,
      D_X = D_X_res$weights),
    mspe = list(y_X_D0 = apo_0$mspe$y_X,
                y_X_D1 = apo_1$mspe$y_X,
                D_X = D_X_res$mspe),
    r2 = list(y_X_D0 = apo_0$r2$y_X,
              y_X_D1 = apo_1$r2$y_X,
              D_X = D_X_res$r2),
    psi_a = psi_a, psi_b = psi_b,
    scores = scores, J = J_list,
    coef_names = coef_names,
    estimator_name = "Average Treatment Effect",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    learners = learners,
    learners_DX = learners_DX,
    cluster_variable = cluster_variable,
    fitted = list(
      y_X_D0 = apo_0$fitted$y_X,
      y_X_D1 = apo_1$fitted$y_X,
      D_X = build_fitted_entry(D_X_res, save_crossval)),
    splits = list(
      subsamples = indxs$subsamples,
      subsamples_byD = indxs$subsamples_byD,
      cv_subsamples = indxs$cv_subsamples,
      cv_subsamples_byD = indxs$cv_subsamples_byD),
    call = cl)

  elapsed <- round(proc.time()[3] - t0, 1)
  if (!is.null(messages$finish) && messages$finish != "") {
    info_msg(sprintf(messages$finish, elapsed),
             silent = silent)
  }#IF

  class(ddml_fit) <- c("ddml_ate", "ddml")
  return(ddml_fit)
}#DDML_ATE
