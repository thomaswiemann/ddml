#' DDML Estimator for the Average Potential Outcome
#'
#' @description \code{ddml_apo} provides a Double/Debiased Machine Learning
#'     estimator for the average potential outcome, allowing for custom 
#'     weights \eqn{\omega(X)}.
#'
#' @details \code{ddml_apo} provides a Double/Debiased Machine Learning
#'     estimator for the target parameter \eqn{\theta_0} in the model given by
#'
#' \eqn{Y = g_0(D, X) + U,}
#'
#' where \eqn{(Y, D, X, U)} is a random vector such that
#'     \eqn{E[U\vert D, X] = 0} and
#'     \eqn{\Pr(D=d\vert X) \in (0, 1)} with probability 1,
#'     and \eqn{g_0} is an unknown nuisance function.
#'
#' In this model, the average potential outcome (APO) for treatment 
#'     level \eqn{d} is defined as
#'
#' \eqn{\theta_0^{\textrm{APO}} \equiv E[\omega(X) g_0(d, X)]},
#'
#' where \eqn{\omega(X)} is a known weighting function. If \eqn{\omega(X) = 1},
#'     this parameter corresponds to the standard Average Potential Outcome (APO)
#'     at treatment level \eqn{d}.
#'
#' The estimating equation is
#'     \eqn{E[m(W; \theta_0, \eta_0)] = 0}, where
#'     \eqn{W = (Y, D, X)}, and the Neyman orthogonal score is
#'
#' \eqn{m(W; \theta, \eta) = \left( \frac{\mathbbm{1}\{D=d\} (Y - \ell(X))}{p(X)} + \ell(X) \right) \omega(X) - \theta}
#'
#'     with nuisance parameters \eqn{\eta = (\ell, p)} taking
#'     true values \eqn{\ell_0(X) = E[Y|D=d, X]} and \eqn{p_0(X) = E[\mathbbm{1}\{D=d\}|X]}.
#'
#' @inheritParams ddml_plm
#' @inheritParams ddml_ate
#' @param D The endogenous variable of interest. Can be discrete or continuous.
#' @param d The treatment level of interest. The default is \code{d = 1}.
#' @param weights A numeric vector of length \code{nobs} specifying the weights 
#'     \eqn{\omega(X)}. If \code{weights = NULL} (the default), a vector of 1s 
#'     is used, which estimates the Average Potential Outcome (APO).
#' @param splits An optional list of sample split objects. For
#'     \code{ddml_apo}, this must be a list with elements \code{subsamples} and
#'     \code{cv_subsamples} (and optionally \code{subsamples_byd} and
#'     \code{cv_subsamples_byd} for stratified splitting). Typically
#'     obtained from a previous fit via \code{fit$splits}.
#' @param ... Deprecated arguments (\code{subsamples_byd},
#'     \code{cv_subsamples_byd}) are still accepted for backward
#'     compatibility but should be replaced with the \code{splits}
#'     argument.
#'
#' @return \code{ddml_apo} returns a list object of class \code{ddml_apo} 
#'     and \code{ddml} with the following properties:
#' \describe{
#'   \item{\code{coefficients}}{A matrix of estimated coefficients.}
#'   \item{\code{ensemble_weights}}{The ensemble weights.}
#'   \item{\code{mspe}}{The in-sample mean squared prediction errors.}
#'   \item{\code{r2}}{The out-of-sample R-squared.}
#'   \item{\code{psi_a}, \code{psi_b}}{The scale and scores associated with
#'     evaluating the Neyman orthogonal score.}
#'   \item{\code{scores}}{A list of evaluating the Neyman orthogonal scores.}
#'   \item{\code{J}}{A list of evaluating the Jacobian of the Neyman
#'     orthogonal scores.}
#'   \item{\code{fitted}}{An optional list of fitted nuisance estimators.
#'     See \code{\link{ddml_plm}} for more information.}
#'   \item{\code{splits}}{The data splitting structure.}
#'   \item{\code{learners}}{The original list of learners.}
#'   \item{\code{learners_DX}}{The original list of learners for the propensity
#'     score.}
#'   \item{\code{ensemble_type}}{The ensemble type.}
#'   \item{\code{nobs}}{The number of observations.}
#'   \item{\code{sample_folds}}{The number of sample folds.}
#'   \item{\code{cv_folds}}{The number of cross-validation folds.}
#'   \item{\code{shortstack}}{A boolean indicating if shortstacking was used.}
#'   \item{\code{cluster_variable}}{The original cluster variable.}
#'   \item{\code{call}}{The match.call().}
#' }
#'     
#' @references Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). 
#'     "Model Averaging and Double Machine Learning." 
#'     \url{https://arxiv.org/abs/2401.01645}
#' @seealso \code{\link{summary.ddml}}
#'     \code{\link{ddml_plm}}, \code{\link{ddml_pliv}},
#'     \code{\link{ddml_fpliv}}, \code{\link{ddml_late}},
#'     \code{\link{ddml_ate}}
#' @export
#'
ddml_apo <- function(y, D, X,
                     d = 1,
                     weights = NULL,
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
  messages <- resolve_messages(dots, "ddml_apo", list(
    y_X = paste0("E[Y|D=", d, ",X]"),
    D_X = paste0("P(D=", d, "|X)")))

  D_ind <- 1 * (D == d)

  validate_inputs(y = y, D = D_ind, X = X,
                  learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type, trim = trim,
                  weights = weights,
                  require_binary_D = TRUE)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX,
                          learners_DX)

  nobs <- length(y)
  is_d <- which(D == d)
  if (is.null(weights)) weights <- rep(1, nobs)

  # Construct sample splits for cross-fitting and cross-validation
  splits <- normalize_splits(splits = splits, by_label = "d", ...)
  validate_fitted_splits_pair(fitted, splits, !shortstack)
  indxs <- get_sample_splits(
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    D = D_ind, stratify = stratify,
    subsamples = splits$subsamples,
    subsamples_byD = splits$subsamples_byd,
    cv_subsamples = splits$cv_subsamples,
    cv_subsamples_byD = splits$cv_subsamples_byd)
  check_subsamples(indxs$subsamples, indxs$subsamples_byD,
                   stratify, D_ind)

  # Prepare paralllelization and estimation start
  t0 <- proc.time()[3]
  mode_str <- if (!is.null(parallel)) {
    p <- parse_parallel(parallel)
    paste0("parallel, ", p$num_cores, " cores")
  } else {
    "sequential"
  }#IFELSE
  if (!is.null(messages$start) && messages$start != "") {
    info_msg(sprintf(messages$start, mode_str), silent = silent)
  }#IF

  # == Reduced-form estimation ======================================

  # P(D = d | X)
  D_X_res <- get_CEF(D_ind, X,
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

  # E[Y | D = d, X] (subsamples_byD[[2]] = D_ind == 1)
  y_X_res <- get_CEF(y[is_d], X[is_d, , drop = FALSE],
                     learners = learners,
                     ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights =
                       custom_ensemble_weights,
                     subsamples = indxs$subsamples_byD[[2]],
                     cv_subsamples =
                       indxs$cv_subsamples_byD[[2]],
                     silent = silent, label = messages$y_X,
                     auxiliary_X = get_auxiliary_X(
                       indxs$aux_indx[[2]], X),
                     parallel = parallel,
                     fitted = fitted$y_X)

  ensb_info <- update_ensemble_info(y_X_res$weights,
                                    D_X_res$oos_fitted)
  ensemble_type <- ensb_info$ensemble_type
  nensb <- ensb_info$nensb

  # == Score construction ===========================================

  g_X <- extrapolate_CEF(
    D = D_ind,
    CEF_res_byD = list(list(
      fit = list(
        ensemble_fitted = y_X_res$oos_fitted,
        auxiliary_fitted = y_X_res$auxiliary_fitted),
      d = 1)),
    aux_indx = indxs$aux_indx[2])[, , 1]
  m_X <- D_X_res$oos_fitted
  is_internal <- is.null(messages$start) ||
    messages$start == ""
  m_X_tr <- trim_propensity_scores(m_X, trim, ensemble_type,
                                   silent = is_internal)

  weights_mat <- matrix(weights, nobs, nensb)
  D_ind_mat <- matrix(D_ind, nobs, nensb)
  y_mat <- matrix(y, nobs, nensb)

  psi_b <- (D_ind_mat * (y_mat - g_X) / m_X_tr + g_X) * weights_mat
  psi_a <- matrix(-1, nobs, nensb)

  # == Target parameter =============================================

  apo <- colMeans(psi_b)

  scores <- lapply(seq_len(nensb), function(j) {
    as.matrix(psi_a[, j] * apo[j] + psi_b[, j])
  })
  J_list <- lapply(seq_len(nensb), function(j) {
    as.matrix(mean(psi_a[, j]))
  })

  coef_names <- "APO"
  coef <- matrix(apo, nrow = 1, ncol = nensb)
  rownames(coef) <- coef_names
  colnames(coef) <- ensemble_type

  # == Output =======================================================

  ddml_fit <- list(
    coefficients = coef,
    ensemble_weights = list(y_X = y_X_res$weights,
                            D_X = D_X_res$weights),
    mspe = list(y_X = y_X_res$mspe,
                D_X = D_X_res$mspe),
    r2 = list(y_X = y_X_res$r2,
              D_X = D_X_res$r2),
    psi_a = psi_a, psi_b = psi_b,
    scores = scores, J = J_list,
    coef_names = coef_names,
    estimator_name = "Average Potential Outcome",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    d = d,
    weights = weights,
    learners = learners,
    learners_DX = learners_DX,
    cluster_variable = cluster_variable,
    fitted = list(
      y_X = build_fitted_entry(y_X_res, save_crossval,
                               include_auxiliary = TRUE),
      D_X = build_fitted_entry(D_X_res, save_crossval)),
    splits = list(
      subsamples = indxs$subsamples,
      subsamples_byd = indxs$subsamples_byD,
      cv_subsamples = indxs$cv_subsamples,
      cv_subsamples_byd = indxs$cv_subsamples_byD),
    call = cl)

  elapsed <- round(proc.time()[3] - t0, 1)
  if (!is.null(messages$finish) && messages$finish != "") {
    info_msg(sprintf(messages$finish, elapsed),
             silent = silent)
  }#IF

  class(ddml_fit) <- c("ddml_apo", "ddml")
  return(ddml_fit)
}#DDML_APO

#' @export
summary.ddml_apo <- function(object, ...) {
  summary.ddml(object, ...)
}
