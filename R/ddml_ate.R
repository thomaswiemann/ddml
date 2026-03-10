#' Estimator for the Average Treatment Effect
#'
#' @family ddml estimators
#'
#' @description Estimator for the average treatment effect and the average
#'     treatment effect on the treated.
#'
#' @details
#' \strong{Parameter of Interest:} \code{ddml_ate} and \code{ddml_att} provide Double/Debiased Machine
#'     Learning estimators for the average treatment effect and the average
#'     treatment effect on the treated, respectively, in the interactive model
#'     given by:
#'
#' \deqn{Y = g_0(D, X) + U,}
#'
#' where \eqn{(Y, D, X, U)} is a random vector such that
#'     \eqn{\operatorname{supp} D = \{0,1\}}, \eqn{E[U\vert D, X] = 0}, and
#'     \eqn{\Pr(D=1\vert X) \in (0, 1)} with probability 1,
#'     and \eqn{g_0} is an unknown nuisance function.
#'
#' In this model, the average treatment effect (ATE) is defined as
#'
#' \deqn{\theta_0^{\textrm{ATE}} \equiv E[g_0(1, X) - g_0(0, X)],}
#'
#' and the average treatment effect on the treated (ATT) is defined as
#'
#' \deqn{\theta_0^{\textrm{ATT}} \equiv E[g_0(1, X) - g_0(0, X)\vert D = 1].}
#'
#' \strong{Neyman Orthogonal Score:} The Neyman orthogonal scores are:
#'
#' \deqn{m^{\textrm{ATE}}(W; \theta, \eta) = \frac{D(Y - \ell_1(X))}{p(X)} - \frac{(1-D)(Y-\ell_0(X))}{1-p(X)} + \ell_1(X) - \ell_0(X) - \theta}
#'
#' \deqn{m^{\textrm{ATT}}(W; \theta, \eta) = \frac{D(Y - \ell_0(X))}{\pi} - \frac{p(X)(1-D)(Y-\ell_0(X))}{\pi(1-p(X))} - \frac{D\theta}{\pi}}
#'
#' where the nuisance parameters are \eqn{\eta = (\ell_0, \ell_1, p, \pi)} taking true values
#'     \eqn{\ell_{d,0}(X) = E[Y|D=d, X]}, \eqn{p_0(X) = E[D|X]}, and \eqn{\pi_0 = E[D]}.
#'
#' \strong{Linear Decomposition:} The scores decompose linearly in \eqn{\theta}:
#'
#' \deqn{m(W; \theta, \eta) = \psi_b(W; \eta) + \psi_a(W; \eta)\theta}
#'
#' where:
#'
#' \deqn{\psi_a^{\textrm{ATE}}(W; \eta) = -1}
#'
#' \deqn{\psi_b^{\textrm{ATE}}(W; \eta) = \frac{D(Y - \ell_1(X))}{p(X)} - \frac{(1-D)(Y-\ell_0(X))}{1-p(X)} + \ell_1(X) - \ell_0(X)}
#'
#' \deqn{\psi_a^{\textrm{ATT}}(W; \eta) = -\frac{D}{\pi}}
#'
#' \deqn{\psi_b^{\textrm{ATT}}(W; \eta) = \frac{D(Y - \ell_0(X))}{\pi} - \frac{p(X)(1-D)(Y-\ell_0(X))}{\pi(1-p(X))}}
#'
#' @inheritParams ddml-intro
#' @inheritParams ddml_apo
#' @param D The binary endogenous variable of interest.
#' @param splits An optional list of sample split objects. For
#'     \code{ddml_ate}/\code{ddml_att}, recommended keys are
#'     \code{subsamples}, \code{subsamples_byD}, \code{cv_subsamples},
#'     and \code{cv_subsamples_byD}.
#' @param ... Additional arguments passed to internal methods.
#'
#' @return \code{ddml_ate} and \code{ddml_att} return objects of S3
#'     class \code{ddml_ate}/\code{ddml_att} and \code{ddml}. See
#'     \code{\link{ddml-intro}} for the common output structure.
#'     Additional pass-through fields: \code{learners},
#'     \code{learners_DX}.
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
                  cluster_variable = cluster_variable,
                  require_binary_D = TRUE)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX, learners_DX)

  nobs <- length(y)

  validate_fitted_splits_pair(fitted, splits, !shortstack)

  indxs <- get_sample_splits(
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    D = D, stratify = stratify,
    subsamples = splits$D_X$subsamples,
    subsamples_byD = list(
      splits$y_X_D0$subsamples,
      splits$y_X_D1$subsamples),
    cv_subsamples = splits$D_X$cv_subsamples,
    cv_subsamples_byD = list(
      splits$y_X_D0$cv_subsamples,
      splits$y_X_D1$cv_subsamples))
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

  shared_splits <- list(subsamples = indxs$subsamples,
                        cv_subsamples = indxs$cv_subsamples)
  apo_splits_1 <- list(
    y_X = list(subsamples = indxs$subsamples_byD[[2]],
               cv_subsamples = indxs$cv_subsamples_byD[[2]]),
    D_X = shared_splits)
  apo_splits_0 <- list(
    y_X = list(subsamples = indxs$subsamples_byD[[1]],
               cv_subsamples = indxs$cv_subsamples_byD[[1]]),
    D_X = shared_splits)

  # Pre-ensembled propensity for ddml_apo delegation.
  # For d=0, flip: P(D=0|X) = 1 - P(D=1|X).
  fitted_D_X_1 <- list(cf_fitted = D_X_res$cf_fitted)
  fitted_D_X_0 <- list(cf_fitted = 1 - D_X_res$cf_fitted)

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
  psi_b <- lapply(seq_len(nensb), function(j) {
    apo_1$psi_b[[j]] - apo_0$psi_b[[j]]
  })
  psi_a <- lapply(seq_len(nensb), function(j) array(-1, dim = c(nobs, 1, 1)))

  # == Target parameter =============================================

  ate <- as.vector(apo_1$coefficients) - as.vector(apo_0$coefficients)

  scores <- array(NA_real_, dim = c(nobs, 1, nensb))
  J <- array(NA_real_, dim = c(1, 1, nensb))
  for (j in seq_len(nensb)) {
    scores[, 1, j] <- -ate[j] + as.vector(psi_b[[j]])
    J[1, 1, j] <- -1
  }#FOR

  coef_names <- "ATE"
  coef <- matrix(ate, nrow = 1, ncol = nensb)
  rownames(coef) <- coef_names
  colnames(coef) <- ensemble_type

  # == Output =======================================================

  elapsed <- round(proc.time()[3] - t0, 1)
  if (!is.null(messages$finish) && messages$finish != "") {
    info_msg(sprintf(messages$finish, elapsed),
             silent = silent)
  }#IF

  ddml(
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
    scores = scores, J = J,
    coef_names = coef_names,
    estimator_name = "Average Treatment Effect",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    cluster_variable = cluster_variable,
    fitted = list(
      y_X_D0 = apo_0$fitted$y_X,
      y_X_D1 = apo_1$fitted$y_X,
      D_X = build_fitted_entry(D_X_res, save_crossval)),
    splits = list(
      y_X_D0 = list(
        subsamples = indxs$subsamples_byD[[1]],
        cv_subsamples = indxs$cv_subsamples_byD[[1]]),
      y_X_D1 = list(
        subsamples = indxs$subsamples_byD[[2]],
        cv_subsamples = indxs$cv_subsamples_byD[[2]]),
      D_X = list(
        subsamples = indxs$subsamples,
        cv_subsamples = indxs$cv_subsamples)),
    call = cl,
    subclass = "ddml_ate",
    learners = learners,
    learners_DX = learners_DX
  )
}#DDML_ATE
