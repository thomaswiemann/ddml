#' Estimator for the Multi-Action Policy Value
#'
#' @family ddml estimators
#'
#' @description Estimator for the expected value of a multi-action
#'     policy, with optional per-level margins.
#'
#' @details
#' \strong{Parameter of Interest:} \code{ddml_policy} provides a
#'     Double/Debiased Machine Learning estimator for the expected
#'     value of a multi-action policy
#'     \eqn{\pi:\operatorname{supp}(X) \to \{d_1, \ldots, d_K\}}
#'     that assigns each unit to one of \eqn{K} treatment levels.
#'     The target parameter is
#'
#' \deqn{\theta_0 = \sum_{k=1}^{K}
#'     E\!\left[\omega_k(X)\, E[Y \mid D = d_k, X]\right],}
#'
#' where the known weight functions are
#'     \eqn{\omega_k(X) = c_k \,\mathbf{1}\{\pi(X) = d_k\}},
#'     and \eqn{c_1, \ldots, c_K} are user-supplied margins.
#'     When all margins equal one (\code{margins = NULL}), the
#'     parameter reduces to the policy value
#'     \eqn{E[Y(\pi(X))]}.
#'
#' Each term in the sum is a weighted average potential outcome
#'     (wAPO), estimated internally via \code{\link{ddml_apo}}.
#'
#' \strong{Nuisance Parameters:} For each treatment level
#'     \eqn{d_k}, the nuisance parameters are
#'     \eqn{\eta_k = (g_k, m_k)} taking true values
#'     \eqn{g_{k,0}(X) = E[Y \mid D = d_k, X]} and
#'     \eqn{m_{k,0}(X) = \Pr(D = d_k \mid X)}.
#'     Only \eqn{K-1} propensity models are estimated;
#'     the last is derived as
#'     \eqn{m_K(X) = 1 - \sum_{k=1}^{K-1} m_k(X)}.
#'
#' \strong{Neyman Orthogonal Score / Moment Equation:} The Neyman
#'     orthogonal score is:
#'
#' \deqn{m(W; \theta, \eta) = \sum_{k=1}^{K}
#'     \omega_k(X) \left[
#'     \frac{\mathbf{1}\{D = d_k\}\,(Y - g_k(X))}{m_k(X)}
#'     + g_k(X) \right] - \theta}
#'
#' \strong{Jacobian:}
#'
#' \deqn{J = -1}
#'
#' See \code{\link{ddml-intro}} for how the influence function
#' and inference are derived from these components.
#'
#' @inheritParams ddml-intro
#' @inheritParams ddml_apo
#' @param D The observed discrete (potentially multi-valued) treatment
#'     variable.
#' @param policy A vector of length \code{nobs} giving the
#'     policy-assigned treatment level for each unit. Values must
#'     be a subset of those observed in \code{D}.
#' @param margins An optional numeric vector of length \eqn{K}
#'     (the number of unique values in \code{policy}) giving
#'     per-level multipliers \eqn{c_k}. If \code{NULL} (the
#'     default), all margins are set to one, yielding the policy
#'     value \eqn{E[Y(\pi(X))]}.
#' @param splits An optional list of sample split objects. For
#'     \code{ddml_policy}, this is typically a named list keyed
#'     by treatment level. Typically obtained from a previous fit
#'     via \code{fit$splits}.
#' @param ... Additional arguments passed to internal methods.
#'
#' @return \code{ddml_policy} returns an object of S3 class
#'     \code{ddml_policy} and \code{ddml}. See
#'     \code{\link{ddml-intro}} for the common output structure.
#'     Additional pass-through fields: \code{learners},
#'     \code{learners_DX}, \code{policy}, \code{margins}.
#'
#' @export
#'
#' @references
#' Dudik M, Langford J, Li L (2011). "Doubly Robust Policy
#'     Evaluation and Learning." Proceedings of the 28th
#'     International Conference on Machine Learning, 1097-1104.
#'
#' Zhou Z, Athey S, Wager S (2023). "Offline Multi-Action Policy
#'     Learning: Generalization and Optimization." Operations
#'     Research, 71(2), 698-722.
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Define a simple policy: assign D=1 if age > median, else D=0
#' policy <- ifelse(X[, "age"] > median(X[, "age"]), 1, 0)
#'
#' # Estimate the policy value using a single base learner, ridge.
#' policy_fit <- ddml_policy(y, D, X,
#'                           policy = policy,
#'                           learners = list(what = mdl_glmnet),
#'                           sample_folds = 2,
#'                           silent = TRUE)
#' summary(policy_fit)
#'
ddml_policy <- function(y, D, X,
                        policy,
                        margins = NULL,
                        learners,
                        learners_DX = learners,
                        sample_folds = 10,
                        ensemble_type = "nnls",
                        shortstack = FALSE,
                        cv_folds = 10,
                        custom_ensemble_weights = NULL,
                        custom_ensemble_weights_DX =
                          custom_ensemble_weights,
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

  # Preliminaries --------------------------------------------------------------

  nobs <- length(y)
  d_levels <- sort(unique(policy))
  K <- length(d_levels)

  if (length(policy) != nobs) {
    stop("'policy' must have the same length as 'y'.",
         call. = FALSE)
  }#IF
  if (!all(d_levels %in% unique(D))) {
    stop("All values in 'policy' must appear in 'D'.",
         call. = FALSE)
  }#IF

  if (is.null(margins)) margins <- rep(1, K)
  if (!is.numeric(margins) || length(margins) != K) {
    stop("'margins' must be a numeric vector of length ", K,
         " (the number of unique policy levels).",
         call. = FALSE)
  }#IF

  exhaustive <- setequal(d_levels, sort(unique(D)))
  n_direct_prop <- if (exhaustive) K - 1L else K

  dots <- list(...)
  prop_levels <- d_levels[seq_len(n_direct_prop)]
  msg_labels <- stats::setNames(
    c(vapply(d_levels, function(d) paste0("E[Y|D=", d, ",X]"),
             ""),
      vapply(prop_levels,
             function(d) paste0("P(D=", d, "|X)"), "")),
    c(paste0("y_d", d_levels),
      paste0("D_d", prop_levels)))
  messages <- resolve_messages(dots, "ddml_policy",
                               as.list(msg_labels))

  validate_inputs(y = y, D = D, X = X,
                  learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type, trim = trim,
                  cluster_variable = cluster_variable,
                  require_binary_D = FALSE)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX,
                          learners_DX)

  validate_fitted_splits_pair(fitted, splits)

  # Shared main folds (stratified on multi-valued D)
  shared_splits <- if (!is.null(splits)) {
    splits$D_X
  } else {
    indxs <- get_sample_splits(
      cluster_variable = cluster_variable,
      sample_folds = sample_folds,
      cv_folds = cv_folds,
      D = D, stratify = stratify)
    list(subsamples = indxs$subsamples,
         cv_subsamples = indxs$cv_subsamples)
  }#IFELSE

  t0 <- proc.time()[3]
  announce_start(messages, parallel, silent)

  # Reduced-form estimation ----------------------------------------------------

  apo_list <- vector("list", K)

  # Levels with directly estimated propensity
  for (k in seq_len(n_direct_prop)) {
    d_k <- d_levels[k]
    apo_list[[k]] <- ddml_apo(
      y = y, D = D, X = X,
      d = d_k,
      weights = margins[k] * (policy == d_k),
      learners = learners, learners_DX = learners_DX,
      sample_folds = sample_folds, cv_folds = cv_folds,
      custom_ensemble_weights = custom_ensemble_weights,
      custom_ensemble_weights_DX =
        custom_ensemble_weights_DX,
      cluster_variable = cluster_variable,
      ensemble_type = ensemble_type,
      shortstack = shortstack, stratify = stratify,
      trim = trim, parallel = parallel, silent = silent,
      splits = list(D_X = shared_splits),
      fitted = list(
        y_X = fitted[[paste0("y_X_d", d_k)]],
        D_X = fitted[[paste0("D_X_d", d_k)]]),
      save_crossval = save_crossval,
      messages = list(
        start = "", finish = "",
        y_X = messages[[paste0("y_d", d_k)]],
        D_X = messages[[paste0("D_d", d_k)]]))
  }#FOR

  # Level K: derive propensity when exhaustive, estimate when not
  if (exhaustive) {
    derived_prop <- 1 - Reduce(
      "+", lapply(apo_list[seq_len(K - 1)],
                  function(a) a$fitted$D_X$cf_fitted))
    d_K <- d_levels[K]
    apo_list[[K]] <- ddml_apo(
      y = y, D = D, X = X,
      d = d_K,
      weights = margins[K] * (policy == d_K),
      learners = learners, learners_DX = learners_DX,
      sample_folds = sample_folds, cv_folds = cv_folds,
      custom_ensemble_weights = custom_ensemble_weights,
      custom_ensemble_weights_DX =
        custom_ensemble_weights_DX,
      cluster_variable = cluster_variable,
      ensemble_type = ensemble_type,
      shortstack = shortstack, stratify = stratify,
      trim = trim, parallel = parallel, silent = silent,
      splits = list(D_X = shared_splits),
      fitted = list(
        y_X = fitted[[paste0("y_X_d", d_K)]],
        D_X = list(cf_fitted = derived_prop)),
      save_crossval = save_crossval,
      messages = list(
        start = "", finish = "",
        y_X = messages[[paste0("y_d", d_K)]],
        D_X = ""))
  }#IF

  ensemble_type <- apo_list[[1]]$ensemble_type
  nensb <- ncol(apo_list[[1]]$coefficients)

  # Target parameter & influence function --------------------------------------

  pv <- rep(0, nensb)
  inf_func <- array(0, dim = c(nobs, 1, nensb))
  for (k in seq_len(K)) {
    pv <- pv + as.vector(apo_list[[k]]$coefficients)
    for (j in seq_len(nensb)) {
      inf_func[, 1, j] <- inf_func[, 1, j] +
        apo_list[[k]]$inf_func[, 1, j]
    }#FOR
  }#FOR

  scores <- array(NA_real_, dim = c(nobs, 1, nensb))
  J <- array(NA_real_, dim = c(1, 1, nensb))
  dinf_dtheta <- array(NA_real_, dim = c(nobs, 1, 1, nensb))
  for (j in seq_len(nensb)) {
    scores[, 1, j] <- inf_func[, 1, j]
    J[1, 1, j] <- -1
    dinf_dtheta[, 1, 1, j] <- -1
  }#FOR

  coef_names <- "Policy value"
  coef <- matrix(pv, nrow = 1, ncol = nensb)
  rownames(coef) <- coef_names
  colnames(coef) <- ensemble_type

  # Output ---------------------------------------------------------------------

  announce_finish(t0, messages, silent)

  ew <- mspe_out <- r2_out <- fit_out <- spl_out <- list()
  for (k in seq_len(K)) {
    ky <- paste0("y_X_d", d_levels[k])
    ew[[ky]] <- apo_list[[k]]$ensemble_weights$y_X
    mspe_out[[ky]] <- apo_list[[k]]$mspe$y_X
    r2_out[[ky]] <- apo_list[[k]]$r2$y_X
    fit_out[[ky]] <- apo_list[[k]]$fitted$y_X
    spl_out[[ky]] <- apo_list[[k]]$splits$y_X
  }#FOR
  for (k in seq_len(n_direct_prop)) {
    kd <- paste0("D_X_d", d_levels[k])
    ew[[kd]] <- apo_list[[k]]$ensemble_weights$D_X
    mspe_out[[kd]] <- apo_list[[k]]$mspe$D_X
    r2_out[[kd]] <- apo_list[[k]]$r2$D_X
    fit_out[[kd]] <- apo_list[[k]]$fitted$D_X
  }#FOR
  spl_out$D_X <- shared_splits

  ddml(
    coefficients = coef,
    ensemble_weights = ew,
    mspe = mspe_out,
    r2 = r2_out,
    inf_func = inf_func, dinf_dtheta = dinf_dtheta,
    scores = scores, J = J,
    coef_names = coef_names,
    estimator_name = "Multi-Action Policy Value",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    cluster_variable = cluster_variable,
    fitted = fit_out,
    splits = spl_out,
    call = cl,
    subclass = "ddml_policy",
    # ddml_policy-specific fields
    policy = policy,
    margins = margins,
    learners = learners,
    learners_DX = learners_DX
  )
}#DDML_POLICY
