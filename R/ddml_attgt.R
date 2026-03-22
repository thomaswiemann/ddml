#' Estimator for Group-Time Average Treatment Effects
#'
#' @family ddml estimators
#'
#' @description Estimator for group-time average treatment effects
#'     on the treated (GT-ATT) in staggered Difference-in-Differences
#'     designs, using cross-fitted AIPW scores.
#'
#' @details
#' \strong{Parameter of Interest:} \code{ddml_attgt} provides a
#'     Double/Debiased Machine Learning estimator for the group-time
#'     average treatment effects on the treated (GT-ATT) in the
#'     staggered adoption model given by:
#'
#' \deqn{\theta_0^{(g,t)} = E[\Delta_g Y_{i,t} | G_i = g]
#'     - E[E[\Delta_g Y_{i,t} | X_i, G_i \ne g, G_i > t] | G_i = g]}
#'
#' where \eqn{\Delta_g Y_{i,t} = Y_{i,t} - Y_{i,g^*}} is the
#'     difference relative to the universal base period
#'     \eqn{g^* = \max\{s : s < g - \text{anticipation}\}}.
#'
#' \strong{Neyman Orthogonal Score:} For each cell \eqn{(g,t)}, the
#'     AIPW score is:
#'
#' \deqn{m^{(g,t)}(W_i; \theta, \eta) =
#'     \frac{D_i(\Delta_g Y_{i,t} - \ell(X_i))}{\pi}
#'     - \frac{q(X_i)(1-D_i)(\Delta_g Y_{i,t}
#'     - \ell(X_i))}{\pi(1 - q(X_i))}
#'     - \frac{D_i}{\pi} \theta}
#'
#' where \eqn{D_i = \mathbb{1}\{G_i = g\}} is the cell-level
#'     treatment indicator, and the nuisance parameters are
#'     \eqn{\eta = (q, \ell, \pi)} taking true values
#'     \eqn{q_0(X) = \Pr(G_i = g | X_i, \{G_i = g\} \cup
#'     \{G_i > t\})}, \eqn{\ell_0(X) = E[\Delta_g Y_{i,t} |
#'     G_i \ne g, G_i > t, X_i]}, and
#'     \eqn{\pi_0 = \Pr(G_i = g)}.
#'
#' \strong{Jacobian:}
#'
#' \deqn{J = -E[D / \pi]}
#'
#' See \code{\link{ddml-intro}} for how the influence function
#' and inference are derived from these components.
#'
#' @inheritParams ddml-intro
#' @param y An \eqn{n \times T} numeric matrix of outcomes.
#'     Row \eqn{i} corresponds to unit \eqn{i}, column \eqn{j}
#'     to time period \code{t[j]}.
#' @param X An \eqn{n \times p} matrix of time-invariant covariates,
#'     or \code{NULL}.
#' @param t A numeric vector of length \eqn{T} giving the time
#'     period labels (must match columns of \code{y}).
#' @param G A numeric vector of length \eqn{n}. Entry \eqn{i} is the
#'     first treatment period for unit \eqn{i}. Use \code{0} or
#'     \code{Inf} for never-treated units.
#' @param learners_qX Optional argument to allow for different
#'     estimators of the cell-level propensity score
#'     \eqn{q^{(g,t)}(X)}. Setup is identical to
#'     \code{learners}.
#' @param custom_ensemble_weights_qX Optional argument to allow for
#'     different custom ensemble weights for \code{learners_qX}.
#'     Setup is identical to \code{custom_ensemble_weights}.
#' @param trim Number in (0, 1) for trimming the estimated
#'     propensity scores at \code{trim} and \code{1-trim}.
#' @param control_group Character. \code{"notyettreated"} (default)
#'     uses never-treated and not-yet-treated units as controls.
#'     \code{"nevertreated"} uses only never-treated units.
#' @param anticipation Non-negative integer. Number of periods before
#'     treatment where anticipation effects may occur. Default 0.
#'
#' @return \code{ddml_attgt} returns an object of S3 class
#'     \code{ddml_attgt} and \code{ddml}. See \code{\link{ddml-intro}}
#'     for the common output structure. Additional pass-through
#'     fields: \code{learners}, \code{learners_qX},
#'     \code{cell_info}, \code{control_group}, \code{anticipation}.
#'
#' @references
#' Callaway B, Sant'Anna P H C (2021). "Difference-in-Differences
#' with multiple time periods." Journal of Econometrics,
#' 225(2), 200-230.
#'
#' Chang N-C (2020). "Double/debiased machine learning for
#' difference-in-differences models." Econometrics Journal,
#' 23(2), 177-191.
#'
#' Ahrens A, Chernozhukov V, Hansen C B, Kozbur D, Schaffer M E,
#' Wiemann T (2026). "An Introduction to Double/Debiased Machine
#' Learning." Journal of Economic Literature, forthcoming.
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 200; T_ <- 4
#' X <- matrix(rnorm(n * 2), n, 2)
#' G <- sample(c(3, 4, Inf), n, replace = TRUE,
#'             prob = c(0.3, 0.3, 0.4))
#' y <- matrix(rnorm(n * T_), n, T_)
#' # Add treatment effect for treated units
#' for (i in seq_len(n)) {
#'   if (is.finite(G[i])) {
#'     for (j in seq_len(T_)) {
#'       if (j >= G[i]) y[i, j] <- y[i, j] + 1
#'     }
#'   }
#' }
#' fit <- ddml_attgt(y, X, t = 1:T_, G = G,
#'                 learners = list(what = ols),
#'                 sample_folds = 2,
#'                 silent = TRUE)
#' summary(fit)
#' }
ddml_attgt <- function(y, X = NULL, t, G,
                     learners,
                     learners_qX = learners,
                     sample_folds = 10,
                     ensemble_type = "nnls",
                     shortstack = FALSE,
                     cv_folds = 10,
                     custom_ensemble_weights = NULL,
                     custom_ensemble_weights_qX =
                       custom_ensemble_weights,
                     cluster_variable = seq_len(nrow(as.matrix(y))),
                     trim = 0.01,
                     control_group = c("notyettreated",
                                       "nevertreated"),
                     anticipation = 0,
                     silent = FALSE,
                     parallel = NULL,
                     fitted = NULL,
                     splits = NULL,
                     save_crossval = TRUE,
                     ...) {
  cl <- match.call()

  # Preliminaries --------------------------------------------------------------

  control_group <- match.arg(control_group)
  y <- as.matrix(y)
  n <- nrow(y)
  T_ <- ncol(y)
  stopifnot(length(t) == T_, length(G) == n)
  if (!is.null(X)) {
    X <- as.matrix(X)
    stopifnot(nrow(X) == n)
  }#IF

  dots <- list(...)
  messages <- resolve_messages(dots, "ddml_attgt", list())

  validate_inputs(learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type,
                  trim = trim,
                  cluster_variable = cluster_variable)
  validate_inputs(learners = learners_qX)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_qX, learners_qX)
  validate_fitted_splits_pair(fitted, splits)

  t0 <- proc.time()[3]
  announce_start(messages, parallel = parallel, silent = silent)

  # Global fold assignment -----------------------------------------------------

  never_treated <- (G == 0 | is.infinite(G) | G > max(t))
  ever_treated <- as.integer(!never_treated)

  # On re-entry (pass-through): recover global subsamples from the
  # dedicated "global" key stored by the previous fit.
  global_subsamples_init <- if (!is.null(splits)) splits[["global"]]$subsamples

  global_indxs <- get_sample_splits(
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    D = G, stratify = TRUE,
    subsamples = global_subsamples_init)
  global_subsamples <- global_indxs$subsamples

  # Convert to fold-vector for fast projection
  global_fold_vec <- integer(n)
  for (k in seq_len(sample_folds)) global_fold_vec[global_subsamples[[k]]] <- k

  # Enumerate (g, t) pairs -----------------------------------------------------

  groups <- sort(unique(G[!never_treated]))
  times <- sort(t)
  gt_list <- list()
  for (g_val in groups) {
    # Universal base period: last time before g - anticipation
    eligible_base <- times[times < g_val - anticipation]
    if (length(eligible_base) == 0) next
    base_t <- max(eligible_base)
    base_col <- match(base_t, t)
    for (tt in times) {
      if (tt == base_t) next
      # Check control pool size before adding the cell
      if (control_group == "nevertreated") {
        n_ctrl <- sum(never_treated)
      } else {
        cutoff <- max(tt, base_t + anticipation)
        n_ctrl <- sum((never_treated | (G > cutoff)) & G != g_val)
      }#IFELSE
      n_treat <- sum(G == g_val)
      if (n_ctrl == 0 || n_treat == 0) next
      gt_list <- c(gt_list, list(list(
        g = g_val, t = tt,
        base_t = base_t,
        base_col = base_col,
        tt_col = match(tt, t))))
    }#FOR
  }#FOR
  C <- length(gt_list)
  if (C == 0) {
    stop("No valid (g,t) cells found. Check that 'G' and ",
         "'t' define at least one post-treatment period.", call. = FALSE)
  }#IF


  # Pre-allocate ---------------------------------------------------------------

  nensb <- NULL
  coef_mat <- NULL
  scores_arr <- NULL
  J_arr <- NULL
  psi_b_full <- NULL
  psi_a_full <- NULL

  # Flat diagnostics lists — keyed by "ATT(g,t):equation"
  all_ensemble_weights <- list()
  all_mspe <- list()
  all_r2 <- list()
  all_fitted <- list()
  all_splits <- list()
  cell_info <- data.frame(
    group = integer(C), time = integer(C),
    base_period = integer(C),
    n_treated = integer(C), n_control = integer(C))

  # Main loop: per-cell ATT estimation -----------------------------------------

  for (idx in seq_len(C)) {
    gtp <- gt_list[[idx]]
    g_val <- gtp$g; tt <- gtp$t
    cell_prefix <- paste0("ATT(", g_val, ",", tt, ")")

    # Identify treated and control units
    treated <- (G == g_val)
    if (control_group == "nevertreated") {
      control <- never_treated
    } else {
      cutoff <- max(tt, gtp$base_t + anticipation)
      control <- never_treated | (G > cutoff)
    }#IFELSE
    keep <- which(treated | control)
    n_cell <- length(keep)
    cell_info[idx, ] <- list(
      g_val, tt, gtp$base_t,
      sum(treated[keep]),
      sum(!treated[keep]))

    info_msg(sprintf(
      "  Cell (%s,%s): n_treat=%d, n_ctrl=%d [%d/%d]",
      g_val, tt, cell_info$n_treated[idx],
      cell_info$n_control[idx], idx, C),
      silent = silent)

    # Difference outcomes
    delta_y <- y[keep, gtp$tt_col] - y[keep, gtp$base_col]
    D_cell <- as.integer(treated[keep])
    X_cell <- if (!is.null(X)) {
      X[keep, , drop = FALSE]
    } else {
      matrix(1, nrow = n_cell, ncol = 1)
    }#IFELSE

    # Project global folds to cell subset
    cell_fold_vec <- global_fold_vec[keep]
    active_folds <- sort(unique(cell_fold_vec))
    cell_subsamples <- lapply(active_folds, function(k) {
      which(cell_fold_vec == k)
    })

    # Reconstruct per-cell fitted/splits on re-entry
    cell_fitted <- if (!is.null(fitted)) {
      list(
        y_X_D0 = fitted[[paste0(cell_prefix, ":y_X_D0")]],
        D_X    = fitted[[paste0(cell_prefix, ":D_X")]],
        D      = fitted[[paste0(cell_prefix, ":D")]])
    }#IF
    cell_splits <- if (!is.null(splits)) {
      list(
        y_X_D0 = splits[[paste0(cell_prefix, ":y_X_D0")]],
        y_X_D1 = splits[[paste0(cell_prefix, ":y_X_D1")]],
        D_X    = splits[[paste0(cell_prefix, ":D_X")]],
        D      = splits[[paste0(cell_prefix, ":D")]])
    } else {
      list(D_X = list(subsamples = cell_subsamples))
    }#IFELSE

    # Run ddml_att for this cell
    fit <- ddml_att(
      y = delta_y, D = D_cell, X = X_cell,
      learners = learners, learners_DX = learners_qX,
      sample_folds = length(cell_subsamples),
      ensemble_type = ensemble_type,
      shortstack = shortstack, cv_folds = cv_folds,
      custom_ensemble_weights = custom_ensemble_weights,
      custom_ensemble_weights_DX = custom_ensemble_weights_qX,
      cluster_variable = cluster_variable[keep],
      trim = trim, silent = TRUE, parallel = parallel,
      save_crossval = save_crossval,
      fitted = cell_fitted,
      splits = cell_splits,
      ...)

    # Collect per-cell diagnostics (flat keying)
    for (eq in names(fit$ensemble_weights)) {
      key <- paste0(cell_prefix, ":", eq)
      all_ensemble_weights[[key]] <- fit$ensemble_weights[[eq]]
      all_mspe[[key]] <- fit$mspe[[eq]]
      all_r2[[key]] <- fit$r2[[eq]]
    }#FOR
    for (eq in names(fit$fitted)) {
      key <- paste0(cell_prefix, ":", eq)
      all_fitted[[key]] <- fit$fitted[[eq]]
    }#FOR
    for (eq in names(fit$splits)) {
      key <- paste0(cell_prefix, ":", eq)
      all_splits[[key]] <- fit$splits[[eq]]
    }#FOR

    # Determine nensb from first fit and allocate
    if (is.null(nensb)) {
      nensb <- ncol(fit$coefficients)
      ens_type <- colnames(fit$coefficients)
      if (is.null(ens_type)) ens_type <- fit$ensemble_type
      coef_mat <- matrix(NA_real_, C, nensb)
      scores_arr <- array(0, dim = c(n, C, nensb))
      J_arr <- array(0, dim = c(C, C, nensb))
      inf_func_full <- array(0, dim = c(n, C, nensb))
      dinf_dtheta_full <- array(0, dim = c(n, C, C, nensb))
    }#IF

    # Embed cell scores into full-sample arrays
    for (j in seq_len(nensb)) {
      att_val <- fit$coefficients[1, j]
      coef_mat[idx, j] <- att_val

      # Cell-level components (p = 1 for ATT)
      cell_scores <- fit$scores[, 1, j]
      cell_J <- fit$J[1, 1, j]
      cell_inf_func <- fit$inf_func[, 1, j]
      cell_dinf_dtheta <- fit$dinf_dtheta[, 1, 1, j]

      # Rescale from cell-level to population-level:
      # J_pop = (n_cell/n) * J_cell, so IF_pop = scores/J_pop
      #       = (n/n_cell) * scores/J_cell = (n/n_cell) * IF_cell
      pop_scale <- n / n_cell

      scores_arr[keep, idx, j] <- cell_scores
      inf_func_full[keep, idx, j] <- cell_inf_func * pop_scale
      # Jacobian: rescale to full-sample average
      J_arr[idx, idx, j] <- (n_cell / n) * cell_J
      # dinf_dtheta: same rescaling (used for HC3 leverage)
      dinf_dtheta_full[keep, idx, idx, j] <- cell_dinf_dtheta * pop_scale
    }#FOR
  }#FOR

  # Coefficient names ----------------------------------------------------------

  coef_names <- paste0("ATT(", cell_info$group, ",", cell_info$time, ")")
  rownames(coef_mat) <- coef_names
  colnames(coef_mat) <- ens_type

  # Store global subsamples for pass-through round-trip
  all_splits[["global"]] <- list(subsamples = global_subsamples)

  # Output ---------------------------------------------------------------------

  announce_finish(t0, messages, silent)

  ddml(
    coefficients = coef_mat,
    ensemble_weights = all_ensemble_weights,
    mspe = all_mspe,
    r2 = all_r2,
    scores = scores_arr,
    J = J_arr,
    inf_func = inf_func_full,
    dinf_dtheta = dinf_dtheta_full,
    nobs = n,
    coef_names = coef_names,
    estimator_name = "Group-Time Average Treatment Effects on the Treated",
    ensemble_type = ens_type,
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    fitted = all_fitted,
    splits = all_splits,
    call = cl,
    subclass = "ddml_attgt",
    # ddml_attgt-specific fields
    learners = learners,
    learners_qX = learners_qX,
    cell_info = cell_info,
    G = G,
    control_group = control_group,
    anticipation = anticipation
  )
}#DDML_ATTGT
