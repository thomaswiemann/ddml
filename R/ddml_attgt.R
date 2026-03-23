#' Estimator for Group-Time Average Treatment Effects
#'
#' @family ddml estimators
#'
#' @description Estimator for group-time average treatment effects
#'     on the treated (GT-ATT) in staggered Difference-in-Differences
#'     designs.
#'
#' @details
#' \strong{Parameter of Interest:} \code{ddml_attgt} provides a
#'     Double/Debiased Machine Learning estimator for the group-time
#'     average treatment effects on the treated (GT-ATT) in the
#'     staggered adoption model. For each group \eqn{g} and time
#'     period \eqn{t}, define the differenced outcome
#'     \eqn{\Delta_g Y_{i,t} = Y_{i,t} - Y_{i,g^*}} where
#'     \eqn{g^*} is the universal base period. The GT-ATT is:
#'
#' \deqn{\theta_0^{(g,t)} = E[\Delta_g Y_{i,t} | G_i = g]
#'     - E[E[\Delta_g Y_{i,t} | X_i, G_i \ne g, G_i > t] | G_i = g]}
#'
#' \strong{Neyman Orthogonal Score:} The Neyman orthogonal score
#'     is:
#'
#' \deqn{m_i^{(g,t)} =
#'     \frac{\mathbf{1}\{G_i = g\} (\Delta_g Y_{i,t}
#'     - \ell^{(g,t)}(X_i))}{\pi^g}
#'     - \frac{q^{(g,t)}(X_i) \mathbf{1}\{G_i \ne g\}
#'     \mathbf{1}\{G_i > t\} (\Delta_g Y_{i,t}
#'     - \ell^{(g,t)}(X_i))}{\pi^g (1 - q^{(g,t)}(X_i))}
#'     - \frac{\mathbf{1}\{G_i = g\}}{\pi^g} \theta}
#'
#' where the nuisance parameters are
#'     \eqn{\eta = (\ell, q, \pi)} taking true values
#'     \eqn{\ell_0^{(g,t)}(X) = E[\Delta_g Y_{i,t} \mid
#'     G_i \ne g, G_i > t, X_i]},
#'     \eqn{q_0^{(g,t)}(X) = \Pr(G_i = g \mid X_i,
#'     \{G_i = g\} \cup \{G_i > t\})},
#'     and \eqn{\pi_0^g = \Pr(G_i = g)}.
#'
#' \strong{Jacobian:}
#'
#' \deqn{J^{(g,t)} = -1}
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
  global_cv_subsamples <- global_indxs$cv_subsamples

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

  # Reduced-form estimation ----------------------------------------------------

  # Cross-fit global pi^g = Pr(G_i = g) per group
  pi_g_fitted <- list()
  for (g_val in groups) {
    g_key <- as.character(g_val)
    D_g <- as.integer(G == g_val)
    pi_fitted_init <- if (!is.null(fitted)) fitted[[paste0("pi_g:", g_key)]]
    pi_g_res <- get_CEF(D_g, matrix(1, n, 1),
                        learners = list(what = ols, args = list(const = FALSE)),
                        ensemble_type = "average",
                        shortstack = FALSE,
                        subsamples = global_subsamples,
                        cv_subsamples = global_cv_subsamples,
                        silent = TRUE,
                        label = paste0("pi(G=", g_key, ")"),
                        fitted = pi_fitted_init)
    pi_g_fitted[[g_key]] <- pi_g_res
  }#FOR

  # Pre-allocate output arrays
  nensb <- NULL
  coef <- NULL
  scores <- NULL
  J <- NULL

  ensemble_weights <- list()
  mspe <- list()
  r2 <- list()
  fitted_list <- list()
  splits_list <- list()
  cell_info <- data.frame(
    group = integer(C), time = integer(C),
    base_period = integer(C),
    n_treated = integer(C), n_control = integer(C))

  # Store global pi^g fitted entries for pass-through
  for (g_key in names(pi_g_fitted)) {
    fitted_list[[paste0("pi_g:", g_key)]] <-
      build_fitted_entry(pi_g_fitted[[g_key]], save_crossval)
  }#FOR
  splits_list[["global"]] <- list(subsamples = global_subsamples)

  # Per-cell reduced-form estimation and score construction
  for (idx in seq_len(C)) {
    gtp <- gt_list[[idx]]
    g_val <- gtp$g; tt <- gtp$t
    cell_prefix <- paste0("ATT(", g_val, ",", tt, ")")
    g_key <- as.character(g_val)

    # Cell membership ---------------------------------------------------------

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

    # Cell-level reduced forms via ddml_att ------------------------------------

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
      ensemble_weights[[key]] <- fit$ensemble_weights[[eq]]
      mspe[[key]] <- fit$mspe[[eq]]
      r2[[key]] <- fit$r2[[eq]]
    }#FOR
    for (eq in names(fit$fitted)) {
      key <- paste0(cell_prefix, ":", eq)
      fitted_list[[key]] <- fit$fitted[[eq]]
    }#FOR
    for (eq in names(fit$splits)) {
      key <- paste0(cell_prefix, ":", eq)
      splits_list[[key]] <- fit$splits[[eq]]
    }#FOR

    # Determine nensb from first fit and allocate
    if (is.null(nensb)) {
      nensb <- ncol(fit$coefficients)
      ens_type <- colnames(fit$coefficients)
      if (is.null(ens_type)) ens_type <- fit$ensemble_type
      coef <- matrix(NA_real_, C, nensb)
      psi_a_arr <- array(0, dim = c(n, C))
      psi_b_arr <- array(0, dim = c(n, C, nensb))
      scores <- array(0, dim = c(n, C, nensb))
      J <- array(0, dim = c(C, C, nensb))
      inf_func <- array(0, dim = c(n, C, nensb))
      dinf_dtheta <- array(0, dim = c(n, C, C, nensb))
    }#IF

    # Score construction -------------------------------------------------------
    # Extract nuisance estimates from ddml_att fitted entry and construct
    # population-level scores.

    # E[DeltaY | D=0, X] extrapolated to all cell units
    y_X_D0_entry <- fit$fitted$y_X_D0
    g_X_D0 <- matrix(NA_real_, n_cell, nensb)
    g_X_D0[D_cell == 0, ] <- as.matrix(y_X_D0_entry$cf_fitted)
    for (k in seq_along(cell_subsamples)) {
      fold_k <- cell_subsamples[[k]]
      d1_in_k <- fold_k[D_cell[fold_k] == 1]
      if (length(d1_in_k) > 0) {
        g_X_D0[d1_in_k, ] <- as.matrix(
          y_X_D0_entry$auxiliary_fitted[[k]])
      }#IF
    }#FOR

    # E[D|X] trimmed propensity
    m_X <- as.matrix(fit$fitted$D_X$cf_fitted)
    m_X_tr <- trim_propensity_scores(m_X, trim, ens_type)

    # Global pi^g = Pr(G_i = g)
    pi_g_cell <- pi_g_fitted[[g_key]]$cf_fitted[keep, 1]

    # Score components
    D_cell_mat <- matrix(D_cell, n_cell, nensb)
    delta_y_mat <- matrix(delta_y, n_cell, nensb)
    pi_g_mat <- matrix(pi_g_cell, n_cell, nensb)

    psi_b_arr[keep, idx, ] <- D_cell_mat *
      (delta_y_mat - g_X_D0) / pi_g_mat -
      m_X_tr * (1 - D_cell_mat) * (delta_y_mat - g_X_D0) /
      (pi_g_mat * (1 - m_X_tr))
    psi_a_arr[keep, idx] <- -D_cell / pi_g_cell
  }#FOR

  # Target parameter & influence function --------------------------------------

  mean_psi_a <- colMeans(psi_a_arr)     # C-vector (= J diagonal)
  J_inv_vec <- 1 / mean_psi_a           # C-vector

  # J and dinf_dtheta are ensemble-independent, they only depend on constants...
  J_diag_idx <- cbind(seq_len(C), seq_len(C))
  dinf_dtheta_common <- t(t(psi_a_arr) * J_inv_vec)
  for (j in seq_len(nensb)) {
    J[, , j][J_diag_idx] <- mean_psi_a
    dinf_dtheta[, , , j] <- dinf_dtheta_common
  }#FOR

  for (j in seq_len(nensb)) {
    coef[, j] <- -colMeans(psi_b_arr[, , j]) / mean_psi_a

    scores[, , j] <- t(t(psi_a_arr) * coef[, j]) + psi_b_arr[, , j]
    inf_func[, , j] <- t(t(scores[, , j]) * J_inv_vec)
  }#FOR

  coef_names <- paste0("ATT(", cell_info$group, ",", cell_info$time, ")")
  rownames(coef) <- coef_names
  colnames(coef) <- ens_type

  # Output ---------------------------------------------------------------------

  announce_finish(t0, messages, silent)

  ddml(
    coefficients = coef,
    ensemble_weights = ensemble_weights,
    mspe = mspe,
    r2 = r2,
    scores = scores,
    J = J,
    inf_func = inf_func,
    dinf_dtheta = dinf_dtheta,
    nobs = n,
    coef_names = coef_names,
    estimator_name = "Group-Time Average Treatment Effects on the Treated",
    ensemble_type = ens_type,
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    fitted = fitted_list,
    splits = splits_list,
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
