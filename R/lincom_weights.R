#' Difference-in-Differences Aggregation Weights for lincom
#'
#' @description Constructs the contrast matrix \eqn{R} and
#'     its influence function matrix \code{inf_func_R} for
#'     standard DiD aggregation types. The output is
#'     designed to be passed directly to
#'     \code{\link{lincom}}.
#'
#' @details Let \eqn{\theta^{(g,t)}_0} denote the GT-ATT
#'     from \code{\link{ddml_attgt}}. Each aggregation type
#'     defines a summary parameter as a weighted average of
#'     GT-ATTs over a subset of post-treatment cells
#'     (\eqn{t \geq g}).
#'
#' \strong{Dynamic} (\code{type = "dynamic"}): aggregates
#' by event time \eqn{e = t - g} (Callaway and Sant'Anna,
#' 2021, eq. 9). For each \eqn{e}:
#'
#' \deqn{\tau_0(e) = \sum_{g \in \mathcal{G}}
#'   \mathbf{1}\{g + e \leq T\}\,
#'   \Pr(G = g \mid G + e \leq T)\,
#'   \theta^{(g,\, g+e)}_0.}
#'
#' \strong{Group} (\code{type = "group"}): aggregates by
#' cohort \eqn{g}. For each \eqn{g}:
#'
#' \deqn{\theta(g) = \frac{1}{|\mathcal{T}_g|}
#'   \sum_{t \in \mathcal{T}_g}
#'   \theta^{(g,t)}_0,}
#'
#' where \eqn{\mathcal{T}_g = \{t : t \geq g\}} and the
#' weights reduce to uniform within each cohort.
#'
#' \strong{Calendar} (\code{type = "calendar"}): aggregates
#' by time period \eqn{t}. For each \eqn{t}:
#'
#' \deqn{\theta(t) = \sum_{g:\, g \leq t}
#'   \frac{P(G = g)}{\sum_{g':\, g' \leq t} P(G = g')}\,
#'   \theta^{(g,t)}_0.}
#'
#' \strong{Simple} (\code{type = "simple"}): a single
#' weighted average across all post-treatment cells:
#'
#' \deqn{\theta_{ATT} = \sum_{(g,t):\, t \geq g}
#'   \frac{P(G = g)}{\sum_{(g',t'):\, t' \geq g'}
#'   P(G = g')}\, \theta^{(g,t)}_0.}
#'
#' The influence function for the estimated weights is
#' derived via the quotient rule and passed to
#' \code{\link{lincom}} as \code{inf_func_R}.
#'
#' @param fit A \code{ddml_attgt} or \code{ddml_rep} object
#'     whose underlying fits are \code{ddml_attgt}.
#' @param type Aggregation type: \code{"dynamic"} (default),
#'     \code{"group"}, \code{"simple"}, or \code{"calendar"}.
#' @param min_e,max_e Event-time range filter (dynamic only).
#'     Cells with event time outside \code{[min_e, max_e]}
#'     are excluded.
#' @param fit_idx Integer index of the fit (ensemble type)
#'     to use for computing the weighting leverage, or
#'     \code{NULL} (default) for all ensemble types.
#'     When \code{NULL}, \code{dinf_dR} is a 4D array.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{R}}{A \eqn{(C \times q)}{C x q} contrast
#'       matrix where \eqn{C} is the number of GT cells.}
#'   \item{\code{inf_func_R}}{An \eqn{(n \times C \times q)}
#'       {n x C x q} array of influence functions for
#'       the contrast matrix \eqn{R}. Slice \code{[,,k]}
#'       contains the IFs for column \eqn{k} of \eqn{R}.}
#'   \item{\code{dinf_dR}}{An \eqn{(n \times q \times q)}
#'       {n x q x q} array of weighting leverage. Each
#'       slice is the constant matrix \eqn{V'V} where
#'       \eqn{V_{g,k} = \sum_{c \in \mathcal{K}_k,\,
#'       g_c = g} (\theta_c - \gamma_k) / S_k}. See
#'       D89 §5.3.}
#'   \item{\code{labels}}{Character vector of length
#'       \eqn{q} naming the aggregated quantities.}
#' }
#'
#' @seealso \code{\link{lincom}}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' n <- 200; T_ <- 4
#' X <- matrix(rnorm(n * 2), n, 2)
#' G <- sample(c(3, 4, Inf), n, replace = TRUE,
#'             prob = c(0.3, 0.3, 0.4))
#' y <- matrix(rnorm(n * T_), n, T_)
#' fit <- ddml_attgt(y, X, t = 1:T_, G = G,
#'                 learners = list(what = ols),
#'                 sample_folds = 2,
#'                 silent = TRUE)
#' w <- lincom_weights_did(fit, type = "dynamic")
#' dyn <- lincom(fit, R = w$R,
#'               inf_func_R = w$inf_func_R,
#'               dinf_dR = w$dinf_dR,
#'               labels = w$labels)
#' summary(dyn)
#' }
#'
#' @references
#' Callaway B, Sant'Anna P H C (2021). "Difference-in-Differences
#' with multiple time periods." Journal of Econometrics,
#' 225(2), 200-230.
#'
#' @export
lincom_weights_did <- function(fit,
                         type = c("dynamic", "group",
                                  "simple", "calendar"),
                         min_e = -Inf, max_e = Inf,
                         fit_idx = NULL) {
  type <- match.arg(type)
  ref <- if (inherits(fit, "ddml_rep")) {
    fit[[1]]
  } else {
    fit
  }#IFELSE
  if (!inherits(ref, "ddml_attgt")) {
    stop("'fit' must be a ddml_attgt or ddml_rep object.", call. = FALSE)
  }#IF

  cell_info <- ref$cell_info
  group <- cell_info$group
  time  <- cell_info$time
  C     <- length(group)
  G     <- ref$G
  n     <- ref$nobs

  # Group probabilities
  glist <- sort(unique(group))
  pg_by_group <- sapply(glist, function(g) mean(G == g))
  pg <- pg_by_group[match(group, glist)]

  # Post-treatment cells
  post <- which(group <= time)

  # Build aggregation groups
  agg <- switch(type,
    dynamic = {
      et <- time - group
      all_et <- sort(unique(et))
      all_et <- all_et[all_et >= min_e & all_et <= max_e]
      lapply(all_et, function(e) {
        list(keepers = which(et == e), label = paste0("e=", e))
      })
    },
    group = {
      lapply(glist, function(g) {
        list(keepers = which(group == g & seq_len(C) %in% post),
             label = paste0("g=", g))
      })
    },
    calendar = {
      tlist <- sort(unique(time[post]))
      lapply(tlist, function(t1) {
        list(keepers = which(time == t1 & seq_len(C) %in% post),
             label = paste0("t=", t1))
      })
    },
    simple = list(list(keepers = post, label = "ATT"))
  )

  q <- length(agg)
  labels <- vapply(agg, `[[`, character(1), "label")
  R <- matrix(0, C, q)
  inf_func_R <- array(0, dim = c(n, C, q))

  for (k in seq_len(q)) {
    keepers <- agg[[k]]$keepers
    if (length(keepers) == 0L) next  # reference period
    w_norm <- pg[keepers] / sum(pg[keepers])
    R[keepers, k] <- w_norm

    # Weight IF via quotient rule:
    #   w_c = pg_c / S,  S = sum(pg[keepers])
    #   phi_i^{w_c} = (1{G==g_c} - pg_c) / S
    #           - w_c * sum_{c'}(1{G==g_{c'}} - pg_{c'}) / S
    S <- sum(pg[keepers])
    if1 <- sapply(keepers, function(kk) {
      (as.numeric(G == group[kk]) - pg[kk]) / S
    })
    if2 <- rowSums(if1) %*% t(w_norm)
    wif_k <- if1 - if2
    for (j in seq_along(keepers)) {
      inf_func_R[, keepers[j], k] <- wif_k[, j]
    }#FOR
  }#FOR

  # Weighting leverage: dinf_dR = V'V (constant across i) ----
  # V_{g,k} = sum_{c in K_k, gc=g} (theta_c - gamma_k) / S_k
  nensb <- ncol(ref$coefficients)
  nG <- length(glist)
  if (is.null(fit_idx)) {
    j_seq <- seq_len(nensb)
  } else {
    j_seq <- fit_idx
  }#IFELSE

  # Helper: compute VtV for a single ensemble column
  compute_VtV <- function(jj) {
    theta <- ref$coefficients[, jj]
    gamma <- as.numeric(crossprod(R, theta))
    V <- matrix(0, nG, q)
    for (k in seq_len(q)) {
      keepers <- agg[[k]]$keepers
      if (length(keepers) == 0L) next
      S_k <- sum(pg[keepers])
      for (gg in seq_len(nG)) {
        cells_gk <- keepers[group[keepers] == glist[gg]]
        if (length(cells_gk) > 0L) {
          V[gg, k] <- sum(theta[cells_gk] - gamma[k]) / S_k
        }#IF
      }#FOR
    }#FOR
    crossprod(V)
  }#COMPUTE_VTV

  if (length(j_seq) == 1L) {
    # Single ensemble -> 3D output (backward compatible)
    VtV <- compute_VtV(j_seq)
    dinf_dR <- array(rep(VtV, each = n), dim = c(n, q, q))
  } else {
    # Multi-ensemble -> 4D output
    dinf_dR <- array(0, dim = c(n, q, q, length(j_seq)))
    for (jj in seq_along(j_seq)) {
      VtV <- compute_VtV(j_seq[jj])
      dinf_dR[, , , jj] <- array(rep(VtV, each = n),
                                  dim = c(n, q, q))
    }#FOR
  }#IFELSE

  colnames(R) <- labels
  list(R = R, inf_func_R = inf_func_R,
       dinf_dR = dinf_dR, labels = labels)
}#LINCOM_WEIGHTS_DID
