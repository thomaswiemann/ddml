# Pairwise cross-validation comparison test.
#
# Compares predictive ability of two learners using their
# OOS residuals. Uses multiplier bootstrap for inference,
# with fold-level demeaning to account for cross-fitting
# dependence.
#
# @param resid1 Numeric vector of OOS residuals from learner 1.
# @param resid2 Numeric vector of OOS residuals from learner 2.
# @param fid Integer vector of fold IDs (same length as resid1).
# @param bootnum Number of bootstrap replications.
#
# @return P-value for H0: learner 1 has equal or larger MSPE
#     than learner 2. Small p-value means learner 1 is better.
cvc_pairwise <- function(resid1, resid2, fid,
                         bootnum = 500) {
  n <- length(resid1)
  zeta <- resid1^2 - resid2^2

  # Demean by fold to remove cross-fitting dependence
  fid_unique <- unique(fid)
  zeta_m <- numeric(n)
  for (k in fid_unique) {
    sel <- which(fid == k)
    zeta_m[sel] <- mean(zeta[sel])
  }#FOR
  zeta_til <- zeta - zeta_m
  zeta_sd <- sqrt(stats::var(zeta_til))
  if (zeta_sd < .Machine$double.eps) return(NA_real_)

  Tx <- sqrt(n) * mean(zeta) / zeta_sd

  # Multiplier bootstrap (no model refitting)
  zeta_scaled <- zeta_til / zeta_sd
  Txb <- vapply(seq_len(bootnum), function(b) {
    (1 / sqrt(n)) * sum(zeta_scaled * stats::rnorm(n))
  }, numeric(1))

  mean(Txb > Tx)
}#CVC_PAIRWISE

# One-vs-many cross-validation comparison test.
#
# Tests whether learner i is dominated by any learner in the
# comparison set. Uses a sup-type statistic with multiplier
# bootstrap.
#
# @param resid_base Numeric vector of OOS residuals from the
#     base learner being tested.
# @param resid_others Matrix (n x K) of OOS residuals from
#     the comparison learners.
# @param fid Integer vector of fold IDs.
# @param bootnum Number of bootstrap replications.
#
# @return P-value. Large p-value means the base learner is
#     not significantly worse than the best alternative.
cvc_one_vs_many <- function(resid_base, resid_others, fid,
                            bootnum = 500) {
  n <- length(resid_base)
  K <- ncol(resid_others)
  resid_base_mat <- matrix(resid_base, n, K)

  zeta <- resid_base_mat^2 - resid_others^2

  # Demean by fold
  fid_unique <- unique(fid)
  zeta_til <- zeta
  for (k in fid_unique) {
    sel <- which(fid == k)
    fold_means <- colMeans(zeta[sel, , drop = FALSE])
    zeta_til[sel, ] <- sweep(
      zeta[sel, , drop = FALSE], 2, fold_means)
  }#FOR

  zeta_m <- colMeans(zeta)
  zeta_sd <- apply(zeta_til, 2, stats::sd)
  zeta_sd[zeta_sd < .Machine$double.eps] <- Inf

  # Sup-type test statistic
  Tx <- max(sqrt(n) * zeta_m / zeta_sd)

  # Multiplier bootstrap
  zeta_scaled <- sweep(zeta_til, 2, zeta_sd, FUN = "/")
  Txb <- vapply(seq_len(bootnum), function(b) {
    bw <- stats::rnorm(n)
    max(colSums(zeta_scaled * bw) / sqrt(n))
  }, numeric(1))

  mean(Txb > Tx)
}#CVC_ONE_VS_MANY

# Model confidence set via cross-validation comparison.
#
# For each learner, tests whether it is dominated by any
# other learner. Learners with large p-values form the
# model confidence set.
#
# @param resid_mat Matrix (n x nlearners) of OOS residuals.
# @param fid Integer vector of fold IDs.
# @param bootnum Number of bootstrap replications.
# @param alpha Significance level for the confidence set.
#
# @return List with: pvalues (numeric vector of p-values per
#     learner), confidence_set (logical vector).
cvc_confidence_set <- function(resid_mat, fid,
                               bootnum = 500,
                               alpha = 0.05) {
  nlearners <- ncol(resid_mat)
  pvalues <- numeric(nlearners)
  for (i in seq_len(nlearners)) {
    pvalues[i] <- cvc_one_vs_many(
      resid_mat[, i], resid_mat[, -i, drop = FALSE],
      fid, bootnum)
  }#FOR
  list(pvalues = pvalues,
       confidence_set = pvalues > alpha)
}#CVC_CONFIDENCE_SET

# Main CVC entry point.
#
# Runs pairwise comparisons and the model confidence set
# procedure on a matrix of OOS residuals.
#
# @param cv_resid Matrix (n x nlearners) of OOS residuals.
# @param subsamples List of sample fold index vectors (used
#     to derive fold IDs).
# @param bootnum Number of bootstrap replications.
# @param alpha Significance level for the confidence set.
#
# @return List with: pairwise (nlearners x nlearners matrix
#     of p-values), confidence_set (logical vector),
#     pvalues (numeric vector from confidence set procedure).
cvc_test <- function(cv_resid, subsamples,
                     bootnum = 500, alpha = 0.05) {
  n <- nrow(cv_resid)
  nlearners <- ncol(cv_resid)

  # Derive fold IDs from subsamples
  fid <- integer(n)
  for (k in seq_along(subsamples)) {
    fid[subsamples[[k]]] <- k
  }#FOR

  # Pairwise comparisons
  pairwise <- matrix(NA_real_, nlearners, nlearners)
  if (!is.null(colnames(cv_resid))) {
    rownames(pairwise) <- colnames(pairwise) <-
      colnames(cv_resid)
  }#IF
  for (i in seq_len(nlearners)) {
    for (j in seq_len(nlearners)) {
      if (i != j) {
        pairwise[i, j] <- cvc_pairwise(
          cv_resid[, i], cv_resid[, j],
          fid, bootnum)
      }#IF
    }#FOR
  }#FOR

  # Model confidence set
  cs <- cvc_confidence_set(cv_resid, fid, bootnum, alpha)

  list(pairwise = pairwise,
       pvalues = cs$pvalues,
       confidence_set = cs$confidence_set)
}#CVC_TEST
