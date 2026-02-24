compute_ddml_variance <- function(scores, J,
                                 cluster_variable = NULL,
                                 type = "HC1") {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  sc <- as.matrix(scores)
  J_j <- as.matrix(J)
  p <- ncol(sc)

  # Cluster aggregation: rowsum scores to cluster level
  if (!is.null(cluster_variable) &&
      !identical(seq_len(nrow(sc)),
                 as.integer(cluster_variable))) {
    sc <- rowsum(sc, cluster_variable)
  }#IF

  n_eff <- nrow(sc)

  if (type == "HC3") {
    PtP_inv <- solve(crossprod(sc))
    h <- rowSums((sc %*% PtP_inv) * sc)
    sc <- sc / (1 - h)
    meat <- crossprod(sc) / n_eff
    J_inv <- solve(J_j)
    V <- J_inv %*% meat %*% t(J_inv) / n_eff
  } else if (type == "HC1") {
    meat <- crossprod(sc) / n_eff
    J_inv <- solve(J_j)
    V <- J_inv %*% meat %*% t(J_inv) *
      n_eff / (n_eff - p) / n_eff
  } else {
    # HC0
    meat <- crossprod(sc) / n_eff
    J_inv <- solve(J_j)
    V <- J_inv %*% meat %*% t(J_inv) / n_eff
  }#IFELSE

  V
}#COMPUTE_DDML_VARIANCE

compute_ddml_inference <- function(coefficients, scores, J,
                                  coef_names = NULL,
                                  ensemble_type,
                                  cluster_variable = NULL,
                                  type = "HC1") {
  nensb <- length(ensemble_type)
  p <- NCOL(scores[[1]])

  inf_results <- array(0, dim = c(p, 4, nensb))

  for (j in seq_len(nensb)) {
    if (is.matrix(coefficients)) {
      theta_j <- coefficients[, j]
    } else {
      theta_j <- coefficients[j]
    }#IFELSE

    V <- compute_ddml_variance(scores[[j]], J[[j]],
                              cluster_variable,
                              type = type)
    se <- sqrt(diag(V))
    t_val <- theta_j / se
    p_val <- 2 * stats::pnorm(abs(t_val),
                               lower.tail = FALSE)

    inf_results[, 1, j] <- theta_j
    inf_results[, 2, j] <- se
    inf_results[, 3, j] <- t_val
    inf_results[, 4, j] <- p_val
  }#FOR

  dimnames(inf_results) <- list(
    coef_names,
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    ensemble_type
  )
  inf_results
}#COMPUTE_DDML_INFERENCE
